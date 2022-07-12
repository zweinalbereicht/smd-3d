use rand::distributions::Distribution;

use rand_distr::Uniform;
use rand_distr::{Exp, Normal};
use rand_pcg::*;

use num_complex::{Complex, Complex64};
use std::collections::HashMap;

use std::f64::consts::PI;
use std::fmt;
use std::process::exit;

use crate::ejectionenvironement::*;
use crate::tools::*;

pub struct EjectionParticle {
    pub radial_position: f64,
    pub angular_position: f64,
    pub lifetime: f64,
    pub status: Status,
}

impl fmt::Display for EjectionParticle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "radial position : {}\nangular position : {}\nlifetime : {}\nstatus : {}",
            self.radial_position,
            self.angular_position,
            self.lifetime,
            match self.status {
                Status::Dead(_) => "dead",
                Status::Bulk(_) => "bulk",
                Status::Absorbed => "absorbec",
            }
        )
    }
}

impl Default for EjectionParticle {
    fn default() -> Self {
        EjectionParticle {
            radial_position: 1.,
            angular_position: 0.,
            lifetime: 0.,
            status: Status::Bulk(PointlikeTarket::NotCrossed),
        }
    }
}

impl EjectionParticle {
    pub fn reset(&mut self, radial_pos: f64, angular_pos: f64) {
        self.radial_position = radial_pos;
        self.angular_position = angular_pos;
        self.lifetime = 0.;
        self.status = Status::Bulk(PointlikeTarket::NotCrossed);
    }
}

// we will put it as a class function for portability
impl EjectionParticle {
    // helper funciton to go to complex
    pub fn position_to_complex(&self) -> Complex<f64> {
        Complex::from_polar(self.radial_position, self.angular_position)
    }

    pub fn move_particle(
        &mut self,
        environement: &EjectionEnvironment,
        rng: &mut Lcg128Xsl64,
        dt: f64,
    ) {
        //println!("{}\n", particle);
        match self.status {
            Status::Dead(_) => {
                //println!("Particle is dead after time {}", particle.lifetime);
                exit(0);
            }
            Status::Absorbed => {
                //println!("in the absorbed phase");
                let exp = Exp::new(1. / environement.desorption_time).unwrap();
                let normal = Normal::new(0., 1.).unwrap();
                let elapsed_time = exp.sample(rng);
                let move_along_boundary = (2. * environement.boundary_coefficient * elapsed_time)
                    .sqrt()
                    * normal.sample(rng)
                    / environement.outer_radius;

                // check if we crossed zero during our move - either we changed sign or went around
                let new_pos = self.angular_position + move_along_boundary;
                if (new_pos.signum() == self.angular_position.signum()) | (new_pos.abs() > 2. * PI)
                {
                    self.status = Status::Bulk(PointlikeTarket::Crossed);
                } else {
                    self.status = Status::Bulk(PointlikeTarket::NotCrossed);
                }

                //move
                // this is a bit weird but its to ensure we stay in the right angle zone
                self.angular_position =
                    Complex::from_polar(environement.outer_radius, new_pos).arg();

                self.lifetime += elapsed_time;
                //eject
                self.radial_position -= environement.ejection_length;
            }

            Status::Bulk(_) => {
                //println!("in the bulk");

                // convert to complex for easier use
                let mut particle_position_complex = self.position_to_complex();

                //implement the jump
                let normal = Normal::new(0., 1.).unwrap();
                let increment = Complex::new(
                    (2. * environement.bulk_coefficient * dt).sqrt() * normal.sample(rng),
                    (2. * environement.bulk_coefficient * dt).sqrt() * normal.sample(rng),
                );

                //move the particle
                particle_position_complex += increment;
                self.lifetime += dt;

                // if we jump outside we get absorbed
                if particle_position_complex.norm() >= environement.outer_radius {
                    self.radial_position = environement.outer_radius;
                    self.angular_position = particle_position_complex.arg();
                    self.status = Status::Absorbed;
                } else {
                    self.angular_position = particle_position_complex.arg();
                    self.radial_position = particle_position_complex.norm();
                }

                // if we are inside a trap, we need to stop the process. Status::Dead stores the
                // index of the absorbing trap
                if let Some(absorbing_trap) = &environement
                    .traps
                    .iter()
                    .position(|&x| single_trap_intersect(particle_position_complex, x))
                {
                    self.status = Status::Dead(*absorbing_trap);
                }
            }
        }
    }

    // we implement the event driven move algorithm where we just pick a random point on the
    // exterior of the biggest possible circle. It is also useful to calculate quicker
    // mfpts.
    pub fn move_particle_event(
        &mut self,
        environement: &EjectionEnvironment,
        rng: &mut Lcg128Xsl64,
        tolerance: f64,
    ) {
        //println!("{}\n", particle);
        match self.status {
            Status::Dead(_) => {
                //println!("Particle is dead after time {}", particle.lifetime);
                exit(0);
            }
            Status::Absorbed => {
                //println!("in the absorbed phase");
                let exp = Exp::new(1. / environement.desorption_time).unwrap();
                let normal = Normal::new(0., 1.).unwrap();
                let elapsed_time = exp.sample(rng);
                let move_along_boundary = (2. * environement.boundary_coefficient * elapsed_time)
                    .sqrt()
                    * normal.sample(rng)
                    / environement.outer_radius;

                // check if we crossed zero during our move - either we changed sign or went around
                let new_pos = self.angular_position + move_along_boundary;
                if (new_pos.signum() == self.angular_position.signum()) | (new_pos.abs() > 2. * PI)
                {
                    self.status = Status::Bulk(PointlikeTarket::Crossed);
                } else {
                    self.status = Status::Bulk(PointlikeTarket::NotCrossed);
                }

                //move
                // this is a bit weird but its to ensure we stay in the right angle zone
                self.angular_position =
                    Complex::from_polar(environement.outer_radius, new_pos).arg();

                // here we only add the mean time, should be sufficient.
                self.lifetime += environement.desorption_time;

                //eject
                self.radial_position -= environement.ejection_length;
            }
            // in here we are going to implement the algortithm that moves the particle to the boundary of a decentered circle, that will intersect
            // both the outer rim and the inner circle. This is based on some simulation algorithms devised by other people, see overleaf for more
            // info
            Status::Bulk(_) => {
                //println!("in the bulk");

                // convert to complex for easier use
                let mut particle_position_complex = self.position_to_complex();

                // find smallest radius for the ball between bonudary and traps, and add a little tolerance.
                let smallest_radius = (environement.outer_radius - self.radial_position).min(
                    environement
                        .traps
                        .iter()
                        .map(|t| (particle_position_complex - t.center).norm() - t.radius)
                        .reduce(f64::min)
                        .unwrap(),
                ) + tolerance;

                //calculate the jump
                let increment =
                    Complex::from_polar(smallest_radius, Uniform::new(-1., 1.).sample(rng) * PI);

                //move the particle
                particle_position_complex += increment;

                //increment time by the mean time needed to exit the biggest circle possible -->
                //this value can be found in getoor and blumenthal.
                self.lifetime += smallest_radius.powi(2) / (4. * environement.bulk_coefficient);

                // if we jump outside we get absorbed
                if particle_position_complex.norm() >= environement.outer_radius {
                    self.radial_position = environement.outer_radius;
                    self.angular_position = particle_position_complex.arg();
                    self.status = Status::Absorbed;
                } else {
                    self.angular_position = particle_position_complex.arg();
                    self.radial_position = particle_position_complex.norm();
                }

                // if we are inside a trap, we need to stop the process. Status::Dead stores the
                // index of the absorbing trap
                if let Some(absorbing_trap) = &environement
                    .traps
                    .iter()
                    .position(|&x| (particle_position_complex - x.center).norm() - x.radius < 0.)
                {
                    self.status = Status::Dead(*absorbing_trap);
                }
            }
        }
    }

    pub fn move_and_record(
        &mut self,
        environement: &EjectionEnvironment,
        rng: &mut Lcg128Xsl64,
        dt: f64,
        particle_size: f64,
        total_territory: &mut HashMap<CartesianGridLocation, Complex64>,
    ) {
        // first move the particle
        self.move_particle(environement, rng, dt);
        // insert into the hashmap the location on the cartesian grid (at the end we will count how
        // many occurrences there are in the grid to compute the mean covered territory)
        if matches!(self.status, Status::Bulk(_)) {
            // TODO implement the deleting of visited squares - let's do a dumb version for now
            // where we check everyone
            let complex_position = self.position_to_complex();

            // for each square on the grid we check if it's not further than the particle size. If
            // it's within reach, we remove it from the set.
            // The important thing is that in the end we will onnly have unvisited guys
            total_territory.drain_filter(|_, grid_point| {
                (complex_position - *grid_point).norm() < particle_size
            });
        }
    }
}
