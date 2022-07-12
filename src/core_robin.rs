use rand::distributions::Distribution;
use rand::Rng;
use rand_distr::{Exp, Normal};
use rand_pcg::*;

use num_complex::Complex;
use std::f64::consts::PI;
use std::fmt;
use std::process::exit;

use crate::tools::*;
// we will track our particle according to radial and angular position

pub struct Particle {
    pub radial_position: f64,
    pub angular_position: f64,
    pub lifetime: f64,
    pub status: Status,
}

impl fmt::Display for Particle {
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

impl Default for Particle {
    fn default() -> Self {
        Particle {
            radial_position: 1.,
            angular_position: 0.,
            lifetime: 0.,
            status: Status::Bulk(PointlikeTarket::NotCrossed),
        }
    }
}

impl Particle {
    pub fn reset(&mut self, radial_pos: f64, angular_pos: f64) {
        self.radial_position = radial_pos;
        self.angular_position = angular_pos;
        self.lifetime = 0.;
    }
}

// stores the info on our environement
// the inner circle wil always be considered to be offseted toward the 'left'
pub struct Environment {
    pub inner_radius: f64,
    pub outer_radius: f64,
    pub center_offset: f64,
    //diffusive coefficents
    pub bulk_coefficient: f64,
    pub boundary_coefficient: f64,
    // probability of being absorbed upon boundary encounter
    pub absorption_probability: f64,
    // time spend on the boundary
    pub desorption_time: f64,
}

// the default is a centered inner target
impl Default for Environment {
    fn default() -> Self {
        Environment {
            inner_radius: 1.,
            outer_radius: 10.,
            center_offset: 0.,
            bulk_coefficient: 1.0,
            boundary_coefficient: 1.0,
            absorption_probability: 0.,
            desorption_time: 1.,
        }
    }
}

impl Particle {
    pub fn move_particle(&mut self, environement: &Environment, rng: &mut Lcg128Xsl64, dt: f64) {
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

                self.angular_position =
                    Complex::from_polar(environement.outer_radius, new_pos).arg();
                self.lifetime += elapsed_time;
            }

            // in here we are going to implement the algortithm that moves the particle to the boundary of a decentered circle, that will intersect
            // both the outer rim and the inner circle. This is based on some simulation algorithms devised by other people, see overleaf for more
            // info
            Status::Bulk(_) => {
                //println!("in the bulk");

                // convert to complex for easier use
                let mut particle_position_complex =
                    Complex::from_polar(self.radial_position, self.angular_position);
                let inner_sphere_center_complex = Complex::new(-environement.center_offset, 0.); // watch out for the minus sign here !

                //implement the jump
                let normal = Normal::new(0., 1.).unwrap();
                let increment = Complex::new(
                    (2. * environement.bulk_coefficient * dt).sqrt() * normal.sample(rng),
                    (2. * environement.bulk_coefficient * dt).sqrt() * normal.sample(rng),
                );

                //move the particle
                particle_position_complex += increment;
                self.lifetime += dt;

                // if we jump outside
                if particle_position_complex.norm() >= environement.outer_radius {
                    let projection_length: f64 =
                        particle_position_complex.norm() - environement.outer_radius;

                    //implement the absorption condition
                    // if we are rejected
                    let reflection_probability = (1.
                        - environement.absorption_probability
                            * (PI * dt / environement.bulk_coefficient).sqrt())
                    .max(0.); // just in case it's bigger than 1

                    if rng.gen_bool(reflection_probability) {
                        self.angular_position = particle_position_complex.arg();
                        self.radial_position =
                            particle_position_complex.norm() - 2. * projection_length;
                    }
                    // if we are absorbed
                    else {
                        self.angular_position = particle_position_complex.arg();
                        self.radial_position = environement.outer_radius;
                        self.status = Status::Absorbed;
                    }
                } else {
                    self.angular_position = particle_position_complex.arg();
                    self.radial_position = particle_position_complex.norm();
                }

                // check if we are dead cuz it can happen with a rebound for instance
                if (particle_position_complex - inner_sphere_center_complex).norm()
                    < environement.inner_radius
                {
                    self.status = Status::Dead(0);
                }
            }
        }
    }
}
