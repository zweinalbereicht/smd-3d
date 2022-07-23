use rand::distributions::Distribution;

use rand_distr::Uniform;
use rand_distr::{Exp, Normal};
use rand_pcg::*;

// algebra package
use na::*;

use num_complex::{Complex, Complex64};
use std::collections::HashMap;

use std::f64::consts::PI;
use std::fmt;
use std::process::exit;

use crate::ejectionenvironement::*;
use crate::tools::*;
use crate::traps::*;

pub struct EjectionParticle {
    pub position: na::Vector3<f64>, // we use these vectors cuz it will be easier to implement rotations.
    pub lifetime: f64,
    pub status: Status,
}

impl fmt::Display for EjectionParticle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "position : {}\nlifetime : {}\nstatus : {}",
            self.position,
            self.lifetime,
            match self.status {
                Status::Dead(_) => "dead",
                Status::Bulk(_) => "bulk",
                Status::Absorbed(_) => "absorbed",
            }
        )
    }
}

// default position is at the center of the sphere.
impl Default for EjectionParticle {
    fn default() -> Self {
        EjectionParticle {
            position: na::vector![0., 0., 0.],
            lifetime: 0.,
            status: Status::Bulk(PointlikeTarket::NotCrossed),
        }
    }
}

impl EjectionParticle {
    pub fn reset(&mut self, position: Vector3<f64>) {
        self.position = position.clone();
        self.lifetime = 0.;
        self.status = Status::Bulk(PointlikeTarket::NotCrossed);
    }
}

// we will put it as a class function for portability
impl EjectionParticle {

    // moves the particle when it's absorbed on a sticky trap
    pub fn move_on_sticky_trap(
        &mut self,
        stickytrap: &StickyTrap,
        environement: &EjectionEnvironment,
        rng: &mut Lcg128Xsl64,
    ) {
        //// prepare all randomness needed
        let exp = Exp::new(1. / stickytrap.desorption_time).unwrap();
        let normal = Normal::new(0., 1.).unwrap();
        let elapsed_time = exp.sample(rng);
        // this move move_along_boundary returns tha engle of the rotation we will apply to
        // our vector. The 4 here comes from the fact that we are moving on a 2D surface
        let move_along_boundary = (4. * stickytrap.boundary_coefficient * elapsed_time).sqrt()
            * normal.sample(rng)
            / stickytrap.radius;

        //move --> it doesnt matter which rotation we take at this point since we are going
        //to be isotropic.
        let axisangle = move_along_boundary * new_sphere_vector(rng);
        let rotation = Rotation3::new(axisangle);

        // rotate vector on the sticky trap surface and go to new position by ejecting
        let mut position_on_stickytrap = self.position - stickytrap.center;
        position_on_stickytrap = (rotation * position_on_stickytrap)
            .scale((stickytrap.radius + environement.ejection_length) / stickytrap.radius);
        self.position = stickytrap.center + position_on_stickytrap;

        // add the elapsed time
        self.lifetime += elapsed_time;
        self.status = Status::Bulk(PointlikeTarket::NotCrossed)
    }

    // move the particle when it's absorbed on the boundary
    pub fn move_on_boundary(&mut self, environement: &EjectionEnvironment, rng: &mut Lcg128Xsl64) {
        //println!("in the absorbed phase");
        let exp = Exp::new(1. / environement.desorption_time).unwrap();
        let normal = Normal::new(0., 1.).unwrap();
        let elapsed_time = exp.sample(rng);
        // this move move_along_boundary returns tha engle of the rotation we will apply to
        // our vector. The 4 here comes from the fact that we are moving on a 2D surface
        let move_along_boundary = (4. * environement.boundary_coefficient * elapsed_time).sqrt()
            * normal.sample(rng)
            / environement.outer_radius;

        //move --> it doesnt matter which rotation we take at this point since we are going
        //to be isotropic.
        let axisangle = move_along_boundary * new_sphere_vector(rng);
        let rotation = Rotation3::new(axisangle);
        self.position = rotation * self.position;

        //eject
        self.position -= self.position.normalize() * environement.ejection_length;
        // add the elapsed time
        self.lifetime += elapsed_time;
        self.status = Status::Bulk(PointlikeTarket::NotCrossed)
    }

    // move the particle in the bulk if it's in the event driven algorithm.
    pub fn move_in_bulk_event(
        &mut self,
        environement: &EjectionEnvironment,
        tolerance: f64,
        rng: &mut Lcg128Xsl64,
    ) {
        // find smallest radius for the ball between bonudary and traps, and add a little tolerance.
        let smallest_radius = environement.find_smallest_radius(self.position, tolerance);
        //calculate the jump by drawing on a sphere and then multiplying by the min radius
        let increment = smallest_radius * new_sphere_vector(rng);

        //move the particle
        self.position += increment;

        //increment time by the mean time needed to exit the biggest circle possible -->
        //this value can be found in getoor and blumenthal.
        // the 6 here comes from the fact that we are 3D.
        self.lifetime += smallest_radius.powi(2) / (6. * environement.bulk_coefficient);

        // if we jump outside we get absorbed
        if self.position.norm() >= environement.outer_radius {
            self.position
                .scale_mut(environement.outer_radius / self.position.norm());
            self.status = Status::Absorbed(Surface::Boundary);
        }
        // if we are inside a trap, we need to stop the process. Status::Dead stores the
        // index of the absorbing trap
        else if let Some(absorbing_trap) = &environement
            .traps
            .iter()
            .position(|&t| single_trap_intersect(self.position, &t))
        {
            self.status = Status::Dead(*absorbing_trap);
        } else if let Some(stickytrap) = &environement
            .stickytraps
            .iter()
            .find(|&t| single_stickytrap_intersect(self.position, &t))
        {
            let position_on_stickytrap = self.position-stickytrap.center;
            self.position = stickytrap.center + position_on_stickytrap.scale(stickytrap.radius/position_on_stickytrap.norm());
            self.status = Status::Absorbed(Surface::StickyTrap(**stickytrap));
        }
    }

    pub fn move_in_bulk_discrete(
        &mut self,
        environement: &EjectionEnvironment,
        timestep: f64,
        rng: &mut Lcg128Xsl64,
    ) {
        //println!("in the bulk");

        //draw the move, the 6 here is because we are in 3D
        let normal = Normal::new(0., 1.).unwrap();
        let jump = (6. * environement.bulk_coefficient * timestep).sqrt()
            * normal.sample(rng)
            * new_sphere_vector(rng);

        //move the particle and increment time
        self.position += jump;
        self.lifetime += timestep;

        // if we jump outside we get absorbed
        if self.position.norm() >= environement.outer_radius {
            self.position
                .scale_mut(environement.outer_radius / self.position.norm());
            self.status = Status::Absorbed(Surface::Boundary);
        }
        // if we are inside a trap, we need to stop the process. Status::Dead stores the
        // index of the absorbing trap
        else if let Some(absorbing_trap) = &environement
            .traps
            .iter()
            .position(|&t| single_trap_intersect(self.position, &t))
        {
            self.status = Status::Dead(*absorbing_trap);
        } else if let Some(stickytrap) = &environement
            .stickytraps
            .iter()
            .find(|&t| single_stickytrap_intersect(self.position, &t))
        {
            let position_on_stickytrap = self.position-stickytrap.center;
            self.position = stickytrap.center + position_on_stickytrap.scale(stickytrap.radius/position_on_stickytrap.norm());
            self.status = Status::Absorbed(Surface::StickyTrap(**stickytrap));
        }
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
            Status::Absorbed(Surface::Boundary) => {
                self.move_on_boundary(environement, rng);
            }
            Status::Absorbed(Surface::StickyTrap(stickytrap)) => {
                self.move_on_sticky_trap(&stickytrap, environement, rng);
            }
            // here we implement the step by step discrete time algorithm.
            Status::Bulk(_) => {
                self.move_in_bulk_discrete(environement, dt, rng);
            }
        }
    }

    // we implement the event driven move algorithm where we just pick a random point on the
    // exterior of the biggest possible sphere. It is also useful to calculate quicker
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
            Status::Absorbed(Surface::Boundary) => {
                self.move_on_boundary(environement, rng);
            }
            Status::Absorbed(Surface::StickyTrap(stickytrap)) => {
                self.move_on_sticky_trap(&stickytrap, environement, rng);
            }
            // here we implement the algortihm to do the biggest centered sphere movement.
            Status::Bulk(_) => {
                self.move_in_bulk_event(environement, tolerance, rng);
            }
        }
    }
}

/* pub fn move_and_record(
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
} */
