use na::{vector, Vector3};
use num_complex::Complex64;
use rand_distr::{Distribution, Exp, Normal};
use std::collections::HashMap;

use crate::core_ejection_length::*;
use crate::ejectionenvironement::EjectionEnvironment;
use crate::traps::*;
use rand_pcg::*;
use std::cmp::Eq;
use std::io;

/* fn _float_from_cmd(msg: &str) -> f64 {
    println!("{}", msg);
    let mut input = String::new();

    io::stdin().read_line(&mut input).expect("entry not valid");
    input.trim().parse::<f64>().expect("entry not valid")
} */

pub fn trap_intersect(point: na::Vector3<f64>, traps: &Vec<Trap>) -> bool {
    traps.iter().any(|t| (t.center - point).norm() < t.radius)
}

pub fn single_trap_intersect(point: na::Vector3<f64>, t: Trap) -> bool {
    (t.center - point).norm() < t.radius
}

//creates a random vector on the unit sphere
pub fn new_sphere_vector(rng: &mut Lcg128Xsl64) -> Vector3<f64> {
    let mut var = 1. / 3.;
    let normal = Normal::new(0., var).unwrap();
    let mut vector_on_sphere = vector![0., 0., 0.];
    for i in 0..vector_on_sphere.len() {
        vector_on_sphere[i] = normal.sample(rng);
    }
    vector_on_sphere.normalize()
}

// keeps track of if pointlike target was crossed
pub enum PointlikeTarket {
    Crossed,
    NotCrossed,
}

// keeps track of the status of the particle
pub enum Status {
    Bulk(PointlikeTarket),
    Absorbed,
    Dead(usize),
}

/* // helps to record visited space
#[derive(Hash, PartialEq, Eq)]
pub struct CartesianGridLocation {
    pub x: i32,
    pub y: i32,
}

impl CartesianGridLocation {
    // find the position on the grid corresponding to the location f the walker
    pub fn new_from_particle(particle: &EjectionParticle, gridsize: f64) -> Self {
        let x_projection = ((particle.radial_position * particle.angular_position.cos()) / gridsize)
            .floor() as i32;
        let y_projection = ((particle.radial_position * particle.angular_position.sin()) / gridsize)
            .floor() as i32;
        CartesianGridLocation {
            x: x_projection,
            y: y_projection,
        }
    }

    // construct position on grid from a complex number
    pub fn new_from_complex(position: &Complex64, gridsize: f64) -> Self {
        let x_projection = ((position.norm() * position.arg().cos()) / gridsize).floor() as i32;
        let y_projection = ((position.norm() * position.arg().sin()) / gridsize).floor() as i32;
        CartesianGridLocation {
            x: x_projection,
            y: y_projection,
        }
    }

    // construct complex from position on grid
    pub fn to_complex(&self, gridsize: f64) -> Complex64 {
        Complex64::new(gridsize * self.x as f64, gridsize * self.y as f64)
    }

    pub fn build_map(
        environement: &EjectionEnvironment,
        gridsize: f64,
    ) -> HashMap<CartesianGridLocation, Complex64> {
        let mut territory_map = HashMap::<CartesianGridLocation, Complex64>::new();
        let mut x_step = 0;
        let mut y_step = 0;
        for _ in arrange::FloatRange::new(
            -environement.outer_radius,
            environement.outer_radius,
            gridsize,
        )
        .range()
        .iter()
        {
            for _ in arrange::FloatRange::new(
                -environement.outer_radius,
                environement.outer_radius,
                gridsize,
            )
            .range()
            .iter()
            {
                // we store in the map the complex number associated to the mesh, if it' not in one
                // of the traps or out of the cirlce.
                let complex_point = Complex64::new(
                    -environement.outer_radius + gridsize * x_step as f64,
                    -environement.outer_radius + gridsize * y_step as f64,
                );

                if (complex_point.norm() < environement.outer_radius)
                    && (!trap_intersect(complex_point, &environement.traps))
                {
                    territory_map.insert(
                        CartesianGridLocation {
                            x: x_step,
                            y: y_step,
                        },
                        complex_point,
                    );
                }
                y_step += 1;
            }
            y_step = 0;
            x_step += 1;
        }
        territory_map
    }
} */
