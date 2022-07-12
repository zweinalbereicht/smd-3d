use std::collections::HashMap;
use std::collections::HashSet;

use crate::core_ejection_length::*;
use crate::core_robin::*;
use crate::ejectionenvironement::*;
use crate::tools::*;
use itertools::Itertools;
use num_complex::Complex64;
use rand_pcg::Lcg128Xsl64;

// computes one first passge time to inner sphere
fn first_passage_time_to_inner_sphere_robin(
    initial_radial_pos: f64,
    initial_angular_pos: f64,
    environement: &Environment,
    rng: &mut Lcg128Xsl64,
    dt: f64,
) -> f64 {
    let mut particle = Particle::default();
    particle.reset(initial_radial_pos, initial_angular_pos);
    // check if we are already on the boundary
    particle.status = match initial_radial_pos >= environement.outer_radius {
        true => Status::Absorbed,
        false => Status::Bulk(PointlikeTarket::NotCrossed),
    };
    while !(matches!(particle.status, Status::Dead(_))) {
        particle.move_particle(&environement, rng, dt);
    }
    //println!("{}", particle);
    particle.lifetime
}
// creates a distribution of fpt
pub fn fpt_distribution_robin(
    initial_radial_pos: f64,
    initial_angular_pos: f64,
    number_of_simulations: u32,
    environement: &Environment,
    mut rng: &mut Lcg128Xsl64,
    dt: f64,
) -> Vec<f64> {
    vec![0.; number_of_simulations as usize]
        .iter()
        .map(|_| {
            first_passage_time_to_inner_sphere_robin(
                initial_radial_pos,
                initial_angular_pos,
                &environement,
                &mut rng,
                dt,
            )
        })
        .collect_vec()
}

// computes one first passge time to inner sphere
fn first_passage_time_to_inner_sphere_ejection(
    initial_radial_pos: f64,
    initial_angular_pos: f64,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    dt: f64,
) -> f64 {
    let mut particle = EjectionParticle::default();
    particle.reset(initial_radial_pos, initial_angular_pos);
    // check if we are already on the boundary
    particle.status = match initial_radial_pos >= environement.outer_radius {
        true => Status::Absorbed,
        false => Status::Bulk(PointlikeTarket::NotCrossed),
    };
    while !(matches!(particle.status, Status::Dead(_))) {
        particle.move_particle(&environement, rng, dt);
    }
    //println!("{}", particle);
    particle.lifetime
}

// mean first passage time to inner sphere
pub fn mean_first_passage_time_to_inner_sphere(
    initial_position: f64,
    initial_angle: f64,
    nb_simulations: u32,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    tolerance: f64,
) -> f64 {
    let mut mfpt = 0.;
    let mut particle = EjectionParticle::default();
    for _ in 0..nb_simulations {
        particle.reset(initial_position, initial_angle);
        particle.status = match particle.radial_position < environement.outer_radius {
            true => Status::Bulk(PointlikeTarket::NotCrossed),
            false => Status::Absorbed,
        };

        // until we touch a target
        while !(matches!(particle.status, Status::Dead(_))) {
            particle.move_particle_event(&environement, rng, tolerance);
        }
        mfpt += particle.lifetime;
    }
    mfpt / (nb_simulations as f64)
}

// computes the are of covered territory before hitting one of the traps.
fn covered_territory(
    initial_radial_pos: f64,
    initial_angular_pos: f64,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    dt: f64,
    particle_size: f64,
    gridsize: f64,
) -> f64 {
    // creates the hashmap that will serve as a territory recorder and fill it.

    let mut total_territory = CartesianGridLocation::build_map(environement, gridsize);

    let true_number_cell = total_territory.len();

    let mut particle = EjectionParticle::default();
    particle.reset(initial_radial_pos, initial_angular_pos);

    // check if we are already on the boundary
    particle.status = match initial_radial_pos >= environement.outer_radius {
        true => Status::Absorbed,
        false => Status::Bulk(PointlikeTarket::NotCrossed),
    };

    // let the particle evolve
    while !(matches!(particle.status, Status::Dead(_))) {
        particle.move_and_record(&environement, rng, dt, particle_size, &mut total_territory);
    }

    // we now need to evaluate how many "cells" of side grid_size have been visited.
    // return the territory covered in the shape of (nb cells)*area of cell
    (true_number_cell - total_territory.len()) as f64 * gridsize.powi(2)
}

// computes a distribution of the first deorption angle for a given starting point.
fn escape_angle_from_outer_sphere(
    initial_radial_pos: f64,
    initial_angular_pos: f64,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    dt: f64,
) -> f64 {
    let mut particle = EjectionParticle::default();
    particle.reset(initial_radial_pos, initial_angular_pos);
    // we need to handle the case where we start on the boundary, done through this match function
    match initial_radial_pos >= environement.outer_radius {
        true => {
            particle.status = Status::Absorbed; // its important we do this to make sure the move is along the boundary
            particle.move_particle(&environement, rng, dt);
            particle.angular_position
        }
        false => {
            let mut absorbed_once = false;
            while !(absorbed_once) {
                match particle.status {
                    Status::Bulk(_) => particle.move_particle(&environement, rng, dt),
                    Status::Dead(_) => particle.reset(initial_radial_pos, initial_angular_pos),
                    Status::Absorbed => {
                        absorbed_once = true;
                        particle.move_particle(&environement, rng, dt);
                    }
                }
            }
            particle.angular_position
        }
    }
}

// creates a distribution of fpt
pub fn fpt_distribution_ejection(
    initial_radial_pos: f64,
    initial_angular_pos: f64,
    number_of_simulations: u32,
    environement: &EjectionEnvironment,
    mut rng: &mut Lcg128Xsl64,
    dt: f64,
) -> Vec<f64> {
    vec![0.; number_of_simulations as usize]
        .iter()
        .map(|_| {
            first_passage_time_to_inner_sphere_ejection(
                initial_radial_pos,
                initial_angular_pos,
                &environement,
                &mut rng,
                dt,
            )
        })
        .collect_vec()
}

pub fn escape_angle_distribution(
    initial_radial_pos: f64,
    initial_angular_pos: f64,
    number_of_simulations: u32,
    environement: &EjectionEnvironment,
    mut rng: &mut Lcg128Xsl64,
    dt: f64,
) -> Vec<f64> {
    vec![0.; number_of_simulations as usize]
        .iter()
        .map(|_| {
            escape_angle_from_outer_sphere(
                initial_radial_pos,
                initial_angular_pos,
                &environement,
                &mut rng,
                dt,
            )
        })
        .collect_vec()
}

pub fn covered_territory_distribution(
    initial_radial_pos: f64,
    initial_angular_pos: f64,
    number_of_simulations: u32,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    dt: f64,
    particle_size: f64,
    gridsize: f64,
) -> Vec<f64> {
    vec![0.; number_of_simulations as usize]
        .iter()
        .map(|_| {
            covered_territory(
                initial_radial_pos,
                initial_angular_pos,
                environement,
                rng,
                dt,
                particle_size,
                gridsize,
            )
        })
        .collect_vec()
}
// splitting probability to outer rim instead of targets
pub fn splitting_boundary_targets(
    initial_position: f64,
    initial_angle: f64,
    nb_simulations: u32,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    tolerance: f64,
) -> f64 {
    let mut splitting = 0.;
    for _ in 0..nb_simulations {
        let mut particle = EjectionParticle::default();
        particle.reset(initial_position, initial_angle);
        particle.status = Status::Bulk(PointlikeTarket::NotCrossed);

        // until we either toush a target or the boundary
        while !(matches!(particle.status, Status::Dead(_) | Status::Absorbed)) {
            particle.move_particle_event(&environement, rng, tolerance);
        }
        if matches!(particle.status, Status::Absorbed) {
            splitting += 1.;
        }
    }
    splitting / (nb_simulations as f64)
}

// splitting probability to one of the inner targets.
pub fn splitting_targets(
    initial_position: f64,
    initial_angle: f64,
    target: u32,
    nb_simulations: u32,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    tolerance: f64,
) -> f64 {
    let mut splitting = 0.;
    let mut particle = EjectionParticle::default();
    for _ in 0..nb_simulations {
        particle.reset(initial_position, initial_angle);
        particle.status = match particle.radial_position < environement.outer_radius {
            true => Status::Bulk(PointlikeTarket::NotCrossed),
            false => Status::Absorbed,
        };

        // until we either touch a target or the boundary
        while !(matches!(particle.status, Status::Dead(_))) {
            particle.move_particle_event(&environement, rng, tolerance);
        }
        let touchedtarget = match particle.status {
            Status::Dead(p) => p as u32,
            _ => (environement.traps.len() + 1) as u32, // just to make sure it doesn't happen
        };

        if target == touchedtarget {
            splitting += 1.;
        }
    }
    splitting / (nb_simulations as f64)
}

pub fn splitting_pointlike_bulk(
    initial_position: f64,
    initial_angle: f64,
    target: i32,
    nb_simulations: u32,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    tolerance: f64,
) -> f64 {
    let mut splitting = 0.;
    let mut particle = EjectionParticle::default();
    for _ in 0..nb_simulations {
        particle.reset(initial_position, initial_angle);
        particle.status = match particle.radial_position < environement.outer_radius {
            true => Status::Bulk(PointlikeTarket::NotCrossed),
            false => Status::Absorbed,
        };

        // until we either touch a target or the pointlike target on the boundary.
        while !(matches!(
            particle.status,
            Status::Dead(_) | Status::Bulk(PointlikeTarket::Crossed)
        )) {
            particle.move_particle_event(&environement, rng, tolerance);
        }
        let touchedtarget = match particle.status {
            Status::Dead(p) => p as i32,
            Status::Bulk(PointlikeTarket::Crossed) => -1, // this is we've been absorbed on the pointlike target
            _ => (environement.traps.len() + 1) as i32,   // just to make sure it doesn't happen
        };

        if target == touchedtarget {
            splitting += 1.;
        }
    }
    splitting / (nb_simulations as f64)
}
