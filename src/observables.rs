use std::collections::HashMap;
use std::collections::HashSet;

use crate::core_ejection_length::*;
use crate::ejectionenvironement::*;
use crate::tools::*;
use itertools::Itertools;
use na::Vector3;
use rand_pcg::Lcg128Xsl64;

// computes one first passage time to inner sphere
// TODO(once we impement the step by step simulation, needed for distributions)
fn first_passage_time_to_inner_sphere(
    initial_position : Vector3<f64>,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    dt: f64,
) -> f64 {
    let mut particle = EjectionParticle::default();
    particle.reset(initial_position);
    // check if we are already on the boundary
    particle.status = match particle.position.norm() >= environement.outer_radius {
        true => Status::Absorbed(Surface::Boundary),
        false => Status::Bulk(PointlikeTarket::NotCrossed),
    };
    while !(matches!(particle.status, Status::Dead(_))) {
        particle.move_particle(&environement, rng, dt);
    }
    //println!("{}", particle);
    particle.lifetime
}


// creates a distribution of fpt
// Note that this cannot be done by using the event driven algorithm, so it' necessarily
// longer...we could implement an event driven version at some point I guess. --> need to think
// about it.
pub fn fpt_distribution(
    initial_position: Vector3<f64>,
    number_of_simulations: u32,
    environement: &EjectionEnvironment,
    mut rng: &mut Lcg128Xsl64,
    dt: f64,
) -> Vec<f64> {
    vec![0.; number_of_simulations as usize]
        .iter()
        .map(|_| {
            first_passage_time_to_inner_sphere(
                initial_position,
                &environement,
                &mut rng,
                dt,
            )
        })
        .collect_vec()
}

// mean first passage time to inner sphere
pub fn mean_first_passage_time_to_inner_sphere(
    initial_position: Vector3<f64>,
    nb_simulations: u32,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    tolerance: f64,
) -> f64 {
    let mut mfpt = 0.;
    let mut particle = EjectionParticle::default();
    for _ in 0..nb_simulations {
        particle.reset(initial_position);
        particle.status = match particle.position.norm() < environement.outer_radius {
            true => Status::Bulk(PointlikeTarket::NotCrossed),
            false => Status::Absorbed(Surface::Boundary),
        };

        // until we touch a target, we use the event algorithm.
        while !(matches!(particle.status, Status::Dead(_))) {
            particle.move_particle_event(&environement, rng, tolerance);
        }
        mfpt += particle.lifetime;
    }
    mfpt / (nb_simulations as f64)
}

// computes the are of covered territory before hitting one of the traps.
// TODO (later)
/* fn covered_territory(
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
*/

// splitting probability to outer rim instead of targets. note that evolving on a sticky traps
// doesn't count as reaching the boundary.
pub fn splitting_boundary_targets(
    initial_position: Vector3<f64>,
    nb_simulations: u32,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    tolerance: f64,
) -> f64 {
    let mut splitting = 0.;
    for _ in 0..nb_simulations {
        let mut particle = EjectionParticle::default();
        particle.reset(initial_position);
        particle.status = Status::Bulk(PointlikeTarket::NotCrossed);

        // until we either touch a target or the boundary
        while !(matches!(particle.status, Status::Dead(_) | Status::Absorbed(Surface::Boundary))) {
            particle.move_particle_event(&environement, rng, tolerance);
        }
        if matches!(particle.status, Status::Absorbed(Surface::Boundary)) {
            splitting += 1.;
        }
    }
    splitting / (nb_simulations as f64)
}

// splitting probability to one of the inner targets.
pub fn splitting_targets(
    initial_position: Vector3<f64>,
    target: u32,
    nb_simulations: u32,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    tolerance: f64,
) -> f64 {
    let mut splitting = 0.;
    let mut particle = EjectionParticle::default();
    for _ in 0..nb_simulations {
        particle.reset(initial_position);
        particle.status = match particle.position.norm() < environement.outer_radius {
            true => Status::Bulk(PointlikeTarket::NotCrossed),
            false => Status::Absorbed(Surface::Boundary),
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

// returns the state of the particle after some user inputted time has elapsed
// Since we are using the event driven algorithm here, this is going to be only approximate.
pub fn stationnary_state(
    initial_position: Vector3<f64>,
    evolve_time: f64,
    environement: &EjectionEnvironment,
    rng: &mut Lcg128Xsl64,
    nb_simulations: u32,
    tolerance: f64,
) -> Vec<String> {
    let mut particle = EjectionParticle::default();
    let mut old_status = particle.status.clone();
    let mut states: Vec<String> = vec![];
    for _ in 0..nb_simulations {
        particle.reset(initial_position);
        while particle.lifetime < evolve_time {
            old_status = particle.status.clone();
            particle.move_particle_event(environement, rng, tolerance);
        }
        states.push(match old_status {
            Status::Bulk(_) => String::from("bulk"),
            Status::Absorbed(_) => String::from("boundary"),
            Status::Dead(_) => String::from("dead"),
        })
    }
    print!("{}", states.len());
    states
}
