#![feature(default_free_fn)]
#![feature(hash_drain_filter)]
mod core_ejection_length;
mod core_robin;
mod ejectionenvironement;
mod observables;
mod tools;
mod traps;
use std::default::default;
use std::env;
use crate::core_ejection_length::EjectionParticle;
use crate::core_robin::*;
use crate::ejectionenvironement::*;
use crate::observables::*;
use crate::traps::*;
use rand::SeedableRng;
use pyo3::prelude::*;

// our user function for the fpt
#[pyfunction]
fn fpt_distribution_robin_python(
    inner_radius: f64,
    outer_radius: f64,
    center_offset: f64,
    bulk_coefficient: f64,
    boundary_coefficient: f64,
    absorption_probability: f64, // this is the kappa in the flux at the boundary
    desorption_time: f64,
    initial_position: f64,
    initial_angle: f64,
    nb_simulations: u32,
    timestep: f64,
) -> PyResult<Vec<f64>> {
    let environement = Environment {
        inner_radius,
        outer_radius,
        center_offset,
        bulk_coefficient,
        boundary_coefficient,
        absorption_probability,
        desorption_time,
    };

    let mut rng = rand_pcg::Pcg64::from_entropy();

    let fpts = fpt_distribution_robin(
        initial_position,
        initial_angle,
        nb_simulations,
        &environement,
        &mut rng,
        timestep,
    );

    Ok(fpts)
}

#[pyfunction(text_signature = "
computes a distribution of fpt in the SMD ejection setup, with multiple traps.

parameters : 

outer_radius,
bulk_coefficient,
boundary_coefficient,
ejection_length,
desorption_time,
traps, [[r1, a1, radius1],[radial2,angle2,radius2],...]
initial_position,
initial_angle,
nb_simulations,
timestep,
")]
fn fpt_distribution_ejection_python(
    outer_radius: f64,
    bulk_coefficient: f64,
    boundary_coefficient: f64,
    ejection_length: f64,
    desorption_time: f64,
    traps: Vec<Vec<f64>>,
    initial_position: f64,
    initial_angle: f64,
    nb_simulations: u32,
    timestep: f64,
) -> PyResult<Vec<f64>> {
    let mut environement = EjectionEnvironment {
        outer_radius,
        bulk_coefficient,
        boundary_coefficient,
        ejection_length,
        desorption_time,
        ..default()
    };

    if traps.len() > 0 {
        // take away the default trap
        environement.traps.pop();
        // add all the new traps
        traps
            .iter()
            .map(|x| parse_trap(x))
            .for_each(|p| environement.add_trap(p.unwrap()));
    }

    let mut rng = rand_pcg::Pcg64::from_entropy();

    let fpts = fpt_distribution_ejection(
        initial_position,
        initial_angle,
        nb_simulations,
        &environement,
        &mut rng,
        timestep,
    );

    Ok(fpts)
}

#[pyfunction(text_signature = "
computes the mfpt in the SMD ejection setup, using an event driven algorithm.

parameters : 

outer_radius,
bulk_coefficient,
boundary_coefficient,
ejection_length,
desorption_time,
traps, [[r1, a1, radius1],[radial2,angle2,radius2],...]
initial_position,
initial_angle,
nb_simulations,
tolerance,
")]
fn mean_fpt_ejection_python(
    outer_radius: f64,
    bulk_coefficient: f64,
    boundary_coefficient: f64,
    ejection_length: f64,
    desorption_time: f64,
    traps: Vec<Vec<f64>>,
    initial_position: f64,
    initial_angle: f64,
    nb_simulations: u32,
    tolerance: f64,
) -> PyResult<f64> {
    let mut environement = EjectionEnvironment {
        outer_radius,
        bulk_coefficient,
        boundary_coefficient,
        ejection_length,
        desorption_time,
        ..default()
    };

    if traps.len() > 0 {
        // take away the default trap
        environement.traps.pop();
        // add all the new traps
        traps
            .iter()
            .map(|x| parse_trap(x))
            .for_each(|p| environement.add_trap(p.unwrap()));
    }

    let mut rng = rand_pcg::Pcg64::from_entropy();

    let mfpt = mean_first_passage_time_to_inner_sphere(
        initial_position,
        initial_angle,
        nb_simulations,
        &environement,
        &mut rng,
        tolerance,
    );

    Ok(mfpt)
}

#[pyfunction(text_signature = "
computes the probability to hit the boundary before any of the inner targets

parameters : 

outer_radius,
ejection_length,
desorption_time,
boundary_coefficient
traps, [[r1, a1, radius1],[radial2,angle2,radius2],...]
initial_position,
initial_angle,
nb_simulations,
tolerance,
")]
fn splitting_boundary_targets_python(
    outer_radius: f64,
    ejection_length: f64,
    desorption_time: f64,
    boundary_coefficient: f64,
    traps: Vec<Vec<f64>>,
    initial_position: f64,
    initial_angle: f64,
    nb_simulations: u32,
    tolerance: f64,
) -> PyResult<f64> {
    let mut environement = EjectionEnvironment {
        outer_radius,
        bulk_coefficient: 1.,
        boundary_coefficient,
        ejection_length,
        desorption_time,
        ..default()
    };
    if traps.len() > 0 {
        // take away the default trap
        environement.traps.pop();
        // add all the new traps
        traps
            .iter()
            .map(|x| parse_trap(x))
            .for_each(|p| environement.add_trap(p.unwrap()));
    }
    let mut rng = rand_pcg::Pcg64::from_entropy();
    let splitting_probabilities = splitting_boundary_targets(
        initial_position,
        initial_angle,
        nb_simulations,
        &environement,
        &mut rng,
        tolerance,
    );
    Ok(splitting_probabilities)
}

#[pyfunction(text_signature = "
computes the probability to hit a specific target before the others

parameters : 

outer_radius,
ejection_length,
desorption_time,
boundary_coefficient
traps, [[r1, a1, radius1],[radial2,angle2,radius2],...]
initial_position,
initial_angle,
target, 
nb_simulations,
tolerance,
")]
fn splitting_targets_python(
    outer_radius: f64,
    ejection_length: f64,
    desorption_time: f64,
    boundary_coefficient: f64,
    traps: Vec<Vec<f64>>,
    initial_position: f64,
    initial_angle: f64,
    target: u32, // the target we are aiming for
    nb_simulations: u32,
    tolerance: f64,
) -> PyResult<f64> {
    let mut environement = EjectionEnvironment::default();
    environement.outer_radius = outer_radius;
    environement.ejection_length = ejection_length;
    environement.desorption_time = desorption_time;
    environement.boundary_coefficient = boundary_coefficient;

    // add the traps to the environement
    if traps.len() > 0 {
        // take away the default trap
        environement.traps.pop();
        // add all the new traps
        traps
            .iter()
            .map(|x| parse_trap(x))
            .for_each(|p| environement.add_trap(p.unwrap()));
    }

    let mut rng = rand_pcg::Pcg64::from_entropy();
    let splitting_probabilities = splitting_targets(
        initial_position,
        initial_angle,
        target,
        nb_simulations,
        &environement,
        &mut rng,
        tolerance,
    );
    Ok(splitting_probabilities)
}

#[pyfunction(text_signature = "
computes the probability to hit a specific bulk target before any of the others, or a pointlike
target located at theta=0 by default.

parameters : 

outer_radius,
ejection_length,
desorption_time,
boundary_coefficient
traps, [[r1, a1, radius1],[radial2,angle2,radius2],...]
initial_position,
initial_angle,
target, // if target == -1 it means we're looking for the pointlike one.
nb_simulations,
tolerance,
")]
fn splitting_pointlike_bulk_python(
    outer_radius: f64,
    ejection_length: f64,
    desorption_time: f64,
    boundary_coefficient: f64,
    traps: Vec<Vec<f64>>,
    initial_position: f64,
    initial_angle: f64,
    target: i32, // the target we are aiming for
    nb_simulations: u32,
    tolerance: f64,
) -> PyResult<f64> {
    let mut environement = EjectionEnvironment::default();
    environement.outer_radius = outer_radius;
    environement.ejection_length = ejection_length;
    environement.desorption_time = desorption_time;
    environement.boundary_coefficient = boundary_coefficient;

    // add the traps to the environement
    if traps.len() > 0 {
        // take away the default trap
        environement.traps.pop();
        // add all the new traps
        traps
            .iter()
            .map(|x| parse_trap(x))
            .for_each(|p| environement.add_trap(p.unwrap()));
    }

    let mut rng = rand_pcg::Pcg64::from_entropy();
    let splitting_probabilities = splitting_pointlike_bulk(
        initial_position,
        initial_angle,
        target,
        nb_simulations,
        &environement,
        &mut rng,
        tolerance,
    );
    Ok(splitting_probabilities)
}

#[pyfunction(text_signature = "
computes a distribution of covered territories

parameters : 

outer_radius,
ejection_length,
desorption_time,
boundary_coefficient
traps, [[r1, a1, radius1],[radial2,angle2,radius2],...]
initial_position,
initial_angle,
nb_simulations,
timestep,
particle_size, // this corresponds to the actual particle size
grid_size, // this corresponds to the size of the grid we use to check if the volume is visited or not.
")]
fn covered_territory_distribution_python(
    outer_radius: f64,
    ejection_length: f64,
    desorption_time: f64,
    boundary_coefficient: f64,
    traps: Vec<Vec<f64>>,
    initial_radial_pos: f64,
    initial_angular_pos: f64,
    nb_simulations: u32,
    dt: f64,
    particle_size: f64,
    gridsize: f64,
) -> PyResult<Vec<f64>> {
    let mut environement = EjectionEnvironment {
        outer_radius,
        bulk_coefficient: 1.,
        boundary_coefficient,
        ejection_length,
        desorption_time,
        ..default()
    };

    // add the traps to the environement
    if traps.len() > 0 {
        // take away the default trap
        environement.traps.pop();
        // add all the new traps
        traps
            .iter()
            .map(|x| parse_trap(x))
            .for_each(|p| environement.add_trap(p.unwrap()));
    }

    let mut rng = rand_pcg::Pcg64::from_entropy();
    let covered_territories = covered_territory_distribution(
        initial_radial_pos,
        initial_angular_pos,
        nb_simulations,
        &environement,
        &mut rng,
        dt,
        particle_size,
        gridsize,
    );

    Ok(covered_territories)
}

#[pyfunction]
fn escape_angle_distribution_python(
    outer_radius: f64,
    bulk_coefficient: f64,
    boundary_coefficient: f64,
    ejection_length: f64, // this is the kappa in the flux at the boundary
    //TODO add the traps parsing
    desorption_time: f64,
    initial_position: f64,
    initial_angle: f64,
    nb_simulations: u32,
    timestep: f64,
) -> PyResult<Vec<f64>> {
    let environement = EjectionEnvironment {
        outer_radius,
        bulk_coefficient,
        boundary_coefficient,
        ejection_length,
        desorption_time,
        ..default()
    };

    let mut rng = rand_pcg::Pcg64::from_entropy();

    let fpts = escape_angle_distribution(
        initial_position,
        initial_angle,
        nb_simulations,
        &environement,
        &mut rng,
        timestep,
    );

    Ok(fpts)
}

/// General callable functions from python module.
#[pymodule]
fn smd_event_driven(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(escape_angle_distribution_python, m)?)?;
    //m.add_function(wrap_pyfunction!(fpt_distribution_robin_python, m)?)?;
    m.add_function(wrap_pyfunction!(fpt_distribution_ejection_python, m)?)?;
    m.add_function(wrap_pyfunction!(splitting_boundary_targets_python, m)?)?;
    m.add_function(wrap_pyfunction!(splitting_targets_python, m)?)?;
    m.add_function(wrap_pyfunction!(splitting_pointlike_bulk_python, m)?)?;
    m.add_function(wrap_pyfunction!(covered_territory_distribution_python, m)?)?;
    m.add_function(wrap_pyfunction!(mean_fpt_ejection_python, m)?)?;
    Ok(())
}
