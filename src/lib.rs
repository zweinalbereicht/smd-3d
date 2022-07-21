#![feature(default_free_fn)]
#![feature(hash_drain_filter)]
mod core_ejection_length;
mod ejectionenvironement;
mod observables;
mod tools;
mod traps;
use crate::core_ejection_length::EjectionParticle;
use crate::ejectionenvironement::*;
use crate::observables::*;
use crate::traps::*;
use na::Vector3;
use pyo3::prelude::*;
use rand::SeedableRng;
use std::default::default;
use std::env;
extern crate nalgebra as na;

#[pyfunction(text_signature = "
computes a distribution of fpt in the SMD ejection setup, with possible multiple traps.

parameters :

outer_radius,
bulk_coefficient,
boundary_coefficient,
ejection_length,
desorption_time,
traps, [[x1,y1,z1, radius1],[x2,y2,z2,radius2],...]
initial_position, [x,y,z]
nb_simulations,
timestep,
")]
fn fpt_distribution_python(
    outer_radius: f64,
    bulk_coefficient: f64,
    boundary_coefficient: f64,
    ejection_length: f64,
    desorption_time: f64,
    traps: Vec<Vec<f64>>,
    initial_position: Vec<f64>,
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
    }
    // add all the new traps
    traps
        .iter()
        .map(|x| parse_trap(x))
        .for_each(|p| environement.add_trap(p.unwrap()));

    let mut rng = rand_pcg::Pcg64::from_entropy();

    let fpts = fpt_distribution(
        Vector3::from_vec(initial_position),
        nb_simulations,
        &environement,
        &mut rng,
        timestep,
    );

    Ok(fpts)
}

#[pyfunction(text_signature = "
computes the mfpt in the 3D SMD ejection setup, using an event driven algorithm.

parameters : 

outer_radius,
bulk_coefficient,
boundary_coefficient,
ejection_length,
desorption_time,
traps, [[x1,y1,z1,radius1],[x2,y2,z2,radius2],...]
initial_position, [x,y,z]
nb_simulations,
tolerance,
")]
fn mfpt_python(
    outer_radius: f64,
    bulk_coefficient: f64,
    boundary_coefficient: f64,
    ejection_length: f64,
    desorption_time: f64,
    traps: Vec<Vec<f64>>,
    initial_position: Vec<f64>,
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

    // fill with the traps
    traps
        .iter()
        .map(|x| parse_trap(x))
        .for_each(|p| environement.add_trap(p.unwrap()));

    let mut rng = rand_pcg::Pcg64::from_entropy();

    let mfpt = mean_first_passage_time_to_inner_sphere(
        Vector3::<f64>::new(
            initial_position[0],
            initial_position[1],
            initial_position[2],
        ),
        nb_simulations,
        &mut environement,
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
traps, [[x1,y1,z1,radius1],[x2,y2,z2,radius2],...]
initial_position, [x,y,z]
nb_simulations,
tolerance,
")]
fn splitting_boundary_targets_python(
    outer_radius: f64,
    ejection_length: f64,
    desorption_time: f64,
    boundary_coefficient: f64,
    traps: Vec<Vec<f64>>,
    initial_position: Vec<f64>,
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

    // fill with the traps
    traps
        .iter()
        .map(|x| parse_trap(x))
        .for_each(|p| environement.add_trap(p.unwrap()));

    let mut rng = rand_pcg::Pcg64::from_entropy();
    let splitting_probabilities = splitting_boundary_targets(
        Vector3::<f64>::new(
            initial_position[0],
            initial_position[1],
            initial_position[2],
        ),
        nb_simulations,
        &mut environement,
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
traps, [[x1,y1,z1,radius1],[x2,y2,z2,radius2],...]
initial_position, [x,y,z]
target, --> starts at 0
nb_simulations,
tolerance,
")]
fn splitting_targets_python(
    outer_radius: f64,
    ejection_length: f64,
    desorption_time: f64,
    boundary_coefficient: f64,
    traps: Vec<Vec<f64>>,
    initial_position: Vec<f64>,
    target: u32, // the target we are aiming for
    nb_simulations: u32,
    tolerance: f64,
) -> PyResult<f64> {
    let mut environement = EjectionEnvironment::default();
    environement.outer_radius = outer_radius;
    environement.ejection_length = ejection_length;
    environement.desorption_time = desorption_time;
    environement.boundary_coefficient = boundary_coefficient;

    // fill with the traps
    traps
        .iter()
        .map(|x| parse_trap(x))
        .for_each(|p| environement.add_trap(p.unwrap()));

    let mut rng = rand_pcg::Pcg64::from_entropy();
    let splitting_probabilities = splitting_targets(
        Vector3::<f64>::new(
            initial_position[0],
            initial_position[1],
            initial_position[2],
        ),
        target,
        nb_simulations,
        &mut environement,
        &mut rng,
        tolerance,
    );
    Ok(splitting_probabilities)
}

/* #[pyfunction(text_signature = "
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
} */

#[pyfunction(text_signature = "
Returns an array of statuses after some elapsed time, stickytraps are included.

parameters : 

outer_radius,
bulk_coefficient,
boundary_coefficient,
ejection_length,
desorption_time,
stickytraps, [[x1,y1,z1,radius1,desorption_time, diffusive coefficient]]
initial_position, [x,y,z]
evolve_time,
nb_simulations,
tolerance,
")]
pub fn stationnary_state_python(
    outer_radius: f64,
    bulk_coefficient: f64,
    boundary_coefficient: f64,
    ejection_length: f64,
    desorption_time: f64,
    stickytraps : Vec<Vec<f64>>,
    initial_position: Vec<f64>,
    evolve_time: f64,
    nb_simulations: u32,
    tolerance: f64,
) -> PyResult<Vec<String>> {
    let mut environement = EjectionEnvironment {
        outer_radius,
        bulk_coefficient,
        boundary_coefficient,
        ejection_length,
        desorption_time,
        ..default()
    };

    // fill with the stickytraps
    stickytraps
        .iter()
        .map(|x| parse_sticky_trap(x))
        .for_each(|p| environement.add_sticky_trap(p.unwrap()));

    let mut rng = rand_pcg::Pcg64::from_entropy();
    let initial_position_vector = Vector3::<f64>::new(
        initial_position[0],
        initial_position[1],
        initial_position[2],
    );

    let stationnary_states = stationnary_state(
        initial_position_vector,
        evolve_time,
        &environement,
        &mut rng,
        nb_simulations,
        tolerance,
    );

    Ok(stationnary_states)
}

#[pyfunction(text_signature = "
small tests on the rust side
")]
pub fn testing_stuff() -> PyResult<Vec<f64>> {
    let mut v = Vector3::from_vec(vec![1., 1., 1.]);
    v.scale_mut(2.);
    Ok(vec![v.x, v.y, v.z])
}

/// General callable functions from python module.
#[pymodule]
fn smd_3d(_py: Python, m: &PyModule) -> PyResult<()> {
    // m.add_function(wrap_pyfunction!(escape_angle_distribution_python, m)?)?;
    // m.add_function(wrap_pyfunction!(fpt_distribution_ejection_python, m)?)?;
    m.add_function(wrap_pyfunction!(splitting_boundary_targets_python, m)?)?;
    m.add_function(wrap_pyfunction!(splitting_targets_python, m)?)?;
    // m.add_function(wrap_pyfunction!(splitting_pointlike_bulk_python, m)?)?;
    // m.add_function(wrap_pyfunction!(covered_territory_distribution_python, m)?)?;
    m.add_function(wrap_pyfunction!(mfpt_python, m)?)?;
    m.add_function(wrap_pyfunction!(fpt_distribution_python, m)?)?;
    m.add_function(wrap_pyfunction!(stationnary_state_python, m)?)?;
    m.add_function(wrap_pyfunction!(testing_stuff, m)?)?;
    Ok(())
}
