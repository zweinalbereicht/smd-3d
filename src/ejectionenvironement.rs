use na::Vector3;

use crate::traps::*;
// stores the info on our environement
// the inner circle wil always be considered to be offseted toward the 'left'
#[derive(Debug)]
pub struct EjectionEnvironment {
    pub outer_radius: f64,
    //diffusive coefficents
    pub bulk_coefficient: f64,
    pub boundary_coefficient: f64,
    // how far we are ejected from boundary
    pub ejection_length: f64,
    // time spend on the boundary
    pub desorption_time: f64,

    //traps
    pub traps: Vec<Trap>,
    //stickytraps
    pub stickytraps: Vec<StickyTrap>,
}

// the default is a centered inner target ejected at a distance of 1.
impl Default for EjectionEnvironment {
    fn default() -> Self {
        EjectionEnvironment {
            outer_radius: 10.,
            bulk_coefficient: 1.0,
            boundary_coefficient: 1.0,
            ejection_length: 1.,
            desorption_time: 1.,
            // prepare empty vector for traps
            traps: vec![],
            stickytraps: vec![],
        }
    }
}

//implemeting the add trap and add sticky trap function
impl EjectionEnvironment {
    pub fn add_trap(&mut self, trap: Trap) {
        self.traps.push(trap);
    }
    pub fn add_sticky_trap(&mut self, stickytrap: StickyTrap) {
        self.stickytraps.push(stickytrap);
    }

    // find the smallest raius from position to the various holes in the system, and the general
    // bounding domain. Surely there's a nicer way to do this but...
    pub fn find_smallest_radius(&self, position: Vector3<f64>, tolerance: f64) -> f64 {
        let smallest_radius = ((self.outer_radius - position.norm()).min(
            self.traps
                .iter()
                .map(|t| (position - t.center).norm() - t.radius)
                .reduce(f64::min)
                .unwrap_or(self.outer_radius + 1.),
        ))
        .min(
            self.stickytraps
                .iter()
                .map(|t| (position - t.center).norm() - t.radius)
                .reduce(f64::min)
                .unwrap_or(self.outer_radius + 1.), // Ensure nothing goes wrong if we have no traps
        );
        smallest_radius + tolerance
    }
}
