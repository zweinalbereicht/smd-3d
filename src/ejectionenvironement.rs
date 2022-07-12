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
    pub traps : Vec<Trap>,
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
            // prepare empty vector for traps and fill it with one centered trap.
            traps:vec![Trap::default()],
        }
    }
}

//implemeting the add trap function
impl EjectionEnvironment{
    pub fn add_trap(&mut self, trap : Trap) {
        self.traps.push(trap);
    }
}
