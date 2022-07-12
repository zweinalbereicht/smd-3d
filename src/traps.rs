use num_complex::Complex;

#[derive(Clone,Copy,Debug)]
pub struct Trap {
    pub radius: f64,
    pub center: Complex<f64>,
}

// the default is a centered target with radius one
impl Default for Trap {
    fn default() -> Self {
        Trap {
            radius: 1.,
            center: Complex::new(0., 0.),
        }
    }
}

// this function parses a vector containing 3 numbers [radial_distance, angle, radius]
pub fn parse_trap(trap : &Vec<f64>) -> Option<Trap> {
    Some(Trap{
        radius : trap[2],
        center : Complex::from_polar(trap[0],trap[1]) 
    })
}
