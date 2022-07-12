#[derive(Clone,Copy,Debug)]
pub struct Trap {
    pub center : na::Vector3<f64>,
    pub radius: f64,
}

// the default is a centered target with radius one
impl Default for Trap {
    fn default() -> Self {
        Trap {
            center : na::vector![0.,0.,0.],
            radius: 1.,
        }
    }
}

// this function parses a vector containing 4 numbers [x_pos, y_pos, z_pos, radius]
pub fn parse_trap(trap : &Vec<f64>) -> Option<Trap> {
    Some(Trap{
        radius : trap[3],
        center : na::vector![trap[0],trap[1],trap[2]]
    })
}
