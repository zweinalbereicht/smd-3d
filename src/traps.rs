// code for traps and sticky traps

#[derive(Clone, Copy, Debug)]
pub struct Trap {
    pub center: na::Vector3<f64>,
    pub radius: f64,
}

// the default is a centered target with radius one
impl Default for Trap {
    fn default() -> Self {
        Trap {
            center: na::vector![0., 0., 0.],
            radius: 1.,
        }
    }
}

// this function parses a vector containing 4 numbers [x_pos, y_pos, z_pos, radius]
pub fn parse_trap(trap: &Vec<f64>) -> Option<Trap> {
    match trap.len() == 5 {
        true => Some(Trap {
            radius: trap[3],
            center: na::vector![trap[0], trap[1], trap[2]],
        }),
        false => None,
    }
}

pub fn trap_intersect(point: na::Vector3<f64>, traps: &Vec<Trap>) -> bool {
    traps.iter().any(|t| (t.center - point).norm() < t.radius)
}

pub fn single_trap_intersect(point: na::Vector3<f64>, t: &Trap) -> bool {
    (t.center - point).norm() < t.radius
}

#[derive(Clone, Copy, Debug)]
pub struct StickyTrap {
    pub center: na::Vector3<f64>,
    pub radius: f64,
    pub desorption_time: f64,
    pub boundary_coefficient: f64,
}

impl Default for StickyTrap {
    fn default() -> Self {
        StickyTrap {
            center: na::vector![0., 0., 0.],
            radius: 1.,
            desorption_time: 1.,
            boundary_coefficient:1.,
        }
    }
}

pub fn parse_sticky_trap(stickytrap: &Vec<f64>) -> Option<StickyTrap> {
    match stickytrap.len() == 6 {
        true => Some(StickyTrap {
            radius: stickytrap[3],
            center: na::vector![stickytrap[0], stickytrap[1], stickytrap[2]],
            desorption_time: stickytrap[4],
            boundary_coefficient: stickytrap[5],
        }),
        false => None,
    }
}

pub fn stickytrap_intersect(point: na::Vector3<f64>, stickytraps: &Vec<StickyTrap>) -> bool {
    stickytraps.iter().any(|t| (t.center - point).norm() < t.radius)
}

pub fn single_stickytrap_intersect(point: na::Vector3<f64>, t: &StickyTrap) -> bool {
    (t.center - point).norm() < t.radius
}
