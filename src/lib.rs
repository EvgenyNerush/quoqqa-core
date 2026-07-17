use std::ops::{Add, Sub};

pub mod pretty_print;

#[derive(Debug, Copy, Clone)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

// Plus
impl Add for Vec3 {
    type Output = Self;
    
    fn add(self, other: Self) -> Self::Output {
        Vec3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

// Minus
impl Sub for Vec3 {
    type Output = Self;
    
    fn sub(self, other: Self) -> Self::Output {
        Vec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

// Multiply vector by a number
pub trait Mult {
    fn mult(self, scalar: f64) -> Self;
}

// Dot-product
pub trait Dot {
    fn dot(self, other: Self) -> f64;
}

// Norm squared
pub trait Norm2 {
    fn norm2(self) -> f64;
}

// Norm
pub trait Norm {
    fn norm(self) -> f64;
}

// Cross-product
pub trait Cross {
    fn cross(self, other: Self) -> Self;
}

//// Implementations ////

impl Mult for Vec3 {
    fn mult(self, scalar: f64) -> Self {
        Vec3 {
            x: scalar * self.x,
            y: scalar * self.y,
            z: scalar * self.z,
        }
    }
}

impl Dot for Vec3 {
    fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

impl Norm2 for Vec3 {
    fn norm2(self) -> f64 {
        self.dot(self)
    }
}

impl Norm for Vec3 {
    fn norm(self) -> f64 {
        self.norm2().sqrt()
    }
}

impl Cross for Vec3 {
    fn cross(self, other: Self) -> Self {
        Vec3 {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
}

//// Pushers ////

// Returns coordinates a time step forward, in three-dimensional space.
pub fn midpoint3r(r: Vec3, v: Vec3, time_step: f64) -> Vec3 {
    Vec3 {
        x: r.x + v.x * time_step,
        y: r.y + v.y * time_step,
        z: r.z + v.z * time_step,
    }
}

// Returns momentum a time step forward, in three-dimensional space.
// Uses midpoint pusher proposed in J.L. Vay, Phys. Plasmas 15, 056701 (2008),
// https://dx.doi.org/10.1063/1.2837054
pub fn vay3p(
    p: Vec3, // momentum
    charge_mass_ratio: f64,
    time_step: f64,
    e: Vec3, // electric field
    b: Vec3, // magnetic field
) -> Vec3 {
    let gamma: f64 = (1.0 + p.norm2()).sqrt();
    let q: f64 = time_step * charge_mass_ratio / 2.0;
    // see (13) in Vay 2008
    let u_13: Vec3 = p + (e + p.cross(b).mult(1.0 / gamma)).mult(q);
    let tau: Vec3 = b.mult(q);
    let u_stroked: Vec3 = u_13 + e.mult(q);
    let sigma: f64 = 1.0 + u_stroked.norm2() - tau.norm2();
    let u_starred: f64 = u_stroked.dot(tau);
    // see (11)
    let gamma_11: f64 = (0.5 * (sigma + (sigma * sigma + 4.0 * (tau.norm2() + u_starred * u_starred)).sqrt())).sqrt();
    let t: Vec3 = tau.mult(1.0 / gamma_11);
    let s: f64 = 1.0 / (1.0 + t.norm2());
    let ut: f64 = u_stroked.dot(t);
    // see (12)
    (u_stroked + t.mult(ut) + u_stroked.cross(t)).mult(s)
}

// Returns momentum a time step forward, in three-dimensional space.
// Uses pusher proposed in A.V. Higuera and J.R. Cary, Physics of Plasmas 24, 052104 (2017),
// https://dx.doi.org/10.1063/1.4979989
pub fn higuera_cary(
    p: Vec3, // momentum
    charge_mass_ratio: f64,
    time_step: f64,
    e: Vec3, // electric field
    b: Vec3, // magnetic field
) -> Vec3 {
    let gamma: f64 = (1.0 + p.norm2()).sqrt();
    let q: f64 = time_step * charge_mass_ratio / 2.0;
    let epsilon = e.mult(q);
    // see Eq. (18) in Higuera 2017
    let u_ = p + epsilon;

    // see (16)
    let beta = b.mult(q);
    let s = 1.0 + u_.norm2() - beta.norm2();
    let beta_u_ = beta.dot(u_);
    let beta_u_u_ = u_.mult(beta_u_);
    // see (20)
    let gamma_new = (0.5 * (s + (s * s + 4.0 * (beta.norm2() + beta_u_ * beta_u_)).sqrt())).sqrt();

    let m = 1.0 / (1.0 + beta.norm2() / (gamma_new * gamma_new));
    let v = beta.cross(u_).mult(-1.0 / gamma_new);
    let w = beta_u_u_.mult(1.0 / (gamma * gamma));
    // see (17)
    let u_new = (u_ + v + w).mult(m);

    let u_new_x_beta = u_new.cross(beta);
    // see (22) and (23), t = (u_f - u_i) / 2
    let t = epsilon + u_new_x_beta.mult(1.0 / gamma_new);
    p + t.mult(2.0) // u_f
}
