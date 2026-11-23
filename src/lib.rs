use std::ops::{Add, AddAssign, Sub, Mul, Div};

pub mod pretty_print;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

// Dot-product
pub trait Dot {
    fn dot(self, other: Self) -> f64;
}

pub trait Norm {
    fn norm2(self) -> f64; // Norm squared
    fn norm(self) -> f64;  // Norm
}

// Cross-product
pub trait Cross {
    fn cross(self, other: Self) -> Self;
}

//// Implementations ////

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

// +=
impl AddAssign for Vec3 {
    fn add_assign(&mut self, other: Self) -> () {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
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
impl Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, scalar: f64) -> Self::Output {
        Vec3 {
            x: scalar * self.x,
            y: scalar * self.y,
            z: scalar * self.z,
        }
    }
}

// Division
impl Div<f64> for Vec3 {
    type Output = Self;

    fn div(self, scalar: f64) -> Self::Output {
        Vec3 {
            x: self.x / scalar,
            y: self.y / scalar,
            z: self.z / scalar,
        }
    }
}

impl Dot for Vec3 {
    fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

impl Norm for Vec3 {
    fn norm2(self) -> f64 {
        self.dot(self)
    }
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

//// Useful functions ////

// Electron Lorentz factor for a given electron momentum normalized to m_e c.
pub fn gamma_e<P: Norm>(p: P) -> f64 {
    (1.0 + p.norm2()).sqrt()
}

// Elecron velocity / c for a given momentum
pub fn vel_e<P: Clone + Norm + Div<f64, Output = P>>(p: P) -> P {
    p.clone() / gamma_e::<P>(p)
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
// Uses Boris pusher proposed in J.P. Boris, in Proc. 4th Conf. Numerical
// Simulation of Plasmas (NRL, Washington D.C., 1970), pp. 3-67,
// and reviewed in C.K. Birdsall and A.B. Langdon, "Plasma Physics via
// Computer Simulation" (2004), Ch. 4.
pub fn boris(
    p: Vec3, // momentum
    charge_mass_ratio: f64,
    time_step: f64,
    e: Vec3, // electric field
    b: Vec3, // magnetic field
) -> Vec3 {
    let q: f64 = time_step * charge_mass_ratio / 2.0;
    let epsilon = e * q;
    // Half electric kick: u⁻ = u_i + (qΔt/2m) E
    let u_minus = p + epsilon;
    // Rotation via magnetic field
    let gamma = gamma_e(u_minus);
    let t = b * (q / gamma);
    let s = t * (2.0 / (1.0 + t.norm2()));
    let u_prime = u_minus + u_minus.cross(t);
    let u_plus = u_minus + u_prime.cross(s);
    // Second half electric kick: u_f = u⁺ + (qΔt/2m) E
    u_plus + epsilon
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
    let gamma: f64 = gamma_e(p);
    let q: f64 = time_step * charge_mass_ratio / 2.0;
    // see (13) in Vay 2008
    let u_13: Vec3 = p + (e + p.cross(b) / gamma) * q;
    let tau: Vec3 = b * q;
    let u_stroked: Vec3 = u_13 + e * q;
    let sigma: f64 = 1.0 + u_stroked.norm2() - tau.norm2();
    let u_starred: f64 = u_stroked.dot(tau);
    // see (11)
    let gamma_11: f64 = (0.5 * (sigma + (sigma * sigma + 4.0 * (tau.norm2() + u_starred * u_starred)).sqrt())).sqrt();
    let t: Vec3 = tau / gamma_11;
    let s: f64 = 1.0 / (1.0 + t.norm2());
    let ut: f64 = u_stroked.dot(t);
    // see (12)
    (u_stroked + t * ut + u_stroked.cross(t)) * s
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
    let q: f64 = time_step * charge_mass_ratio / 2.0;
    let epsilon = e * q;
    // see Eq. (18) in Higuera 2017
    let u_ = p + epsilon;

    // see (16)
    let beta = b * q;
    let s = 1.0 + u_.norm2() - beta.norm2();
    let beta_u_ = beta.dot(u_);
    let beta_beta_u_ = beta * beta_u_;
    // see (20)
    let gamma_new = (0.5 * (s + (s * s + 4.0 * (beta.norm2() + beta_u_ * beta_u_)).sqrt())).sqrt();

    let m = 1.0 / (1.0 + beta.norm2() / (gamma_new * gamma_new));
    let v = beta.cross(u_) / (-gamma_new);
    let w = beta_beta_u_ / (gamma_new * gamma_new);
    // see (17)
    let u_new = (u_ + v + w) * m;

    let u_new_x_beta = u_new.cross(beta);
    // see (22) and (23), t = (u_f - u_i) / 2
    let t = epsilon + u_new_x_beta / gamma_new;
    p + t * 2.0 // u_f
}

// sin(x) / x, see Sec. 5 in Decyk 2023
fn sinc(x: f64) -> f64 {
    if x == 0.0 { 1.0 } else { x.sin() / x }
}

// (cos(x) - 1) / x^2, see Sec. 5 in Decyk 2023;
// uses the exact identity cos(x) - 1 = -2 sin^2(x/2)
fn cosd(x: f64) -> f64 {
    -0.5 * sinc(0.5 * x).powi(2)
}

// (sin(x) - x) / x^3, see Sec. 5 in Decyk 2023
fn sind(x: f64) -> f64 {
    if x.abs() < 0.1 {
        let x2 = x * x;
        -1.0 / 6.0 * (1.0 - x2 / 20.0 * (1.0 - x2 / 42.0 * (1.0 - x2 / 72.0)))
    } else {
        (x.sin() - x) / (x * x * x)
    }
}

// Returns the Lorentz factor gamma_a = time_step / (proper time of the step), held fixed during
// the step of the analytic Boris pusher. One Newton iteration, see (29), solves time(tau) = time_step
// with time(tau) from (55) and its derivative gamma(tau) from (54) in Decyk 2023.
fn gamma_a_newton(
    p: Vec3,            
    e_tilde_par: Vec3,  // electric field parallel to B, times charge_mass_ratio
    e_tilde_perp: Vec3, // electric field perpendicular to B, times charge_mass_ratio
    omega: Vec3,        // magnetic field times charge_mass_ratio
    time_step: f64,
) -> f64 {
    let omega_norm = omega.norm();
    let gamma_0 = gamma_e(p);
    let e_par_p = e_tilde_par.dot(p); 
    let e_par2 = e_tilde_par.norm2();
    let e_perp_p = e_tilde_perp.dot(p);
    let e_perp_p_omega = e_tilde_perp.dot(p.cross(omega)); 
    let e_perp2 = e_tilde_perp.norm2();

    // see (54)
    let gamma_of_tau = |d_tau: f64, gamma_a: f64| -> f64 {
        let e_perp_du = gamma_a * e_perp2 + e_perp_p_omega; 
        let x = omega_norm * d_tau;
        gamma_0 + e_par_p * d_tau + 0.5 * gamma_a * e_par2 * d_tau * d_tau
            + e_perp_p * sinc(x) * d_tau
            - e_perp_du * cosd(x) * d_tau * d_tau
    };
    // see (55)
    let time_of_tau = |d_tau: f64, gamma_a: f64| -> f64 {
        let e_perp_du = gamma_a * e_perp2 + e_perp_p_omega;
        let x = omega_norm * d_tau;
        gamma_0 * d_tau
            + 0.5 * e_par_p * d_tau * d_tau
            + gamma_a * e_par2 * d_tau * d_tau * d_tau / 6.0
            - e_perp_p * cosd(x) * d_tau * d_tau
            - e_perp_du * sind(x) * d_tau * d_tau * d_tau
    };

    // see (29), starting from gamma_a = gamma_0
    let d_tau_1 = time_step / gamma_0;
    let d_tau_2 = d_tau_1
        - (time_of_tau(d_tau_1, gamma_0) - time_step) / gamma_of_tau(d_tau_1, gamma_0);
    time_step / d_tau_2
}

// Returns momentum a time step forward, in three-dimensional space.
// Uses analytic Boris pusher (analyic relativistic in papers' terms, AR, version) proposed in V.K. Decyk, W.B. Mori and F. Li,
// Comput. Phys. Commun. 282, 108559 (2023),
// https://doi.org/10.1016/j.cpc.2022.108559
pub fn boris_analytic(
    p: Vec3, // momentum
    charge_mass_ratio: f64,
    time_step: f64,
    e: Vec3, // electric field
    b: Vec3, // magnetic field
) -> Vec3 {
    let omega: Vec3 = b * charge_mass_ratio;
    let omega_norm = omega.norm();

    let omega_unit: Vec3 = if omega_norm != 0.0 {
        omega / omega_norm
    } else {
        Vec3 { x: 0.0, y: 0.0, z: 0.0 }
    };
    let e_tilde = e * charge_mass_ratio;
    let e_tilde_par = omega_unit * e_tilde.dot(omega_unit); // see (15)
    let e_tilde_perp = e_tilde - e_tilde_par; // see (16)

    let gamma_a = gamma_a_newton(p, e_tilde_par, e_tilde_perp, omega, time_step);
    let phi: f64 = omega_norm * time_step / gamma_a;

    // duration of the perpendicular half kick, see (35); it is time_step / 2 for B = 0, see (53)
    let kick_time_perp = if omega_norm != 0.0 {
        gamma_a / omega_norm * (phi / 2.0).tan()
    } else {
        time_step / 2.0
    };
    let u_ekick: Vec3 = e_tilde_par * (time_step / 2.0) + e_tilde_perp * kick_time_perp;

    // see (32), the first half-step, electric kick
    let p1: Vec3 = p + u_ekick;

    // see (33), rotate the particle with the magnetic field
    let p2: Vec3 = p1 * phi.cos()
        + p1.cross(omega_unit) * phi.sin()
        + omega_unit * ((1.0 - phi.cos()) * p1.dot(omega_unit));

    // see (34), the second half-step, electric kick again
    p2 + u_ekick
}
