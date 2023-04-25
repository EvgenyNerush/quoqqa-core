const sqrt = @import("std").math.sqrt;
const core = @import("core.zig");
const Vec3 = core.Vec3;

// Returns momentum a time step forward, in three-dimensional space.
// Uses pusher proposed in A.V. Higuera and J.R. Cary, Physics of Plasmas 24, 052104 (2017),
// https://dx.doi.org/10.1063/1.4979989
pub fn higueraCary(
    p: Vec3, // momentum
    charge_mass_ratio: f64,
    time_step: f64,
    e: Vec3, // electric field
    b: Vec3, // magnetic field
) Vec3 {
    const q = charge_mass_ratio * time_step / 2;
    const gamma = sqrt(1 + p.x * p.x + p.y * p.y + p.z * p.z);
    const epsilon = core.mult(Vec3, q, e);
    // see Eq. (18) in Higuera17
    const u_ = core.plus(Vec3, p, epsilon);

    // see (16)
    const beta = core.mult(Vec3, q, b);
    const s = 1 + core.norm2(Vec3, u_) - core.norm2(Vec3, beta);
    const beta_u_ = core.dot(Vec3, beta, u_);
    const beta_u_u_ = core.mult(Vec3, beta_u_, u_);
    // see (20)
    const gamma_new = sqrt(0.5 * (s + sqrt(s * s + 4 * (core.norm2(Vec3, beta) + core.norm2(Vec3, beta_u_u_)))));

    const m = 1 / (1 + core.norm2(Vec3, beta) / (gamma_new * gamma_new));
    const v = core.mult(Vec3, -1 / gamma_new, core.cross(Vec3, beta, u_));
    const w = core.mult(Vec3, 1 / (gamma * gamma), beta_u_u_);
    // see (17)
    const u_new = core.mult(Vec3, m, core.plus3(Vec3, u_, v, w));

    const u_new_x_beta = core.cross(Vec3, u_new, beta);
    // see (22) and (23), t = (u_f - u_i) / 2
    const t = core.plus(Vec3, epsilon, core.mult(Vec3, 1 / gamma_new, u_new_x_beta));
    return core.plus(Vec3, p, core.mult(Vec3, 2, t)); // u_f
}
