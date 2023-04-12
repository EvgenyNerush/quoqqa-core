const sqrt = @import("std").math.sqrt;
const Vec3 = @import("core.zig").Vec3;

// Returns momentum a time step forward, in three-dimensional space.
// Uses midpoint pusher proposed in J.L. Vay, Phys. Plasmas 15, 056701 (2008),
// https://dx.doi.org/10.1063/1.2837054
pub fn vay3P(
    p: Vec3, // momentum
    charge_mass_ratio: f64,
    time_step: f64,
    e: Vec3, // electric field
    b: Vec3, // magnetic field
) Vec3 {
    const gamma: f64 = sqrt(1 + p.x * p.x + p.y * p.y + p.z * p.z);
    const q: f64 = time_step * charge_mass_ratio / 2;
    const u_13 = Vec3{ // see (13) in Vay08
        .x = p.x + q * (e.x + (p.y * b.z - p.z * b.y) / gamma),
        .y = p.y + q * (e.y + (p.z * b.x - p.x * b.z) / gamma),
        .z = p.z + q * (e.z + (p.x * b.y - p.y * b.x) / gamma),
    };
    const tau = Vec3{ .x = q * b.x, .y = q * b.y, .z = q * b.z };
    const u_stroked = Vec3{ .x = u_13.x + q * e.x, .y = u_13.y + q * e.y, .z = u_13.z + q * e.z };
    const sigma: f64 = 1 +
        u_stroked.x * u_stroked.x + u_stroked.y * u_stroked.y + u_stroked.z * u_stroked.z -
        (tau.x * tau.x + tau.y * tau.y + tau.z * tau.z);
    const u_starred: f64 = u_stroked.x * tau.x + u_stroked.y * tau.y + u_stroked.z * tau.z;
    const gamma_11: f64 = // see (11)
        sqrt((sigma + sqrt(sigma * sigma + 4 *
        (tau.x * tau.x + tau.y * tau.y + tau.z * tau.z + u_starred * u_starred))) / 2);
    const t = Vec3{ .x = tau.x / gamma_11, .y = tau.y / gamma_11, .z = tau.z / gamma_11 };
    const s: f64 = 1 / (1 + t.x * t.x + t.y * t.y + t.z * t.z);
    const ut: f64 = u_stroked.x * t.x + u_stroked.y * t.y + u_stroked.z * t.z;
    const u_12 = Vec3{ // see (12)
        .x = s * (u_stroked.x + ut * t.x + (u_stroked.y * t.z - u_stroked.z * t.y)),
        .y = s * (u_stroked.y + ut * t.y + (u_stroked.z * t.x - u_stroked.x * t.z)),
        .z = s * (u_stroked.z + ut * t.z + (u_stroked.x * t.y - u_stroked.y * t.x)),
    };
    return u_12;
}
