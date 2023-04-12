const Vec3 = @import("core.zig").Vec3;

// Returns coordinates a time step forward, in three-dimensional space.
pub fn midpoint3R(r: Vec3, v: Vec3, time_step: f64) Vec3 {
    const r1 = Vec3{ .x = r.x + v.x * time_step, .y = r.y + v.y * time_step, .z = r.z + v.z * time_step };
    return r1;
}
