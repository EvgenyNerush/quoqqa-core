const Vec3 = @import("core.zig").Vec3;

// Returns coordinates a time step forward, in three-dimensional space.
pub fn midpoint3R(r: Vec3, v: Vec3, time_step: f64) Vec3 {
    const r1 = Vec3{ .x = r.x + v.x * time_step, .y = r.y + v.y * time_step, .z = r.z + v.z * time_step };
    return r1;
}

// Version which uses SIMD instructions.
// Returns coordinates a time step forward, in three-dimensional space.
pub fn simdMidpoint3R(r: Vec3, v: Vec3, time_step: f64) Vec3 {
    const _r = @Vector(3, f64){ r.x, r.y, r.z };
    const _v = @Vector(3, f64){ v.x, v.y, v.z };
    const _time_step = @splat(3, time_step);
    const _r1 = _r + _v * _time_step;
    return Vec3{ .x = _r1[0], .y = _r1[1], .z = _r1[2] };
}
