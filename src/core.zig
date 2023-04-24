pub const Vec3 = struct { x: f64, y: f64, z: f64 };

// Multiply vector by a number
pub fn mult(comptime T: type, a: f64, v: T) T {
    return switch (T) {
        Vec3 => Vec3{ .x = a * v.x, .y = a * v.y, .z = a * v.z },
        // SIMD vectors of length 3
        @Vector(3, f64) => {
            const s = @splat(3, a);
            s * a;
        },
        else => unreachable,
    };
}

pub fn plus(comptime T: type, a: T, b: T) T {
    return switch (T) {
        Vec3 => Vec3{ .x = a.x + b.x, .y = a.y + b.y, .z = a.z + b.z },
        @Vector(3, f64) => a + b,
        else => unreachable,
    };
}

pub fn plus3(comptime T: type, a: T, b: T, c: T) T {
    return plus(T, plus(T, a, b), c);
}

// Dot-product
pub fn dot(comptime T: type, a: T, b: T) f64 {
    return switch (T) {
        Vec3 => a.x * b.x + a.y * b.y + a.z * b.z,
        @Vector(3, f64) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
        else => unreachable,
    };
}

// Norm squared
pub fn norm2(comptime T: type, a: T) f64 {
    return dot(T, a, a);
}

// Cross-product
pub fn cross(comptime T: type, a: T, b: T) T {
    return switch (T) {
        Vec3 => Vec3{ .x = a.y * b.z - a.z * b.y, .y = a.z * b.x - a.x * b.z, .z = a.x * b.y - a.y * b.x },
        @Vector(3, f64) => {
            const s = @Vector(3, f64){ a[1], a[2], a[0] };
            const t = @Vector(3, f64){ b[2], b[0], b[1] };
            const u = @Vector(3, f64){ a[2], a[0], a[1] };
            const v = @Vector(3, f64){ b[1], b[2], b[0] };
            s * t - u * v;
        },
        else => unreachable,
    };
}
