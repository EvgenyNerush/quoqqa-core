pub const Vec3 = struct { x: f64, y: f64, z: f64 };

// Dot-product for SIMD vectors of length 3
pub fn dot(a: @Vector(3, f64), b: @Vector(3, f64)) f64 {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Norm squared for SIMD vectors of length 3
pub fn norm2(a: @Vector(3, f64)) f64 {
    return dot(a, a);
}

// Cross-product for SIMD vectors of length 3
pub fn cross(a: @Vector(3, f64), b: @Vector(3, f64)) @Vector(3, f64) {
    const s = @Vector(3, f64){ a[1], a[2], a[0] };
    const t = @Vector(3, f64){ b[2], b[0], b[1] };
    const u = @Vector(3, f64){ a[2], a[0], a[1] };
    const v = @Vector(3, f64){ b[1], b[2], b[0] };
    return s * t - u * v;
}
