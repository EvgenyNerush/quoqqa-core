// This is a reference modelling which demonstrates a single-thread performance
// of Vay's and midpointR pushers. It can be then easily compared with other
// codes/algorithms of interest.

// The rotation of a big number of particles in a constant magnetic field is simulated.
// The test should print "8.0..." as a result.
// On Xeon(R) X5550 the test takes about 175 seconds.

const std = @import("std");
const Vec3 = @import("src/core.zig").Vec3;
const vay3P = @import("src/vay.zig").vay3P;
const midpoint3R = @import("src/midpointR.zig").midpoint3R;
const print = std.debug.print;

const Particle = struct { x: f64, y: f64, z: f64, px: f64, py: f64, pz: f64 };

const p_initial = Vec3{ .x = 10, .y = 0, .z = 0 }; // initial momentum
const gamma_initial: f64 = std.math.sqrt(1 + p_initial.x * p_initial.x + p_initial.y * p_initial.y + p_initial.z * p_initial.z);
const charge_mass_ratio: f64 = -1.0; // electrons
const magnetic_field = Vec3{ .x = 0, .y = 0, .z = 10 };
const electric_field = Vec3{ .x = 0, .y = 0, .z = 0 };
const period: f64 = 2 * std.math.pi * std.math.sqrt(magnetic_field.x * magnetic_field.x +
    magnetic_field.y * magnetic_field.y +
    magnetic_field.z * magnetic_field.z) / gamma_initial;
const time: f64 = 10 * period;
const n_time_steps: i32 = 1000;
const dt: f64 = time / @intToFloat(f64, n_time_steps);

fn push(particle: Particle) Particle {
    const p0 = Vec3{ .x = particle.px, .y = particle.py, .z = particle.pz };
    const p = vay3P(p0, charge_mass_ratio, dt, electric_field, magnetic_field);

    const r0 = Vec3{ .x = particle.x, .y = particle.y, .z = particle.z };
    const gamma: f64 = std.math.sqrt(1 + p.x * p.x + p.y * p.y + p.y * p.y);
    const v = Vec3{ .x = p.x / gamma, .y = p.y / gamma, .z = p.z / gamma };
    const r = midpoint3R(r0, v, dt);
    return Particle{ .x = r.x, .y = r.y, .z = r.z, .px = p.x, .py = p.y, .pz = p.z };
}

fn push_many(num_steps: i32, initial_particle: Particle) Particle {
    var n = num_steps;
    var particle = initial_particle;
    while (n > 0) : (n -= 1) {
        particle = push(particle);
    }
    return particle;
}

pub fn main() !void {
    print("hi!\ndoing stuff...\n", .{});

    const n_particles: usize = 1_000_000;

    const allocator = std.heap.page_allocator;
    const ArrayList = std.ArrayList;
    var particles = try ArrayList(Particle).initCapacity(allocator, n_particles);
    defer particles.deinit();

    var i: usize = 0;
    while (i < n_particles) : (i += 1) {
        const x = @intToFloat(f64, i);
        try particles.append(Particle{ .x = x, .y = 0, .z = 0, .px = p_initial.x, .py = p_initial.y, .pz = p_initial.z });
    }

    const time1 = std.time.milliTimestamp();
    for (particles.items) |_, index| {
        particles.items[index] = push_many(n_time_steps, particles.items[index]);
    }
    const time2 = std.time.milliTimestamp();
    print("computation takes {d} ms\n", .{time2 - time1});

    var sum: f64 = 0.0;
    for (particles.items) |_, index| {
        sum += particles.items[index].px;
    }
    print("result = {e}\n", .{sum / @intToFloat(f64, n_particles)});
    print("bye!\n", .{});
}
