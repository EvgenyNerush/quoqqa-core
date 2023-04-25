// This is a reference modelling which demonstrates a single-thread performance
// of Vay's and midpointR pushers. It can be then easily compared with other
// codes/algorithms of interest. Run it with the following command:
//
// zig run performance-test.zig -O ReleaseFast
//
// The test should print as a result "x = 1." (or 0.9...) and "px = 10.00..." (or "9.99...").
// On Xeon(R) X5550 the test takes about 110 seconds for Vay's pusher.

const std = @import("std");
const Vec3 = @import("src/core.zig").Vec3;
// place here the momentum pusher you want to test
//const pusher3P = @import("src/vay.zig").vay3P;
const pusher3P = @import("src/higuera-cary.zig").higueraCary;
// place here the coordinate pusher you want to test
const pusher3R = @import("src/midpointR.zig").midpoint3R;
const print = std.debug.print;

const Particle = struct { x: f64, y: f64, z: f64, px: f64, py: f64, pz: f64 };

// The rotation of a big number of particles in a constant magnetic field is simulated. The average
// coordinate and momentum after even number of periods is computed.
const p_initial = Vec3{ .x = 10, .y = 0, .z = 0 }; // initial momentum
const gamma_initial: f64 = std.math.sqrt(1 + p_initial.x * p_initial.x + p_initial.y * p_initial.y + p_initial.z * p_initial.z);
const charge_mass_ratio: f64 = -1.0; // electrons
const magnetic_field = Vec3{ .x = 0, .y = 0, .z = 10 };
const electric_field = Vec3{ .x = 0, .y = 0, .z = 0 };
const period: f64 = 2 * std.math.pi * gamma_initial /
    std.math.sqrt(magnetic_field.x * magnetic_field.x +
    magnetic_field.y * magnetic_field.y +
    magnetic_field.z * magnetic_field.z);
const time: f64 = 10 * period;
const n_time_steps: i32 = 1000;
const dt: f64 = time / @intToFloat(f64, n_time_steps);

const n_particles: usize = 1_000_000;

fn step(particle: Particle) Particle {
    const p0 = Vec3{ .x = particle.px, .y = particle.py, .z = particle.pz };
    const p = pusher3P(p0, charge_mass_ratio, dt, electric_field, magnetic_field);

    const r0 = Vec3{ .x = particle.x, .y = particle.y, .z = particle.z };
    const gamma: f64 = std.math.sqrt(1 + p.x * p.x + p.y * p.y + p.z * p.z);
    const v = Vec3{ .x = p.x / gamma, .y = p.y / gamma, .z = p.z / gamma };
    const r = pusher3R(r0, v, dt);
    return Particle{ .x = r.x, .y = r.y, .z = r.z, .px = p.x, .py = p.y, .pz = p.z };
}

fn halfStepR(timestep_sign: f64, particle: Particle) Particle {
    const r0 = Vec3{ .x = particle.x, .y = particle.y, .z = particle.z };
    const gamma: f64 = std.math.sqrt(1 + particle.px * particle.px + particle.py * particle.py + particle.pz * particle.pz);
    const v = Vec3{ .x = particle.px / gamma, .y = particle.py / gamma, .z = particle.pz / gamma };
    const r = pusher3R(r0, v, dt / 2 * timestep_sign);
    return Particle{ .x = r.x, .y = r.y, .z = r.z, .px = particle.px, .py = particle.py, .pz = particle.pz };
}

fn steps(num_steps: i32, initial_particle: Particle) Particle {
    var n: i32 = 0;
    var particle = initial_particle;
    // As we use midpoint pushers, thus for coordinates a forward half step (+dt/2) and
    // a backward half-step (-dt/2) are done before and after the main loop
    particle = halfStepR(1, particle);
    while (n < num_steps) : (n += 1) {
        particle = step(particle);
    }
    particle = halfStepR(-1, particle);
    return particle;
}

pub fn main() !void {
    print("hi!\ndoing stuff...\n", .{});

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();

    const ArrayList = std.ArrayList;
    var particles = try ArrayList(Particle).initCapacity(allocator, n_particles);
    defer particles.deinit();

    var i: usize = 0;
    while (i < n_particles) : (i += 1) {
        const x = 2 * @intToFloat(f64, i) / @intToFloat(f64, n_particles - 1);
        try particles.append(Particle{ .x = x, .y = 0, .z = 0, .px = p_initial.x, .py = p_initial.y, .pz = p_initial.z });
    }

    const time1 = std.time.milliTimestamp();
    for (particles.items) |*particle| {
        particle.* = steps(n_time_steps, particle.*);
    }
    const time2 = std.time.milliTimestamp();
    print("computation takes {d} ms\n", .{time2 - time1});

    var sum_x: f64 = 0;
    var sum_px: f64 = 0;
    for (particles.items) |particle| {
        sum_x += particle.x;
        sum_px += particle.px;
    }
    print("x = {e}\n", .{sum_x / @intToFloat(f64, n_particles)});
    print("px = {e}\n", .{sum_px / @intToFloat(f64, n_particles)});
    print("bye!\n", .{});
}
