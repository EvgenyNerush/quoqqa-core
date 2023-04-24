// In a constant magnetic field the electron momentum change its direction
// to the opposite after a half period of the rotation. It is checked for
// different time steps here. Note that it is intended that the momentum pusher
// has single step error less than O((time step)^3).

const std = @import("std");
const core = @import("src/core.zig");
const Vec3 = core.Vec3;

// place here the momentum pusher you want to test
const pusher3P = @import("src/vay.zig").vay3P;

fn test_it(n_steps: i32) bool {
    const b = Vec3{ .x = 0.7, .y = 1.1, .z = 2.9 }; // just some magnetic field
    const p_initial = Vec3{ .x = -b.z, .y = 0, .z = b.x }; // the initial momentum; it is perpendicular to `b`

    const b_ampl = std.math.sqrt(core.norm2(Vec3, b));
    const gamma_initial = std.math.sqrt(1 + core.norm2(Vec3, p_initial));
    const period = 2 * std.math.pi * gamma_initial / b_ampl;
    const time_step = 0.5 * period / @intToFloat(f64, n_steps);
    const e = Vec3{ .x = 0, .y = 0, .z = 0 }; // electric field

    var i: i32 = 0;
    var p = p_initial;
    while (i < n_steps) : (i += 1) {
        p = pusher3P(p, -1, time_step, e, b);
    }

    const err = core.plus(Vec3, p_initial, p);
    const relative_error = std.math.sqrt(core.norm2(Vec3, err) / core.norm2(Vec3, p_initial));

    // true if the error is smaller than approximately 2 * (the estimate for error of midpoint methods)
    const is_ok = relative_error < std.math.pow(f64, time_step * b_ampl / gamma_initial, 2);
    return is_ok;
}

pub fn main() !void {
    const print = std.debug.print;

    const n_tests: i32 = 100;
    print("hi!\ndoing tests...\n", .{});

    var prng = std.rand.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.os.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    const rand = prng.random();

    var i: i32 = 0;
    while (i < n_tests) : (i += 1) {
        const n_steps = rand.intRangeAtMost(i32, 5, 200);
        if (!test_it(n_steps)) {
            print("test failed for n_steps = {d},\n", .{n_steps});
            break;
        }
    }
    if (i == n_tests) {
        print("successfully passed {d} tests,\n", .{n_tests});
    }
    print("bye!\n", .{});
}
