// Electron momentum is parallel to magnetic and electric fields, which amplitude
// change in time such that the electron should stop at the end.
// The residual electron momentum is compared with twice of the error caused by the terms
// out of accuracy of midpoint methods.

const std = @import("std");
const core = @import("src/core.zig");
const Vec3 = core.Vec3;

// place here the momentum pusher you want to test
//const pusher3P = @import("src/vay.zig").vay3P;
const pusher3P = @import("src/higuera-cary.zig").higueraCary;

fn test_it(n_steps: i32) bool {
    const p_initial = Vec3{ .x = 0.35, .y = 0.55, .z = 0.42 }; // the initial momentum
    const b = p_initial; // just some magnetic field parallel to `p_initial`
    const e_base = p_initial; // the electric field at maximum

    const time = 0.5 * std.math.pi;
    const time_step = time / @intToFloat(f64, n_steps);

    var i: i32 = 0;
    var p = p_initial;
    while (i < n_steps) : (i += 1) {
        const t = time_step * (0.5 + @intToFloat(f64, i));
        const e = core.mult(Vec3, std.math.cos(t), e_base);
        p = pusher3P(p, -1, time_step, e, b);
    }

    const err = p; // because analytically `p` should be zero
    const relative_error = std.math.sqrt(core.norm2(Vec3, err) / core.norm2(Vec3, p_initial));

    const is_ok = relative_error < 2 * std.math.pow(f64, time_step, 2) * time / 6;
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
        const n_steps = rand.intRangeAtMost(i32, 2, 200);
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
