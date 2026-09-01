// Tests that the gyroradius in pure const B preserves.
// Higuera-Cary pusher fails at low number of timesteps in a period.
// For smaller timestep it passes the test.

use std::f64::consts::PI;
use rand::random_range;
use quoqqa_core::{*};
use quoqqa_core::pretty_print::{*};

fn test(pusher: fn(Vec3, f64, f64, Vec3, Vec3) -> Vec3,
    b: Vec3,
    p0: Vec3,
    time: f64,
    dt: f64) -> (bool, f64) {

    let r0 = p0.cross(b) / b.norm2();

    let n_steps = (time / dt) as i32;
    let mut p = p0;
    let mut r = midpoint3r(r0, vel_e(p), -0.5 * dt);
    for _i in 0..n_steps {
        r = midpoint3r(r, vel_e(p), dt);
        let e = Vec3{ x: 0.0, y: 0.0, z: 0.0 }; // zero electric field
        p = pusher(p, -1.0, dt, e, b);
    }
    r = midpoint3r(r, vel_e(p), 0.5 * dt);

    let radius0 = r0.norm();
    let radius = r.cross(b).norm() / b.norm();

    let accuracy = 0.1;
    let relative_error = (radius - radius0).abs() / radius0;

    let is_ok = relative_error < accuracy;
    (is_ok, relative_error)
}

fn main() {
    let n_tests: i32 = 100;

    println!("Hi! Doing tests for pushers:");
    let mut ok: bool;
    let mut score: f64;

    for (pusher_name, function_name, pusher) in
    [("Boris", "boris", boris as fn(Vec3, f64, f64, Vec3, Vec3) -> Vec3),
    ("Vay", "vay3p", vay3p as fn(Vec3, f64, f64, Vec3, Vec3) -> Vec3),
    ("Higuera-Cary", "higuera_cary", higuera_cary as fn(Vec3, f64, f64, Vec3, Vec3) -> Vec3)] {
        print!("Testing {} pusher {}...", blue_bg(pusher_name), italic(function_name));
        ok = true;
        score = 0.0;
        for _i in 0..n_tests {
            let ampl = 5.0;
            let b = Vec3 { x: random_range(0.0..ampl), y: random_range(0.0..ampl), z: random_range(0.0..ampl)};
            let pampl = 4.0;
            let p0 = Vec3 { x: random_range(0.0..pampl), y: random_range(0.0..pampl), z: random_range(0.0..pampl)};
            let period = 2.0 * PI * gamma_e(p0) / b.norm();
            let n_periods = 10;
            let time = random_range(0.0..(n_periods as f64 * period));
            let dt = 2.0 * PI * 0.1;
            let (is_ok, error) = test(pusher, b, p0, time, dt);
            score += error;
            if !is_ok {
                println!(" test {} for period/dt = {:.2}; err = {:.3}", red_bg("failed!"), period / dt, error);
                ok = false;
                break;
            }
        }
        if ok {
            println!(" {} score={:}", green_bg("passed!"), score / (n_tests as f64));
        }
    }

    println!("All done. Bye!");
}
