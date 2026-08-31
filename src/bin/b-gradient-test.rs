// Magnetic field magnitude has non-zero gradient perpendicular to the field.
// Test checks that the sum of momentum and vector potential conserves.

use std::f64::consts::PI;
use rand::random_range;
use quoqqa_core::{*};
use quoqqa_core::pretty_print::{*};

fn test(pusher: fn(Vec3, f64, f64, Vec3, Vec3) -> Vec3, time: f64, n_steps: i32) -> (bool, f64) {

    let g = 0.9; // gradient of magnetic field magnitude along the y axis
    let p_initial = Vec3{ x: 1.5, y: -0.2, z: 0.55 };
    let b0 = gamma_e(p_initial); // magnetic field at the origin, along the z axis,
                                 // such that p_perp rotation with omega approx 1

    let time_step = time / (n_steps as f64);

    let mut p = p_initial;
    let mut r = Vec3 { x: 0.0, y: 0.0, z: 0.0 };
    r = midpoint3r(r, vel_e(p), -0.5 * time_step);
    for _i in 0..n_steps {
        r = midpoint3r(r, vel_e(p), time_step);
        let e = Vec3{ x: 0.0, y: 0.0, z: 0.0 }; // zero electric field
        let b = Vec3 { x: 0.0, y: 0.0, z: b0 + g * r.y };
        p = pusher(p, -1.0, time_step, e, b);
    }
    r = midpoint3r(r, vel_e(p), 0.5 * time_step);

    let h = p.x + b0 * r.y + g * r.y.powi(2) / 2.0;
    let h0 = p_initial.x;
    let relative_error = ((h - h0) / h0).abs();

    // to estimate h, 2x estimate for p error is used: B^3/gamma^2 (Delta t)^2 t / 24
    let estimate = b0.abs() * time_step.powi(2) * time / 12.0;

    let is_ok = relative_error < estimate;
    (is_ok, relative_error / estimate)
}

fn main() {
    let n_tryes: i32 = 100;

    println!("Hi! Doing tests for pushers:");
    let mut ok: bool;
    let mut score: f64;

    for (pusher_name, function_name, pusher) in
    [("Vay", "vay3p", vay3p as fn(Vec3, f64, f64, Vec3, Vec3) -> Vec3),
    ("Higuera-Cary", "higuera_cary", higuera_cary as fn(Vec3, f64, f64, Vec3, Vec3) -> Vec3)] {
        print!("Testing {} pusher {}...", blue_bg(pusher_name), italic(function_name));
        ok = true;
        score = 0.0;
        for _i in 0..n_tryes {
            let n_timesteps = random_range(5..200);
            let time = 4.0 * PI * random_range(0.0..1.0);
            let (is_ok, error) = test(pusher, time, n_timesteps);
            score += error;
            if !is_ok {
                println!(" test {} for time = {:.2}, n_timesteps = {}; err/est = {:.5}", red_bg("failed!"), time, n_timesteps, error);
                ok = false;
                break;
            }
        }
        if ok {
            println!(" {} score={:.7}", green_bg("passed!"), score / (n_tryes as f64));
        }
    }

    println!("All done. Bye!");
}
