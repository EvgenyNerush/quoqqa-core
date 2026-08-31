// In a constant magnetic field the electron momentum change its direction
// to the opposite after a half period of the rotation. It is checked for
// different time steps here. Note that it is assumed here that the momentum
// pusher has single step error less than O((time step)^3).

use std::f64::consts::PI;
use rand::random_range;
use quoqqa_core::{Vec3, Norm, boris, vay3p, higuera_cary};
use quoqqa_core::pretty_print::{*};

fn test(pusher: fn(Vec3, f64, f64, Vec3, Vec3) -> Vec3, n_steps: i32) -> (bool, f64) {
    let b = Vec3{ x: 0.7, y: 1.1, z: 2.9 }; // just some magnetic field
    let p_initial = Vec3{ x: -b.z, y: 0.0, z: b.x }; // the initial momentum perpendicular to `b`

    let b_ampl = b.norm2().sqrt();
    let gamma_initial = (1.0 + p_initial.norm2()).sqrt();
    let period = 2.0 * PI * gamma_initial / b_ampl;
    let time_step = 0.5 * period / (n_steps as f64);
    let e = Vec3{ x: 0.0, y: 0.0, z: 0.0 }; // zero electric field

    let mut p = p_initial;
    for _i in 0..n_steps {
        p = pusher(p, -1.0, time_step, e, b);
    }

    let err = p_initial + p;
    let relative_error = (err.norm2() / p_initial.norm2()).sqrt();

    // true if the error is smaller than approximately 2 * (the estimate for error of midpoint methods)
    let estimate = (time_step * b_ampl / gamma_initial).powf(2.0);
    let is_ok = relative_error < estimate;
    (is_ok, relative_error / estimate)
}

fn main() {
    let n_tryes: i32 = 100;

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
        for _i in 0..n_tryes {
            let n_timesteps = random_range(2..200);
            let (is_ok, error) = test(pusher, n_timesteps);
            score += error;
            if !is_ok {
                println!(" test {} for n_steps = {}", red_bg("failed!"), n_timesteps);
                ok = false;
                break;
            }
        }
        if ok {
            println!(" {} score={:.5}", green_bg("passed!"), score / (n_tryes as f64));
        }
    }

    println!("All done. Bye!");
}
