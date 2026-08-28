// Magnetic and electric fields are parallel to each other. The electric field amplitude
// changes in time such that the electron momentum parallel to it is zero at the end.
// The residual "parallel" component of the electron momentum is compared with
// twice of the error caused by the terms out of the accuracy of the midpoint methods.

use std::f64::consts::PI;
use rand::random_range;
use quoqqa_core::{Vec3, Dot, vay3p, higuera_cary};
use quoqqa_core::pretty_print::{*};

fn test(pusher: fn(Vec3, f64, f64, Vec3, Vec3) -> Vec3, n_steps: i32) -> (bool, f64) {

    let e_base = Vec3{ x: 0.35, y: 0.55, z: 0.42 }; // the electric field at maximum
    let b = e_base; // just some magnetic field parallel to `e_base`
    let p_parallel = e_base;
    let p_perp = Vec3{ x: -e_base.z, y: 0.0, z: e_base.x }; // perpendicular to `e_base`
    let p_initial = p_parallel + p_perp;

    let time = 0.5 * PI;
    let time_step = time / (n_steps as f64);

    let mut p = p_initial;
    for i in 0..n_steps {
        let t = time_step * (0.5 + (i as f64));
        let e = e_base * t.cos();
        p = pusher(p, -1.0, time_step, e, b);
    }

    let err = p.dot(e_base); // analytically parallel part of `p` should be zero at `time`
    let relative_error = (err / p_initial.dot(e_base)).abs();

    let estimate = 2.0 * time_step.powi(2) * time / 6.0;
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
            let n_timesteps = random_range(2..200);
            let (is_ok, relative_error) = test(pusher, n_timesteps);
            score += relative_error;
            if !is_ok {
                println!(" test {} for n_steps = {}", red_bg("failed"), n_timesteps);
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
