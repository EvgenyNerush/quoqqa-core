// Tests motion in crossed E and B. Compares gyroradius with predicted value.

use std::f64::consts::PI;
use rand::random_range;
use quoqqa_core::{*};
use quoqqa_core::pretty_print::{*};

// check r_perp at n_periods + 1/2 in K'
const N_PERIODS: i32 = 10;

fn test(pusher: fn(Vec3, f64, f64, Vec3, Vec3) -> Vec3, n_steps: i32) -> (bool, f64) {

    let e = Vec3 { x: 0.0, y: 0.0, z: 0.5 };
    let b = Vec3 { x: 0.0, y: 1.0, z: 0.0 };
    let p0 = Vec3 { x: 0.0, y: 0.0, z: 0.0 }; // zero velocity in lab frame
    let r0 = Vec3 { x: 0.0, y: 0.0, z: 0.0 };

    // drifting reference frame K'
    let v_ = e.cross(b) / b.norm2();
    let gamma = 1.0 / (1.0 - v_.norm2()).sqrt();
    let b_value_ = (b - v_.cross(e)).norm() * gamma;
    let omega_ = b_value_ / gamma;
    let r_value_ = v_.norm() / omega_;

    let time = gamma * ((N_PERIODS as f64 + 0.5) * 2.0 * PI / omega_); 
    let dt = time / (n_steps as f64);

    let mut p = p0;
    let mut r = midpoint3r(r0, vel_e(p), -0.5 * dt);
    for _i in 0..n_steps {
        r = midpoint3r(r, vel_e(p), dt);
        p = pusher(p, -1.0, dt, e, b);
    }
    r = midpoint3r(r, vel_e(p), 0.5 * dt);

    let accuracy = 0.1;
    let relative_error = (r.z + 2.0 * r_value_).abs() / (2.0 * r_value_);

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
            let n_steps = random_range((N_PERIODS * 20)..(N_PERIODS * 50));
            let (is_ok, error) = test(pusher, n_steps);
            score += error;
            if !is_ok {
                println!(" test {} for n_periods = {}, n_steps = {}; err = {:.3}", red_bg("failed!"), N_PERIODS, n_steps, error);
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
