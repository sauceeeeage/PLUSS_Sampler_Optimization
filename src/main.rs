use std::env;

mod progress;
mod iteration;
mod utils;
mod gemm_sampler;
mod gemm_sampler_rayon;
mod gemm_sampler_spawn;
mod chunk;
mod chunk_dispatcher;

fn main() {
    let args: Vec<String> = env::args().collect();

    let method = &args[1];
    if method.eq_ignore_ascii_case("acc") {
        gemm_sampler_rayon::acc();
        println!();
        gemm_sampler_spawn::acc();
        println!();
        gemm_sampler::acc();
        println!();
    } else if method.eq_ignore_ascii_case("speed") {
        for _ in 0..5{
            gemm_sampler_rayon::speed();
        }
        println!();
        for _ in 0..5{
            gemm_sampler_spawn::speed();
        }
        println!();
        for _ in 0..5{
            gemm_sampler::speed();
        }
        println!();
    }
}