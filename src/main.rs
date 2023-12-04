use std::env;

mod progress;
mod iteration;
mod utils;
mod gemm_sampler;
mod gemm_sampler_rayon;
mod gemm_sampler_spawn;
mod chunk;
mod chunk_dispatcher;
mod unsafe_utils;
mod gemm_profiler;

#[global_allocator]
static ALLOC : mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut global_pool = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();
    let method = &args[1];
    if method.eq_ignore_ascii_case("acc") {
        gemm_sampler_rayon::acc(&mut global_pool);
        println!();
        // gemm_sampler_spawn::acc();
        // println!();
        // gemm_sampler::acc();
        println!();
        // gemm_profiler::acc(128, 128, 128, 4);
    } else if method.eq_ignore_ascii_case("speed") {
        println!("RUST RAYON: ");
        for _ in 0..10{
            gemm_sampler_rayon::speed(&mut global_pool);
        }
        // println!();
        // // for _ in 0..3{
        // //     gemm_sampler_spawn::speed();
        // // }
        // // println!();
        // for _ in 0..3{
        //     gemm_sampler::speed();
        // }
        // println!();
    }
}