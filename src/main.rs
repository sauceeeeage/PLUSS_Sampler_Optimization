mod progress;
mod iteration;
mod pluss_aet;
mod utils;
mod gemm_sampler;
mod gemm_sampler_rayon;
mod gemm_sampler_spawn;
mod chunk;
mod chunk_dispatcher;

fn main() {
    for _ in 0..3{
        gemm_sampler_rayon::main();
        gemm_sampler_spawn::main();
        gemm_sampler::main();
    }
}