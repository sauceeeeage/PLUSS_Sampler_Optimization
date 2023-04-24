mod progress;
mod iteration;
mod pluss_aet;
mod utils;
mod gemm_sampler;
mod chunk;
mod chunk_dispatcher;

fn main() {
    gemm_sampler::main();
}