mod jacobi_sampler;
mod gemm_sampler;
mod chunk_dispatcher;

fn main() {
    let mut progress: [i32; 4] = Default::default();
    assert_eq!(progress, [0, 0, 0, 0]);
}
