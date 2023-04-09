mod jacobi_sampler;
mod gemm_sampler;
mod chunk_dispatcher;

extern "C" {
    // fn test() -> i32;
    fn something() -> i32;
}

fn main() {
    unsafe {
        println!("Hello, world! {}", something());
    }
}
