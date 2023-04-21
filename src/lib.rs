mod gemm_sampler;
mod chunk_dispatcher;
mod test;
mod iteration;
mod pluss_aet;
mod progress;
mod chunk;
mod utils;
use statrs::distribution::NegativeBinomial;

#[cfg(test)]
mod tests {
    use statrs::distribution::Discrete;
    use super::*;

    #[test]
    fn utils_test() {
        let r = NegativeBinomial::new(3.5, 0.5).unwrap();
        let a = r.pmf(0);
        println!("a = {}", a);
    }
}