# PLUSS_Sampler_Optimization

This project aims to optimize PLUSS Sampler(https://cnc-workshop.github.io/cnc2021/slides/CnC_Liu.pdf)
using Rust with parallel techniques to improve the performance and memory safety of the original C++ GEMM sampler.

### HOW TO RUN THE CODE:
1. download and install GSL library from https://www.gnu.org/software/gsl/
2. configure the path of GSL library and g++ in the Makefile
3. in the terminal, type `sh run.sh $METHOD` to run the code
   1. `$METHOD` here will be either `acc` or `speed`
   2. `acc` will run the code with the accuracy test
   3. `speed` will run the code with the speed test
4. check the output.txt file for the result
