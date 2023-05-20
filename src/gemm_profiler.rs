use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use lazy_static::lazy_static;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

type Histogram = HashMap<i64, f64>;

const DS: i64 = 8;
const CLS: i64 = 64;
static THREAD_NUM: AtomicUsize = AtomicUsize::new(1);

lazy_static! {
    // static ref pri_array: Mutex<Vec<Histogram>> = {
    //     let mut tmp: Vec<Histogram> = Default::default();
    //     for _ in 0..THREAD_NUM.load(Ordering::Relaxed) {
    //         tmp.push(Histogram::new());
    //     }
    //     Mutex::new(tmp)
    // };
    // static ref reuse_array: Mutex<Vec<HashMap<String, HashMap<u64, u128>>>> = {
    //     let mut tmp: Vec<HashMap<String, HashMap<u64, u128>>> = Default::default();
    //     for _ in 0..THREAD_NUM.load(Ordering::Relaxed) {
    //         let mut inner: HashMap<String, HashMap<u64, u128>> = Default::default();
    //         tmp.push(inner);
    //     }
    //     Mutex::new(tmp)
    // };
    static ref pri_array: Vec<Mutex<Histogram>> = {
        let mut tmp: Vec<Mutex<Histogram>> = Default::default();
        for _ in 0..THREAD_NUM.load(Ordering::Relaxed) {
            tmp.push(Mutex::new(Histogram::new()));
        }
        tmp
    };
    static ref reuse_array: Vec<Mutex<HashMap<String, HashMap<u64, i128>>>> = {
        let mut tmp: Vec<Mutex<HashMap<String, HashMap<u64, i128>>>> = Default::default();
        for _ in 0..THREAD_NUM.load(Ordering::Relaxed) {
            let mut inner: HashMap<String, HashMap<u64, i128>> = Default::default();
            tmp.push(Mutex::new(inner));
        }
        tmp
    };
}

fn get_addr(arr_size: i64, index0: i64, index1: i64) -> u64 {
    let addr: i64 = (index0 * arr_size) + (index1 * 1);
    (addr * DS / CLS) as u64 // get the cache line address with a 64-bytes granularity
}

fn profiler(tid: usize, arr_size: usize, counter: i128, arr_name: String, arr_counter: (usize, usize)) {
    let index0 = arr_counter.0;
    let index1 = arr_counter.1;
    let addr = get_addr(arr_size as i64, index0 as i64, index1 as i64);
    let mut local_reuse_arr = reuse_array[tid].lock().unwrap(); // t1, t2, t3, ...
    let cloned_arr_name = arr_name.clone();
    let specific_named_reuse_arr = local_reuse_arr.get_mut(&arr_name); // A, B, or C
    let previous_counter: Option<i128>;
    if specific_named_reuse_arr.is_some() {
        let trace_length_hash: &mut HashMap<u64, i128> = specific_named_reuse_arr.unwrap();
        previous_counter = trace_length_hash.get(&addr).cloned();
        trace_length_hash.insert(addr, counter); // update the counter
    } else {
        let mut tmp: HashMap<u64, i128> = Default::default();
        tmp.insert(addr, counter);
        local_reuse_arr.insert(arr_name, tmp);
        previous_counter = None;
    }
    let mut reuse_distance: i128 = -1;
    let mut local_pri_arr = pri_array[tid].lock().unwrap();
    // println!("---------------------------------");
    // println!("previous_counter: {:?}", previous_counter);
    // println!("counter: {:?}, {}, {}", counter, tid, rayon::current_thread_index().unwrap());
    if previous_counter.is_some() { // += 1 for reuse distance
        let previous_counter = previous_counter.unwrap();
        reuse_distance = counter - previous_counter;
        // if counter == 52945 {
        //     println!("counter 52945: counter: {:?}, previous counter: {:?}, arr_name: {:?}, arr_counter: {:?}, reuse: {:?}, addr: {:?}, index0: {:?}, index1: {:?}", counter, previous_counter, cloned_arr_name, arr_counter, reuse_distance, addr, index0, index1);
        // }
        // if reuse_distance == 12862 {
        //     println!("reuse 12862: counter: {:?}, previous counter: {:?}, arr_name: {:?}, arr_counter: {:?}, addr: {:?}, index0: {:?}, index1: {:?}", counter, previous_counter, cloned_arr_name, arr_counter, addr, index0, index1);
        // }
    }
    if local_pri_arr.contains_key(&(reuse_distance as i64)) {
        let reuse_count = local_pri_arr.get_mut(&(reuse_distance as i64)).unwrap();
        *reuse_count += 1.0; // one more time to access this reuse distance
    } else {
        local_pri_arr.insert(reuse_distance as i64, 1.0); // first time to access this reuse distance
    }
}

fn gemm_init(
    i: usize,
    j: usize,
    k: usize,
    C: &mut Vec<Vec<f64>>,
    A: &mut Vec<Vec<f64>>,
    B: &mut Vec<Vec<f64>>)
    -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>) { // (C, A, B)
    for row in 0..i {
        let mut tmp: Vec<f64> = Default::default();
        for col in 0..j {
            tmp.push((((row * col + 1) as f64 % i as f64) / i as f64) as f64);
        }
        C.push(tmp);
    }
    for row in 0..i {
        let mut tmp: Vec<f64> = Default::default();
        for col in 0..k {
            tmp.push((((row * (col + 1)) as f64 % k as f64) / k as f64) as f64);
        }
        A.push(tmp);
    }
    for row in 0..k {
        let mut tmp: Vec<f64> = Default::default();
        for col in 0..j {
            tmp.push((((row * (col + 2)) as f64 % j as f64) / j as f64) as f64);
        }
        B.push(tmp);
    }
    (C.clone(), A.clone(), B.clone())
}

fn kernel_print(kernel: &Vec<Vec<f64>>) {
    for row in kernel {
        for col in row {
            print!("{} ", col);
        }
        println!();
    }
}

fn kernel_gemm(
    i: usize,
    j: usize,
    k: usize,
    alpha: Option<f64>,
    beta: Option<f64>,
    C: &mut Vec<Vec<f64>>,
    A: &Vec<Vec<f64>>,
    B: &Vec<Vec<f64>>)
    -> Vec<Vec<f64>> { // C
    let alpha = alpha.unwrap_or(1.5);
    let beta = beta.unwrap_or(1.2);
    let mut counter: i128 = 0;
    for c0 in 0..i { // row
        for c1 in 0..j { // col
            C[c0][c1] *= beta;
            profiler(0, j, counter, "C".to_string(), (c0, c1));
            counter += 1;
            profiler(0, j, counter, "C".to_string(), (c0, c1));
            counter += 1;
            for c2 in 0..k {
                C[c0][c1] += (alpha * A[c0][c2]) * B[c2][c1];
                profiler(0, k, counter, "A".to_string(), (c0, c2));
                counter += 1;
                profiler(0, k, counter, "B".to_string(), (c2, c1));
                counter += 1;
                profiler(0, k, counter, "C".to_string(), (c0, c1));
                counter += 1;
                profiler(0, k, counter, "C".to_string(), (c0, c1));
                counter += 1;
            }
        }
    }
    C.clone()
}

fn para_opt_gemm_kernel(
    i: usize,
    j: usize,
    k: usize,
    alpha: Option<f64>,
    beta: Option<f64>,
    C: &mut Vec<Vec<f64>>,
    A: &Vec<Vec<f64>>,
    B: &Vec<Vec<f64>>,
    thread_num: usize) {
    let alpha = alpha.unwrap_or(1.5);
    let beta = beta.unwrap_or(1.2);
    let tmp = C.clone().len();
    assert_eq!(tmp % thread_num, 0);
    // panic!("{}", tmp % thread_num);
    C.par_chunks_mut(tmp / thread_num).enumerate().for_each(|(tid, c0_chunk)| {
        let mut counter: i128 = 0;
        // println!("tid: {:?}", tid);
        for c0 in 0..c0_chunk.len() { // row
            for c1 in 0..c0_chunk[0].len() { // col
                c0_chunk[c0][c1] *= beta;
                profiler(rayon::current_thread_index().unwrap(), c0_chunk[0].len(), counter, "C".to_string(), (c0, c1));
                counter += 1;
                profiler(rayon::current_thread_index().unwrap(), c0_chunk[0].len(), counter, "C".to_string(), (c0, c1));
                counter += 1;
                for c2 in 0..k {
                    c0_chunk[c0][c1] += (alpha * A[c0][c2]) * B[c2][c1]; //((alpha * A[c0][c2]) * B[c2][c1])
                    profiler(rayon::current_thread_index().unwrap(), k, counter, "A".to_string(), (c0, c2));
                    counter += 1;
                    profiler(rayon::current_thread_index().unwrap(), k, counter, "B".to_string(), (c2, c1));
                    counter += 1;
                    profiler(rayon::current_thread_index().unwrap(), k, counter, "C".to_string(), (c0, c1));
                    counter += 1;
                    profiler(rayon::current_thread_index().unwrap(), k, counter, "C".to_string(), (c0, c1));
                    counter += 1;
                }
            }
        }
    });
}

// fn para_kernel_gemm(
//     i: usize,
//     j: usize,
//     k: usize,
//     alpha: Option<f64>,
//     beta: Option<f64>,
//     C: &mut Mutex<Vec<Vec<f64>>>,
//     A: &Vec<Vec<f64>>,
//     B: &Vec<Vec<f64>>) {
//     let alpha = alpha.unwrap_or(1.5);
//     let beta = beta.unwrap_or(1.2);
//     (0..i).into_par_iter().for_each(|c0| { // row
//         let mut counter: i128 = 0;
//         let mut C = C.lock().unwrap();
//         for c1 in 0..j { // col
//             C[c0][c1] *= beta; // C * beta => C
//             profiler(rayon::current_thread_index().unwrap(), C.len(), counter, "C".to_string(), (c0, c1));
//             counter += 1;
//             profiler(rayon::current_thread_index().unwrap(), C.len(), counter, "C".to_string(), (c0, c1));
//             counter += 1;
//             for c2 in 0..k {
//                 C[c0][c1] += (alpha * A[c0][c2]) * B[c2][c1]; // C + alpha * A * B => C
//                 profiler(rayon::current_thread_index().unwrap(), C.len(), counter, "C".to_string(), (c0, c1));
//                 counter += 1;
//                 profiler(rayon::current_thread_index().unwrap(), A.len(), counter, "A".to_string(), (c0, c2));
//                 counter += 1;
//                 profiler(rayon::current_thread_index().unwrap(), B.len(), counter, "B".to_string(), (c2, c1));
//                 counter += 1;
//                 profiler(rayon::current_thread_index().unwrap(), C.len(), counter, "C".to_string(), (c0, c1));
//                 counter += 1;
//             }
//         }
//     });
// }

fn main() {
    let mut args = std::env::args().collect::<Vec<_>>();
    if args.len() != 5 {
        println!("Usage: {} <i> <j> <k> <thread_num>", args[0]);
        std::process::exit(1);
    }
    let mut C: Vec<Vec<f64>> = Default::default();
    let mut A: Vec<Vec<f64>> = Default::default();
    let mut B: Vec<Vec<f64>> = Default::default();
    let i = args[1].parse::<usize>().unwrap();
    let j = args[2].parse::<usize>().unwrap();
    let k = args[3].parse::<usize>().unwrap();
    THREAD_NUM.store(args[4].parse::<usize>().unwrap(), Ordering::Relaxed);
    rayon::ThreadPoolBuilder::new().num_threads(THREAD_NUM.load(Ordering::Relaxed)).build_global().unwrap();
    let (mut C, A, B) = gemm_init(i, j, k, &mut C, &mut A, &mut B);
    // kernel_print(&C);
    // kernel_print(&A);
    // kernel_print(&B);
    // let seq_start = Instant::now();
    // let SEQ_C = kernel_gemm(i, j, k, Some(1.5), Some(1.2), &mut C, &A, &B);
    // let seq_end = Instant::elapsed(&seq_start);
    // println!("Sequential time: {:?}", seq_end);
    // kernel_print(&SEQ_C);
    let par_start = Instant::now();
    para_opt_gemm_kernel(i, j, k, Some(1.5), Some(1.2), &mut C, &A, &B, 4);
    let par_end = Instant::elapsed(&par_start);
    println!("Parallel time: {:?}", par_end);
    // kernel_print(&C.lock().unwrap());
    for i in 0..THREAD_NUM.load(Ordering::Relaxed) {
        println!("{:?}", pri_array[i].lock().unwrap());
    }
}

// #[test]
// fn test() {
//     let mut C: Vec<Vec<f64>> = Default::default();
//     let mut A: Vec<Vec<f64>> = Default::default();
//     let mut B: Vec<Vec<f64>> = Default::default();
//     let i = 200;
//     let j = 200;
//     let k = 200;
//     let (mut C, A, B) = gemm_init(i, j, k, &mut C, &mut A, &mut B);
//     // println!("{:?}", A);
//     // println!("{:?}", B);
//     // println!("{:?}", C);
//     let mut D = C.clone();
//     let seq_start = Instant::now();
//     kernel_gemm(i, j, k, Some(1.5), Some(1.2), &mut D, &A, &B);
//
//     let seq_end = Instant::elapsed(&seq_start);
//     println!("Sequential time: {:?}", seq_end);
//     // println!("{:?}", C);
//     let mut C = Mutex::new(C);
//     let par_start = Instant::now();
//     // para_opt_gemm_kernel(i, j, k, Some(1.5), Some(1.2), &mut C, &A, &B, 4);
//     para_kernel_gemm(i, j, k, Some(1.5), Some(1.2), &mut C, &A, &B);
//
//     let par_end = Instant::elapsed(&par_start);
//     println!("Parallel time: {:?}", par_end);
//     // println!("{:?}", D);
// }