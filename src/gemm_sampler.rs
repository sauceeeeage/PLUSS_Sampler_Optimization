// #[link(name = "pluss")]
// extern "C" {
//     fn pluss_timer_start();
//     fn pluss_timer_stop();
//     fn pluss_timer_print();
// }

use std::collections::HashMap;
use std::iter::Map;
use std::process::id;
use std::ptr::null_mut;
use std::sync::{Arc, Mutex};
// use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::channel;
use std::thread;
use std::time::Duration;
use std::time::Instant;
// use derivative::Derivative;
// use rgsl::rng::Rng;
use rayon::prelude::*;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::chunk::Chunk;
use crate::chunk_dispatcher::chunk_dispatcher;
use crate::iteration::Iteration;
use crate::progress::{Progress, self};
use crate::{pluss_aet, utils};
use tracing::{debug, error, info, instrument, span, trace, warn, Level, dispatcher};

// #[path = "utils.rs"]
// mod utils;

// #[path = "chunk_dispatcher.rs"]
// mod chunk_dispatcher;
/*
 * -DTHREAD_NUM=$(TNUM) -DCHUNK_SIZE=4 -DDS=8 -DCLS=64
 */

const THREAD_NUM: usize = 4;
///because the thread_num in gemm is 4
const CHUNK_SIZE: usize = 4;
const DS: usize = 8;
const CLS: usize = 64;

static mut max_iteration_count: i64 = 0;

fn get_addr(idx0: i64, idx1: i64) -> u64 {
    ///don't know if to change the input para to i64, i32 or usize
    let addr_C0: i64 = (idx0 * 128) + (idx1 * 1);
    (addr_C0 * (DS as i64) / (CLS as i64)) as u64
}

fn distance_to(x: u64, y: u64) -> u64 {
    if x > y {
        x - y
    } else {
        y - x
    }
}

fn update_and_clear_array(array: &mut Arc<Mutex<[HashMap<u64, i64>; THREAD_NUM]>>) {
    let mut array = array.lock().unwrap();
    for i in 0..array.len() {
        utils::pluss_cri_noshare_histogram_update(i, -1, array[i].len() as f64, None);
        array[i].clear();
    }
}

//Generate sampler without sampling
fn sampler() {
    //Declare components will be used in Parallel RI search
    ///progress is not needed in this simpler version
    // array<Progress *, THREAD_NUM> progress = { nullptr };
    let mut progress: [Arc<Mutex<Option<Progress>>>; THREAD_NUM] = Default::default();
    let mut idle_threads: Arc<Mutex<[i32; THREAD_NUM]>> = Default::default();
    let mut subscripts: Vec<i32> = Default::default();
    let mut dispatcher: Arc<Mutex<chunk_dispatcher>> =
        Arc::new(Mutex::new(chunk_dispatcher::new_with_default()));
    let mut tid_to_run: i32 = 0;
    let mut start_tid: i32 = 0;
    let mut working_threads: i32 = THREAD_NUM as i32;
    let mut addr: Arc<Mutex<u64>> = Arc::new(Mutex::new(0));
    let mut loop_cnt: u64 = 0;
    let mut count: Arc<Mutex<[i64; THREAD_NUM + 1]>> = Default::default(); //should be init as the C++ version does: {0}
    let mut LAT_A: Arc<Mutex<[HashMap<u64, i64>; THREAD_NUM]>> = Default::default();
    let mut LAT_B: Arc<Mutex<[HashMap<u64, i64>; THREAD_NUM]>> = Default::default();
    let mut LAT_C: Arc<Mutex<[HashMap<u64, i64>; THREAD_NUM]>> = Default::default();
    //srand(time(NULL)); in c++, but seems like not used in the code

    //Generate parallel code for (c0,0,\<,128,(c0 + 1))
    //Threads are scheduled using static scheduling
    //Threads are interleaved using uniform interleaving
    loop_cnt += 1;
    dispatcher
        .lock()
        .unwrap()
        .new_with_para(CHUNK_SIZE as i32, 128, 0, 1);
    // dispatcher.print();

    let mut local_idle_threads = idle_threads.lock().unwrap();
    for tid in 0..THREAD_NUM {
        local_idle_threads[tid] = 1;
    }
    std::mem::drop(local_idle_threads);

    // #pragma omp parallel for num_threads(THREAD_NUM) private(addr) TODO: !!!!! need to parallel this for loop

   /*  (0..THREAD_NUM).into_par_iter().for_each(
        // should be par_iter_mut or similar things...
        |tid: usize| {
            // println!("This is thread {}.", tid);
            let mut idle_threads = Arc::clone(&idle_threads);
            let mut dispatcher = Arc::clone(&dispatcher);
            // let mut new_progress: [Arc<Mutex<Option<Progress>>>; THREAD_NUM] = Default::default();
            // let this_progress = Arc::clone(&current_progress);
            // let mut current_progress: Arc<Mutex<Option<Progress>>> = Default::default();
            let mut current_progress: Arc<Mutex<Option<Progress>>> = Arc::clone(&progress[tid]);
            let mut addr = Arc::clone(&addr);
            let mut count = Arc::clone(&count);
            let mut LAT_C = Arc::clone(&LAT_C);
            let mut LAT_A = Arc::clone(&LAT_A);
            let mut LAT_B = Arc::clone(&LAT_B);

            fun_name(
                &mut idle_threads,
                tid,
                &mut dispatcher,
                &mut current_progress,
                &mut addr,
                &mut LAT_C,
                &mut count,
                &mut LAT_A,
                &mut LAT_B,
            );
        }
    ); */

    let mut handles = vec![];
    for tid in 0..THREAD_NUM {
        // println!("This is thread {}.xxxx", tid);
        let mut dispatcher = Arc::clone(&dispatcher);
        // let mut new_progress: [Arc<Mutex<Option<Progress>>>; THREAD_NUM] = Default::default();
        // let this_progress = Arc::clone(&current_progress);
        let mut idle_threads = Arc::clone(&idle_threads);
        let mut current_progress: Arc<Mutex<Option<Progress>>> = Default::default();
        let mut addr = Arc::clone(&addr);
        let mut count = Arc::clone(&count);
        let mut LAT_C = Arc::clone(&LAT_C);
        let mut LAT_A = Arc::clone(&LAT_A);
        let mut LAT_B = Arc::clone(&LAT_B);

        let handle = thread::spawn(move || {
            // should be par_iter_mut or similar things...
            fun_name(
                    &mut idle_threads,
                    tid,
                    &mut dispatcher,
                    &mut current_progress,
                    &mut addr,
                    &mut LAT_C,
                    &mut count,
                    &mut LAT_A,
                    &mut LAT_B,
                );
            // fun_name(&mut idle_threads, tid, &mut dispatcher, &mut progress, &mut addr, &mut LAT_C, &mut count, &mut LAT_A, &mut LAT_B);
        });

        handles.push(handle);
    } /* end of for(tid) */

    for handle in handles {
        handle.join().unwrap();
    }

    idle_threads.lock().unwrap().fill(0);

    //update and clear the lists
    update_and_clear_array(&mut LAT_C);
    update_and_clear_array(&mut LAT_A);
    update_and_clear_array(&mut LAT_B);

    //delete the progress array
    for i in 0..progress.len() {
        let mut progress_i = progress[i].lock().unwrap();
        if progress_i.is_some() {
            *progress_i = None;
        }
    }

    unsafe {
        max_iteration_count = count.lock().unwrap().iter().sum();
        // println!("count: {:?}", count);
        // println!("max iteration traversed: {}", max_iteration_count);
    }
}

fn fun_name(
    idle_threads: &mut Arc<Mutex<[i32; THREAD_NUM]>>,
    tid: usize,
    dispatcher: &mut Arc<Mutex<chunk_dispatcher>>,
    current_progress: &mut Arc<Mutex<Option<Progress>>>,
    addr: &mut Arc<Mutex<u64>>,
    LAT_C: &mut Arc<Mutex<[HashMap<u64, i64>; THREAD_NUM]>>,
    count: &mut Arc<Mutex<[i64; THREAD_NUM + 1]>>,
    LAT_A: &mut Arc<Mutex<[HashMap<u64, i64>; THREAD_NUM]>>,
    LAT_B: &mut Arc<Mutex<[HashMap<u64, i64>; THREAD_NUM]>>,
) {
    let mut idle_threads = idle_threads.lock().unwrap();
    let mut dispatcher = dispatcher.lock().unwrap();
    let mut current_progress = current_progress.lock().unwrap();
    let mut addr = addr.lock().unwrap();
    let mut LAT_C = LAT_C.lock().unwrap();
    let mut count = count.lock().unwrap();
    let mut LAT_A = LAT_A.lock().unwrap();
    let mut LAT_B = LAT_B.lock().unwrap();

    loop {
        if idle_threads[tid] == 1 && dispatcher.has_next_chunk(true) {
            let c: Chunk = dispatcher.get_next_static_chunk(tid as u32);
            // println!("c.first: {}, c.second: {}", c.first(), c.second());
            let mut parallel_iteration_vector: Vec<i32> = Vec::new();
            parallel_iteration_vector.push(c.first());
            parallel_iteration_vector.push(0);
            // println!("parallel_iteration_vector: {:?}", parallel_iteration_vector);
            if current_progress.is_some() {
                // println!("current_progress.as_mut().unwrap().refs: {}", current_progress.as_mut().unwrap().refs);
                current_progress.as_mut().unwrap().refs = "C0".parse().unwrap();
                current_progress.as_mut().unwrap().iteration = parallel_iteration_vector;
                current_progress.as_mut().unwrap().chunk = c;
            } else {
                let p = Progress::new_with_ref("C0".parse().unwrap(), parallel_iteration_vector, c);
                // println!("p.refs: {}", p.refs);
                *current_progress = Some(p);
            }
            idle_threads[tid] = 0;
        } /* end of chunk availability check */
        //UNIFORM INTERLEAVING
        if current_progress.is_none() || !current_progress.as_ref().unwrap().is_in_bound() {
            //FIXME: don't know if the second condition needs ! or not
            // println!("in the else if, and current_progress.is_none() is: {}", current_progress.is_none());
            // println!("in the else if, and current_progress.as_ref().unwrap().is_in_bound() is: {}", current_progress.as_ref().unwrap().is_in_bound());
            idle_threads[tid] = 1;
            break;
        }
        // println!("tid: {}, current_progress.as_ref().unwrap().to_string: {}",tid, current_progress.as_ref().unwrap().to_string());
        if current_progress.as_ref().unwrap().refs == "C0" {
            // println!("in the C0 if, and current_progress.as_ref().unwrap().refs is: {}", current_progress.as_ref().unwrap().refs);
            *addr = get_addr(
                current_progress.as_ref().unwrap().iteration[0] as i64,
                current_progress.as_ref().unwrap().iteration[1] as i64,
            );
            // if (LAT_C[tid_to_run].find(addr) != LAT_C[tid_to_run].end()) {
            if LAT_C[tid].contains_key(&(*addr)) {
                // long reuse = [count][tid_to_run] - LAT_C[tid_to_run][addr];
                let reuse: i64 = count[tid] - *LAT_C[tid].get(&(*addr)).unwrap();
                // println!("reuse in LAT_C: {}", reuse);
                // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
            }
            // LAT_C[tid_to_run][addr] = count[tid_to_run];
            LAT_C[tid].insert(*addr, count[tid]);
            // count[tid_to_run]++;
            count[tid] += 1;
            // println!("count[tid]: {}", count[tid]);
            current_progress
                .as_mut()
                .unwrap()
                .increment_with_ref("C1".parse().unwrap());
            continue;
        } /* end of check to C0 */

        // if (progress[tid_to_run]->ref == "C1") {
        //     addr = GetAddress_C1(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
        if current_progress.as_ref().unwrap().refs == "C1" {
            // println!("in the C1 if, and current_progress.as_ref().unwrap().refs is: {}", current_progress.as_ref().unwrap().refs);
            *addr = get_addr(
                current_progress.as_ref().unwrap().iteration[0] as i64,
                current_progress.as_ref().unwrap().iteration[1] as i64,
            );
            // if (LAT_C[tid_to_run].find(addr) != LAT_C[tid_to_run].end()) {
            if LAT_C[tid].contains_key(&(*addr)) {
                // long reuse = count[tid_to_run] - LAT_C[tid_to_run][addr];
                let reuse: i64 = count[tid] - *LAT_C[tid].get(&(*addr)).unwrap();
                // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
            }
            // LAT_C[tid_to_run][addr] = count[tid_to_run];
            LAT_C[tid].insert(*addr, count[tid]);
            // count[tid_to_run]++;
            count[tid] += 1;
            // println!("count[tid]: {}", count[tid]);
            // CASE 2
            current_progress.as_mut().unwrap().iteration.push(0);
            current_progress
                .as_mut()
                .unwrap()
                .increment_with_ref("A0".parse().unwrap());
            continue;
        } /* end of check to C1 */
        if current_progress.as_ref().unwrap().refs == "A0" {
            // println!("in the A0 if, and current_progress.as_ref().unwrap().refs is: {}", current_progress.as_ref().unwrap().refs);
            *addr = get_addr(
                current_progress.as_ref().unwrap().iteration[0] as i64,
                current_progress.as_ref().unwrap().iteration[2] as i64,
            );
            if LAT_A[tid].contains_key(&(*addr)) {
                let reuse: i64 = count[tid] - *LAT_A[tid].get(&(*addr)).unwrap();
                // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
            }
            LAT_A[tid].insert(*addr, count[tid]);
            count[tid] += 1;
            // println!("count[tid]: {}", count[tid]);
            // CASE 2
            current_progress
                .as_mut()
                .unwrap()
                .increment_with_ref("B0".parse().unwrap());
            continue;
        } /* end of check to A0 */
        if current_progress.as_ref().unwrap().refs == "B0" {
            // println!("in the B0 if, and current_progress.as_ref().unwrap().refs is: {}", current_progress.as_ref().unwrap().refs);
            *addr = get_addr(
                current_progress.as_ref().unwrap().iteration[2] as i64,
                current_progress.as_ref().unwrap().iteration[1] as i64,
            ); // don't really understand why the iteration is 2 and 1 here and also in other places
            if LAT_B[tid].contains_key(&(*addr)) {
                let reuse: i64 = count[tid] - *LAT_B[tid].get(&(*addr)).unwrap();
                /* Compare c2*/
                /* With c2*/
                /* With c1*/
                /* Compare c1*/
                /* With c2*/
                /* With c1*/
                //B[c2][c1] is carried by (c1,0,\<,128,(c1 + 1))

                // if (distance_to(reuse,0) > distance_to(reuse,(((1)*((128-0)/1)+1)*((128-0)/1)+1))) {
                if distance_to(reuse as u64, 0) > distance_to(reuse as u64, 16513) {
                    //don't really understand why 16513 is used here
                    // pluss_cri_share_histogram_update(tid_to_run,THREAD_NUM-1,reuse,1);
                    utils::pluss_cri_share_histogram_update(
                        tid as i32,
                        (THREAD_NUM - 1) as i32,
                        reuse,
                        1 as f64,
                    );
                } else {
                    // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                    utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
                }
            }
            LAT_B[tid].insert(*addr, count[tid]);
            count[tid] += 1;
            // println!("count[tid]: {}", count[tid]);
            current_progress
                .as_mut()
                .unwrap()
                .increment_with_ref("C2".parse().unwrap());
            continue;
        } /* end of check to B0 */
        if current_progress.as_ref().unwrap().refs == "C2" {
            // println!("in the C2 if, and current_progress.as_ref().unwrap().refs is: {}", current_progress.as_ref().unwrap().refs);
            *addr = get_addr(
                current_progress.as_ref().unwrap().iteration[0] as i64,
                current_progress.as_ref().unwrap().iteration[1] as i64,
            );
            if LAT_C[tid].contains_key(&(*addr)) {
                let reuse: i64 = count[tid] - *LAT_C[tid].get(&(*addr)).unwrap();
                // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
            }
            LAT_C[tid].insert(*addr, count[tid]);
            count[tid] += 1;
            // println!("count[tid]: {}", count[tid]);
            current_progress
                .as_mut()
                .unwrap()
                .increment_with_ref("C3".parse().unwrap());
            continue;
        } /* end of check to C2 */
        if current_progress.as_ref().unwrap().refs == "C3" {
            // println!("in the C3 if, and current_progress.as_ref().unwrap().refs is: {}", current_progress.as_ref().unwrap().refs);
            *addr = get_addr(
                current_progress.as_ref().unwrap().iteration[0] as i64,
                current_progress.as_ref().unwrap().iteration[1] as i64,
            );
            if LAT_C[tid].contains_key(&(*addr)) {
                let reuse: i64 = count[tid] - *LAT_C[tid].get(&(*addr)).unwrap();
                // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
            }
            LAT_C[tid].insert(*addr, count[tid]);
            count[tid] += 1;
            if (current_progress.as_ref().unwrap().iteration[2] + 1) < 128 {
                // println!("going to A0 from C3, CASE 3");
                current_progress.as_mut().unwrap().iteration[2] =
                    current_progress.as_ref().unwrap().iteration[2] + 1;
                current_progress
                    .as_mut()
                    .unwrap()
                    .increment_with_ref("A0".parse().unwrap());
                continue;
            } /* end of check to C3 */
            if (current_progress.as_ref().unwrap().iteration[1] + 1) < 128 {
                // println!("going to C0 from C3, CASE 1");
                current_progress.as_mut().unwrap().iteration[1] =
                    current_progress.as_ref().unwrap().iteration[1] + 1;
                current_progress.as_mut().unwrap().iteration.pop(); //may need to unwrap???
                current_progress
                    .as_mut()
                    .unwrap()
                    .increment_with_ref("C0".parse().unwrap());
                continue;
            } /* end of check to C3 */
            //CASE 1
            current_progress.as_mut().unwrap().iteration[0] =
                current_progress.as_ref().unwrap().iteration[0] + 1;
            if current_progress.as_ref().unwrap().is_in_bound() {
                // println!("going to C0 from C3, CASE 1");
                current_progress.as_mut().unwrap().iteration.pop();
                current_progress.as_mut().unwrap().iteration.pop();
                current_progress.as_mut().unwrap().iteration.push(0);
                current_progress
                    .as_mut()
                    .unwrap()
                    .increment_with_ref("C0".parse().unwrap());
                continue;
            } /* end of check to C3 */
        } /* end of check to C3 */
        if idle_threads[tid] == 0 {
            idle_threads[tid] = 1;
        }
        if idle_threads[tid] == 1 && !dispatcher.has_next_chunk(true) {
            // println!("breaking from the while(true) loop");
            break;
        } /* end of break condition check */
    }
    /* end of while(true) */
}

pub(crate) fn main() {
    let start = Instant::now();
    sampler();
    utils::pluss_cri_distribute(THREAD_NUM as i32);
    pluss_aet::pluss_aet();
    let end = start.elapsed();
    println!("Time elapsed in parallel is: {:?}", end);
    utils::pluss_cri_noshare_print_histogram();
    utils::pluss_cri_share_print_histogram();
    utils::pluss_print_histogram();
    utils::pluss_print_mrc();
    unsafe {
        println!("max iteration traversed: {}", max_iteration_count);
    }
}
