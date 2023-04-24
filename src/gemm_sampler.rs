extern "C"{
    fn reset();
    fn getCurrentChunkRange() -> String;
    fn hasNextChunk(isStatic: bool) -> bool;
}

#[link(name = "pluss")]
extern "C" {
    fn pluss_timer_start();
    fn pluss_timer_stop();
    fn pluss_timer_print();
}


use std::collections::HashMap;
use std::iter::Map;
use std::process::id;
use std::ptr::null_mut;
use std::sync::{Arc, Mutex};
// use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::Duration;
use std::time::Instant;
use std::sync::mpsc::channel;
// use derivative::Derivative;
// use rgsl::rng::Rng;
use std::time::{SystemTime, UNIX_EPOCH};
use rayon::prelude::*;

use crate::chunk::Chunk;
use crate::progress::Progress;
use crate::iteration::Iteration;
use crate::chunk_dispatcher::chunk_dispatcher;
use crate::{pluss_aet, utils};

// #[path = "utils.rs"]
// mod utils;

// #[path = "chunk_dispatcher.rs"]
// mod chunk_dispatcher;
/*
 * -DTHREAD_NUM=$(TNUM) -DCHUNK_SIZE=4 -DDS=8 -DCLS=64
 */

const THREAD_NUM: usize = 4; ///because the thread_num in gemm is 4
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

fn update_and_clear_array(array: &mut [HashMap<u64, i64>]) {
    for i in 0..array.len() {
        // pluss_cri_noshare_histogram_update(i, -1, LAT_C[i].len());
        array[i].clear();
    }
}

//Generate sampler without sampling
fn sampler() {
    //Declare components will be used in Parallel RI search
    ///progress is not needed in this simpler version
    // array<Progress *, THREAD_NUM> progress = { nullptr };
    let mut progress: [Option<Progress>; THREAD_NUM] = Default::default();
    let mut idle_threads: [i32; THREAD_NUM] = Default::default();
    let mut subscripts: Vec<i32> = Default::default();
    // let mut dispatcher = chunk_dispatcher::chunk_dispatcher::new_with_default();
    let mut dispatcher = chunk_dispatcher::new_with_default();
    let mut start_tid: i32 = 0;
    let mut working_threads: i32 = THREAD_NUM as i32;
    let mut addr: u64 = 0;
    let mut loop_cnt: u64 = 0;
    let mut count: [i64; THREAD_NUM + 1] = Default::default(); //should be init as the C++ version does: {0}
    let mut LAT_A: [HashMap<u64, i64>; THREAD_NUM] = Default::default();
    let mut LAT_B: [HashMap<u64, i64>; THREAD_NUM] = Default::default();
    let mut LAT_C: [HashMap<u64, i64>; THREAD_NUM] = Default::default();
    //srand(time(NULL)); in c++, but seems like not used in the code

    //Generate parallel code for (c0,0,\<,128,(c0 + 1))
    //Threads are scheduled using static scheduling
    //Threads are interleaved using uniform interleaving
    loop_cnt += 1;
    dispatcher.new( CHUNK_SIZE as i32, 128, 0, 1);


    for tid in 0..THREAD_NUM {
        idle_threads[tid] = 1;
    }

    // #pragma omp parallel for num_threads(THREAD_NUM) private(addr) TODO: !!!!! need to parallel this for loop
    // for (tid = 0; tid < THREAD_NUM; tid++) {
    for tid in 0..THREAD_NUM {
        loop {
            // if (idle_threads[tid] == 1 && dispatcher.hasNextChunk(1)) {
            if idle_threads[tid] == 1 && dispatcher.has_next_chunk(true) {
                // Chunk c = dispatcher.getNextStaticChunk(tid);
                // vector<int> parallel_iteration_vector;
                // parallel_iteration_vector.push_back(c.first);
                // parallel_iteration_vector.push_back(0 );
                // if (progress[tid]) {
                //     progress[tid]->ref = "C0";
                //     progress[tid]->iteration = parallel_iteration_vector;
                //     progress[tid]->chunk = c;
                // } else {
                //     Progress *p = new Progress("C0", parallel_iteration_vector, c);
                //     progress[tid] = p;
                // }
                // idle_threads[tid] = 0;
                let c: Chunk = dispatcher.get_next_static_chunk(tid as u32);
                let mut parallel_iteration_vector: Vec<i32> = Vec::new();
                parallel_iteration_vector.push(c.first());
                parallel_iteration_vector.push(0);
                if progress[tid].is_some() {
                    // *progress[tid].unwrap().refs = "C0".parse().unwrap();
                    // *progress[tid].unwrap().iteration = parallel_iteration_vector;
                    // let mut a = (progress[tid].unwrap()).as_ref();
                    // a.chunk = c;
                    // let mut p: &'static Progress = Box::leak(progress[tid].unwrap().to_owned().into());
                    progress[tid].as_mut().unwrap().refs = "C0".parse().unwrap();
                    progress[tid].as_mut().unwrap().iteration = parallel_iteration_vector;
                    progress[tid].as_mut().unwrap().chunk = c;
                } else {
                    let p = Progress::new_with_ref("C0".parse().unwrap(), parallel_iteration_vector, c);
                    progress[tid] = Some(p);
                }
                idle_threads[tid] = 0;
            } /* end of chunk availability check */
            //UNIFORM INTERLEAVING
            // if (!progress[tid] || !progress[tid]->isInBound()) {
            //     idle_threads[tid] = 1;
                // break;
            // }
            if progress[tid].is_none() || !progress[tid].as_ref().unwrap().is_in_bound() {
                idle_threads[tid] = 1;
                break;
            }
            // if (progress[tid_to_run]->ref == "C0") {
            // addr = GetAddress_C0(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
            if progress[tid].as_ref().unwrap().refs == "C0" {
                let addr = get_addr(
                    progress[tid].as_ref().unwrap().iteration[0] as i64,
                    progress[tid].as_ref().unwrap().iteration[1] as i64,
                );
                // if (LAT_C[tid_to_run].find(addr) != LAT_C[tid_to_run].end()) {
                if LAT_C[tid].contains_key(&(addr)) {
                    // long reuse = count[tid_to_run] - LAT_C[tid_to_run][addr];
                    let reuse: i64 = count[tid] - LAT_C[tid].get(&(addr)).unwrap();
                    // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                    utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
                }
                // LAT_C[tid_to_run][addr] = count[tid_to_run];
                LAT_C[tid].insert(addr, count[tid]);
                // count[tid_to_run]++;
                count[tid] += 1;
                progress[tid].as_mut().unwrap().increment_with_ref("C1".parse().unwrap());
                continue;
            } /* end of check to C0 */

            // if (progress[tid_to_run]->ref == "C1") {
            //     addr = GetAddress_C1(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
            if progress[tid].as_ref().unwrap().refs == "C1" {
                let addr = get_addr(
                    progress[tid].as_ref().unwrap().iteration[0] as i64,
                    progress[tid].as_ref().unwrap().iteration[1] as i64,
                );
                // if (LAT_C[tid_to_run].find(addr) != LAT_C[tid_to_run].end()) {
                if LAT_C[tid].contains_key(&(addr)) {
                    // long reuse = count[tid_to_run] - LAT_C[tid_to_run][addr];
                    let reuse: i64 = count[tid] - LAT_C[tid].get(&(addr)).unwrap();
                    // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                    utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
                }
                // LAT_C[tid_to_run][addr] = count[tid_to_run];
                LAT_C[tid].insert(addr, count[tid]);
                // count[tid_to_run]++;
                count[tid] += 1;
                // CASE 2
                progress[tid].as_mut().unwrap().iteration.push(0);
                progress[tid].as_mut().unwrap().increment_with_ref("A0".parse().unwrap());
                continue;
            } /* end of check to C1 */
            if progress[tid].as_ref().unwrap().refs == "A0" {
                let addr = get_addr(
                    progress[tid].as_ref().unwrap().iteration[0] as i64,
                    progress[tid].as_ref().unwrap().iteration[2] as i64,
                );
                if LAT_A[tid].contains_key(&(addr)) {
                    let reuse: i64 = count[tid] - LAT_A[tid].get(&(addr)).unwrap();
                    // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                    utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
                }
                LAT_A[tid].insert(addr, count[tid]);
                count[tid] += 1;
                // CASE 2
                progress[tid].as_mut().unwrap().increment_with_ref("B0".parse().unwrap());
                continue;
            } /* end of check to A0 */
            if progress[tid].as_ref().unwrap().refs == "B0" {
                let addr = get_addr(
                    progress[tid].as_ref().unwrap().iteration[2] as i64,
                    progress[tid].as_ref().unwrap().iteration[1] as i64,
                ); // don't really understand why the iteration is 2 and 1 here and also in other places
                if LAT_B[tid].contains_key(&(addr)) {
                    let reuse: i64 = count[tid] - LAT_B[tid].get(&(addr)).unwrap();
                    /* Compare c2*/
                    /* With c2*/
                    /* With c1*/
                    /* Compare c1*/
                    /* With c2*/
                    /* With c1*/
                    //B[c2][c1] is carried by (c1,0,\<,128,(c1 + 1))

                    // if (distance_to(reuse,0) > distance_to(reuse,(((1)*((128-0)/1)+1)*((128-0)/1)+1))) {
                    if distance_to(reuse as u64, 0) > distance_to(reuse as u64, 16513) { //don't really understand why 16513 is used here
                        // pluss_cri_share_histogram_update(tid_to_run,THREAD_NUM-1,reuse,1);
                        // utils::pluss_cri_share_histogram_update(tid, THREAD_NUM - 1, reuse, 1); TODO: !!!!!!!!!!!!
                    } else {
                        // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                        utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
                    }
                }
                LAT_B[tid].insert(addr, count[tid]);
                count[tid] += 1;
                progress[tid].as_mut().unwrap().increment_with_ref("C2".parse().unwrap());
                continue;
            } /* end of check to B0 */
            if progress[tid].as_ref().unwrap().refs == "C2" {
                let addr = get_addr(
                    progress[tid].as_ref().unwrap().iteration[0] as i64,
                    progress[tid].as_ref().unwrap().iteration[1] as i64,
                );
                if LAT_C[tid].contains_key(&(addr)) {
                    let reuse: i64 = count[tid] - LAT_C[tid].get(&(addr)).unwrap();
                    // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                    utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None);
                }
                LAT_C[tid].insert(addr, count[tid]);
                count[tid] += 1;
                progress[tid].as_mut().unwrap().increment_with_ref("C3".parse().unwrap());
                continue;
            } /* end of check to C2 */
            if progress[tid].as_ref().unwrap().refs == "C3" {
                let addr = get_addr(
                    progress[tid].as_ref().unwrap().iteration[0] as i64,
                    progress[tid].as_ref().unwrap().iteration[1] as i64,
                );
                if LAT_C[tid].contains_key(&(addr)) {
                    let reuse: i64 = count[tid] - LAT_C[tid].get(&(addr)).unwrap();
                    // pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                    utils::pluss_cri_noshare_histogram_update(tid, reuse, 1 as f64, None); 
                }
                // //CASE 3
                // if ((progress[tid_to_run]->iteration[2] + 1)<128) {
                //     progress[tid_to_run]->iteration[2] = (progress[tid_to_run]->iteration[2] + 1);
                //     progress[tid_to_run]->increment("A0");
                //     continue;
                // } /* end of check to C3 */
                // //CASE 1
                // if ((progress[tid_to_run]->iteration[1] + 1)<128) {
                //     progress[tid_to_run]->iteration[1] = (progress[tid_to_run]->iteration[1] + 1);
                //     progress[tid_to_run]->iteration.pop_back();
                //     progress[tid_to_run]->increment("C0");
                //     continue;
                // } /* end of check to C3 */
                // //CASE 1
                // progress[tid_to_run]->iteration[0] = (progress[tid_to_run]->iteration[0] + 1);
                // if (progress[tid_to_run]->isInBound()) {
                //     progress[tid_to_run]->iteration.pop_back();
                //     progress[tid_to_run]->iteration.pop_back();
                //     progress[tid_to_run]->iteration.emplace_back(0);
                //     progress[tid_to_run]->increment("C0");
                //     continue;
                // } /* end of check to C3 */
                if (progress[tid].as_ref().unwrap().iteration[2] + 1) < 128 {
                    progress[tid].as_mut().unwrap().iteration[2] = progress[tid].as_ref().unwrap().iteration[2] + 1;
                    progress[tid].as_mut().unwrap().increment_with_ref("A0".parse().unwrap());
                    continue;
                } /* end of check to C3 */
                if (progress[tid].as_ref().unwrap().iteration[1] + 1) < 128 {
                    progress[tid].as_mut().unwrap().iteration[1] = progress[tid].as_ref().unwrap().iteration[1] + 1;
                    progress[tid].as_mut().unwrap().iteration.pop(); //may need to unwrap???
                    progress[tid].as_mut().unwrap().increment_with_ref("C0".parse().unwrap());
                    continue;
                } /* end of check to C3 */
                //CASE 1
                progress[tid].as_mut().unwrap().iteration[0] = progress[tid].as_ref().unwrap().iteration[0] + 1;
                if progress[tid].as_ref().unwrap().is_in_bound() {
                    progress[tid].as_mut().unwrap().iteration.pop();
                    progress[tid].as_mut().unwrap().iteration.pop();
                    progress[tid].as_mut().unwrap().iteration.push(0);
                    progress[tid].as_mut().unwrap().increment_with_ref("C0".parse().unwrap());
                    continue;
                } /* end of check to C3 */
            } /* end of check to C3 */
            if idle_threads[tid] == 0 {
                idle_threads[tid] = 1;
            }
            if idle_threads[tid] == 1 && !dispatcher.has_next_chunk(true){
                break;
            } /* end of break condition check */
        } /* end of while(true) */
    } /* end of for(tid) */

    // //reset both lists so they can be reused for later parallel loop
    // idle_threads.fill(0);
    // //Addresses in C with no data reuse will be marked as -1
    // for (unsigned i = 0; i < LAT_C.size(); i++) {
    //     pluss_cri_noshare_histogram_update(i,-1,LAT_C[i].size());
    //     LAT_C[i].clear();
    // }
    // //Addresses in A with no data reuse will be marked as -1
    // for (unsigned i = 0; i < LAT_A.size(); i++) {
    //     pluss_cri_noshare_histogram_update(i,-1,LAT_A[i].size());
    //     LAT_A[i].clear();
    // }
    // //Addresses in B with no data reuse will be marked as -1
    // for (unsigned i = 0; i < LAT_B.size(); i++) {
    //     pluss_cri_noshare_histogram_update(i,-1,LAT_B[i].size());
    //     LAT_B[i].clear();
    // }
    // for (unsigned i = 0; i < progress.size(); i++) {
    //     if (progress[i]) {
    //         delete progress[i];
    //         progress[i] = nullptr;
    //     }
    // }

    idle_threads.fill(0);

    //update and clear the lists
    update_and_clear_array(&mut LAT_C);
    update_and_clear_array(&mut LAT_A);
    update_and_clear_array(&mut LAT_B);

    //delete the progress array
    for i in 0..progress.len() {
        if progress[i].is_some() {
            progress[i] = None;
        }
    }
    unsafe {
        max_iteration_count = count.iter().sum();
    }
}

fn main(){
    unsafe{
        pluss_timer_start();
    }
    sampler();
    utils::pluss_cri_distribute(THREAD_NUM as i32);
    pluss_aet::pluss_aet();
    unsafe {
        pluss_timer_stop();
        pluss_timer_print();
    }
    utils::pluss_cri_noshare_print_histogram();
    utils::pluss_cri_share_print_histogram();
    utils::pluss_print_histogram();
    utils::pluss_print_mrc();
    unsafe {
        println!("max iteration traversed: {}", max_iteration_count);
    }
}

