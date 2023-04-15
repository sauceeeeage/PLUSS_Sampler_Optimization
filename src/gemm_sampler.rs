extern "C"{
    fn reset();
    fn getCurrentChunkRange() -> String;
    fn hasNextChunk(isStatic: bool) -> bool;
}

use std::collections::HashMap;
use std::iter::Map;
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
use crate::chunk::Chunk;
use crate::progress::Progress;
use crate::iteration::Iteration;
use crate::chunk_dispatcher::chunk_dispatcher;

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

fn get_addr_C0(idx0: i64, idx1: i64) -> i64 {
    ///don't know if to change the input para to i64, i32 or usize
    let addr_C0: i64 = (idx0 * 128) + (idx1 * 1);
    addr_C0 * (DS as i64) / (CLS as i64)
}

fn get_addr_B0(idx0: i64, idx1: i64) -> i64 {
    let addr_B0: i64 = (idx0 * 128) + (idx1 * 1);
    addr_B0 * (DS as i64) / (CLS as i64)
}

fn get_addr_A0(idx0: i64, idx1: i64) -> i64 {
    let addr_A0: i64 = (idx0 * 1) + (idx1 * 128);
    addr_A0 * (DS as i64) / (CLS as i64)
}

fn get_addr_C1(idx0: i64, idx1: i64) -> i64 {
    let addr_C1: i64 = (idx0 * 128) + (idx1 * 1);
    addr_C1 * (DS as i64) / (CLS as i64)
}

fn get_addr_C2(idx0: i64, idx1: i64) -> i64 {
    let addr_C2: i64 = (idx0 * 128) + (idx1 * 1);
    addr_C2 * (DS as i64) / (CLS as i64)
}

fn get_addr_C3(idx0: i64, idx1: i64) -> i64 {
    let addr_C3: i64 = (idx0 * 128) + (idx1 * 1);
    addr_C3 * (DS as i64) / (CLS as i64)
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
    let mut LAT_A: [HashMap<u64, i64>; THREAD_NUM];
    let mut LAT_B: [HashMap<u64, i64>; THREAD_NUM];
    let mut LAT_C: [HashMap<u64, i64>; THREAD_NUM];
    //srand(time(NULL)); in c++, but seems like not used in the code

    //Generate parallel code for (c0,0,\<,128,(c0 + 1))
    //Threads are scheduled using static scheduling
    //Threads are interleaved using uniform interleaving
    loop_cnt += 1;
    // dispatcher = ChunkDispatcher(CHUNK_SIZE,((128-0)),0,1); TODO: !!!!!
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
                let mut parallel_iteration_vector: Vec<i32> = Default::default();
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

            if let Some(progress_tid) = progress[tid].as_ref() {
                if progress_tid.refs == "C0" {
                    let addr = get_addr_C0(
                        progress_tid.iteration[0] as i64,
                        progress_tid.iteration[1] as i64,
                    );

                    progress_tid.refs = "1".to_string();
                    progress_tid.refs = "1".to_string();

            
                    if let Some(reuse) = LAT_C[tid].get(&(addr as u64)) {
                        let reuse = count[tid as usize] - reuse;
                        pluss_cri_noshare_histogram_update(tid, reuse, 1);
                    }
            
                    LAT_C[tid].insert(addr as u64, count[tid as usize]);
            
                    count[tid as usize] += 1;
                    progress[tid as usize].as_mut().unwrap().increment_with_ref("C1".to_string());
                    continue;
                } // end of check to C0
            }
        }

    }

    
}

