
use std::collections::HashMap;
use std::iter::Map;
use std::ptr::null_mut;
use std::sync::{Arc, Mutex};
// use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::Duration;
use std::time::Instant;
use std::sync::mpsc::channel;
// use rgsl::rng::Rng;
use std::time::{SystemTime, UNIX_EPOCH};
#[path = "chunk_dispatcher.rs"]
mod chunk_dispatcher;
/*
 * -DTHREAD_NUM=$(TNUM) -DCHUNK_SIZE=4 -DDS=8 -DCLS=64
 */

const THREAD_NUM: usize = 4;
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
    let mut idle_threads: [i32; THREAD_NUM] = Default::default();
    let mut subscripts: Vec<i32> = Vec::new();
    // let dispatcher: ChunkDispatcher; TODO: !!!!!
    let mut tid_to_run: i32 = 0;
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
    for tid in 0..THREAD_NUM {
        idle_threads[tid] = 1;
    }

    // #pragma omp parallel for num_threads(THREAD_NUM) private(addr) TODO: !!!!! need to parallel this for loop
    // for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
    for tid in 0..THREAD_NUM {
        loop {

        }
    }
}

