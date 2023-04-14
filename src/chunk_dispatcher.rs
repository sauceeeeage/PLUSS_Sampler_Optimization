// // For parallel OPT
// #include <shared_mutex>
// #include <atomic>
// // For CRI
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_cdf.h> // cdf of gaussian distribution

use std::collections::HashMap;
use std::thread;
use std::sync::mpsc::channel;
use std::sync::mpsc::{Sender, Receiver};
use rgsl::rng::Rng;
use std::time::{SystemTime, UNIX_EPOCH};
use std::time::Instant;
use std::time::Duration;
use std::sync::{Arc, Mutex};
use std::thread::Thread;
// use derivative::Derivative;

type Histogram = HashMap<i64, f64>; /// maybe can be replaced by a struct
const THREAD_NUM: usize = 4; ///because the thread_num in gemm is 4

pub(crate) struct Chunk {
    lb: i32,
    ub: i32,
}

impl Chunk {
    pub fn new(lb: i32, ub: i32) -> Chunk {
        Chunk {
            lb,
            ub,
        }
    }
    pub fn first(&self) -> i32 {
        self.lb
    }
    pub fn second(&self) -> i32 {
        self.ub
    }
}

// #[derive(Derivative)]
// #[derivative(Debug, Default)]
pub(crate) struct chunk_dispatcher{
    // int lb = 0;
    // int ub = 0;
    // int chunk_size = 0;
    // int trip = 0;
    // int start = 0;
    // int last = 0;
    // int step = 1;
    // int avail_chunk = 0;
    // array<int, THREAD_NUM> per_thread_start_point; // use for static sheduling only
    // TODO: maybe it's not i32 but i64
    lb: i32,
    ub: i32,
    chunk_size: i32,
    trip: i32,
    start: i32,
    last: i32,
    step: i32, // FIXME: step start from 1
    avail_chunk: i32,
    per_thread_start_point: [i32; THREAD_NUM],
}

trait chunk_dispatcher_trait {
    // ChunkDispatcher(int chunk_size, int trip, int start_point=0, int step=1) {
    fn new(chunk_size: i32, trip: i32, start_point: i32, step: i32) -> Self;
    fn new_with_default() -> Self;
    fn init() -> Self;
    fn reset(&self) -> Self;
    fn get_current_chunk_range(&self) -> String;
    fn has_next_chunk(&self, is_static: bool) -> bool;
    fn get_current_static_chunk_range(&self, tid: u32) -> String;
    // bool hasNextStaticChunk(unsigned tid) {
    fn has_next_static_chunk(&self, tid: u32) -> bool;
    // Chunk getNextChunk(unsigned tid) {
    fn get_next_chunk(&self, tid: u32) -> Chunk;
    // Chunk getNextStaticChunk(unsigned tid) {
    fn get_next_static_chunk(&self, tid: u32) -> Chunk;
    // unsigned getStaticTid(int i) {
    fn get_static_tid(&self, i: i32) -> u32;
    // int getStaticChunkID(int i) {
    fn get_static_chunk_id(&self, i: i32) -> i32;
    // int getStaticThreadLocalPos(int i) {
    fn get_static_thread_local_pos(&self, i: i32) -> i32;
    // void setStartPoint(int i) {
    fn set_start_point(&self, i: i32);
    // Chunk getStaticStartChunk(int i, int tid) {
    fn get_static_start_chunk(&self, i: i32, tid: u32) -> Chunk;
    // Chunk getStartChunk(int i) {
    fn get_start_chunk(&self, i: i32) -> Chunk;
    // void getNextKChunksFrom(int k, Chunk c, vector<Chunk> &next_chunks) {
    fn get_next_k_chunks_from(&self, k: i32, c: Chunk, next_chunks: Vec<Chunk>);
    // void getPrevKChunksFrom(int k, Chunk c, vector<Chunk> &prev_chunks) {
    fn get_prev_k_chunks_from(&self, k: i32, c: Chunk, prev_chunks: Vec<Chunk>);
    // void moveToNextChunk() {
    fn move_to_next_chunk(&self);
    // void moveToNextStaticChunk(int tid) {
    fn move_to_next_static_chunk(&self, tid: u32);
}

impl chunk_dispatcher {

    fn init(&mut self) {
        // void init() {
        // // the range of the parallel loop will in [start, last]
        // this->avail_chunk = (this->trip % this->chunk_size) == 0 ? this->trip / this->chunk_size : this->trip / this->chunk_size + 1;
        //// the lb and ub of the first chunk (dynamic scheduling only)
        // 		if (this->step > 0) {
        // 			this->lb = this->start;
        // 			this->ub = (this->lb + (this->chunk_size-1)*this->step) <= this->last ? (this->lb + (this->chunk_size-1)*this->step) : this->last;
        // 		} else {
        // 			this->ub = this->start;
        // 			this->lb = (this->ub + (this->chunk_size-1)*this->step) >= this->trip ? (this->ub + (this->chunk_size-1)*this->step) : this->last;
        // 		}
        // // assign the start chunk of each thread (static scheduling only)
        // for (int t = 0; t < THREAD_NUM; t++) {
        // this->per_thread_start_point[t] = (this->start + (this->chunk_size*this->step) * t);
        // }
        self.avail_chunk = if self.trip % self.chunk_size == 0 {
            self.trip / self.chunk_size
        } else {
            self.trip / self.chunk_size + 1
        };

        if self.step > 0 {
            self.lb = self.start;
            self.ub = if self.lb + (self.chunk_size - 1) * self.step <= self.last {
                self.lb + (self.chunk_size - 1) * self.step
            } else {
                self.last
            };
        } else {
            self.ub = self.start;
            self.lb = if self.ub + (self.chunk_size - 1) * self.step >= self.trip {
                self.ub + (self.chunk_size - 1) * self.step
            } else {
                self.last
            };
        }

        let mut per_thread_start_point = [0; THREAD_NUM];
        for t in 0..THREAD_NUM {
            self.per_thread_start_point[t] = self.start + (self.chunk_size * self.step) * (t as i32);
        }
    }

    pub fn new(&mut self, chunk_size: i32, trip: i32, start_point: i32, step: i32) {
        // default para for start_point is 0 and step is 1

        // ChunkDispatcher(int chunk_size, int trip, int start_point=0, int step=1) {
        // it is possible that the chunk_size is greater than the trip count
        // this->chunk_size = chunk_size;
        // this->trip = trip;
        // this->start = start_point;
        // this->step = step;
        // this->last = start_point + (trip - 1) * step; // the last iteration
        //
        // init();
        // }

        chunk_dispatcher {
            lb: 0,
            ub: 0,
            chunk_size,
            trip,
            start: start_point,
            last: start_point + (trip - 1) * step,
            step,
            avail_chunk: 0,
            per_thread_start_point: [0; THREAD_NUM],
        };
        self.init(); //FIXME: not sure if this is correct, if put it before the assignment
    }

    pub fn new_with_default() -> Self {
        chunk_dispatcher {
            lb: 0,
            ub: 0,
            chunk_size: 0,
            trip: 0,
            start: 0,
            last: 0,
            step: 1,
            avail_chunk: 0,
            per_thread_start_point: [0; THREAD_NUM],
        }
    }

    pub fn has_next_chunk(&self, is_static: bool) -> bool {
        // for dynamic, the next chunk is available if there avail_chunk is not 0
        // bool hasNextChunk(bool isStatic) {
        //     bool ret = false;
        //     if (isStatic) {
        //         // static
        //         // hasNextChunk return false if hasNextStaticChunk(tid) return false for all tid
        //         for (int tid = 0; tid < THREAD_NUM; tid++) {
        //             ret |= hasNextStaticChunk(tid);
        //         }
        //         return ret;
        //     }
        //     // dynamic
        //     // the lb of the next chunk exceed the parallel loop bound
        //     if (step > 0)
        //     return this->lb <= this->last;
        //     else
        //     return this->ub >= this->last;
        // }
        return if is_static {
            for tid in 0..THREAD_NUM {
                if self.has_next_static_chunk(tid as u32) {
                    return true;
                }
            }
            false
        } else {
            if self.step > 0 {
                self.lb <= self.last
            } else {
                self.ub >= self.last
            }
        }
    }

    pub fn has_next_static_chunk(&self, tid: u32) -> bool {
        // for static, the next chunk is available for a given tid if its next
        // available chunk still in the loop range
        // bool hasNextStaticChunk(unsigned tid) {
        //     if (step > 0)
        //     return this->per_thread_start_point[tid] <= this->last;
        //     else
        //     return this->per_thread_start_point[tid] >= this->last;
        // }
        if self.step > 0 {
            self.per_thread_start_point[tid as usize] <= self.last
        } else {
            self.per_thread_start_point[tid as usize] >= self.last
        }
    }

    pub fn get_next_static_chunk(&mut self, tid: u32) -> Chunk {
        // for static
        // Chunk getNextStaticChunk(unsigned tid) {
        //     int retlb = 0, retub = 0;
        //     if (step > 0) {
        //         retlb = this->per_thread_start_point[tid];
        //         retub = (retlb + (chunk_size-1)*step) < this->last ? (retlb + (chunk_size-1)*step) : this->last;
        //     } else {
        //         retub = this->per_thread_start_point[tid];
        //         retlb = (retub + (chunk_size-1)*step) > this->last ? (retub + (chunk_size-1)*step) : this->last;
        //     }
        //     Chunk curr = make_pair(retlb, retub);
        //     this->per_thread_start_point[tid] += (chunk_size*THREAD_NUM)*step;
        //     return curr;
        // }
        let mut retlb: i32 = 0;
        let mut retub: i32 = 0;
        if self.step > 0 {
            retlb = self.per_thread_start_point[tid as usize];
            retub = if retlb + (self.chunk_size - 1) * self.step < self.last {
                retlb + (self.chunk_size - 1) * self.step
            } else {
                self.last
            };
        } else {
            retub = self.per_thread_start_point[tid as usize];
            retlb = if retub + (self.chunk_size - 1) * self.step > self.last {
                retub + (self.chunk_size - 1) * self.step
            } else {
                self.last
            };
        }
        let curr = Chunk::new(retlb, retub);
        self.per_thread_start_point[tid as usize] += self.chunk_size * THREAD_NUM as i32 * self.step;
        curr
    }
    
}