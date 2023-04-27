use std::collections::HashMap;
use std::thread;
use std::sync::mpsc::channel;
use std::sync::mpsc::{Sender, Receiver};
use std::time::{SystemTime, UNIX_EPOCH};
use std::time::Instant;
use std::time::Duration;
use std::sync::{Arc, Mutex};
use std::thread::Thread;
use crate::chunk::Chunk;
// use derivative::Derivative;

type Histogram = HashMap<i64, f64>; /// maybe can be replaced by a struct
const THREAD_NUM: usize = 4; ///because the thread_num in gemm is 4

// #[derive(Derivative)]
// #[derivative(Debug, Default)]
pub(crate) struct chunk_dispatcher{
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
        if self.chunk_size != 0 {
            self.avail_chunk = if self.trip % self.chunk_size == 0 {
                self.trip / self.chunk_size
            } else {
                self.trip / self.chunk_size + 1
            };
        }

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

        for t in 0..THREAD_NUM {
            self.per_thread_start_point[t] = self.start + (self.chunk_size * self.step) * (t as i32);
        }
    }

    pub fn new(&mut self, chunk_size: i32, trip: i32, start_point: i32, step: i32) {
        // default para for start_point is 0 and step is 1

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

    pub fn new_with_para(&mut self, chunk_size: i32, trip: i32, start_point: i32, step: i32) {
        // // it is possible that the chunk_size is greater than the trip count
        self.chunk_size = chunk_size;
        self.trip = trip;
        self.start = start_point;
        self.step = step;
        self.last = start_point + (trip - 1) * step; // the last iteration
        self.avail_chunk = if trip % chunk_size == 0 {
            trip / chunk_size
        } else {
            trip / chunk_size + 1
        };

        if step > 0 {
            self.lb = start_point;
            self.ub = if self.lb + (chunk_size - 1) * step <= self.last {
                self.lb + (chunk_size - 1) * step
            } else {
                self.last
            };
        } else {
            self.ub = start_point;
            self.lb = if self.ub + (chunk_size - 1) * step >= self.trip {
                self.ub + (chunk_size - 1) * step
            } else {
                self.last
            };
        }

        for t in 0..THREAD_NUM {
            self.per_thread_start_point[t] = start_point + (chunk_size * step) * (t as i32);
        }
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
        if self.step > 0 {
            self.per_thread_start_point[tid as usize] <= self.last
        } else {
            self.per_thread_start_point[tid as usize] >= self.last
        }
    }

    pub fn get_next_static_chunk(&mut self, tid: u32) -> Chunk {
        // for static
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

    pub fn print(&self) {
        println!("lb: {}, ub: {}, chunk_size: {}, trip: {}, start: {}, last: {}, step: {}, avail_chunk: {}, per_thread_start_point: {:?}",
            self.lb, self.ub, self.chunk_size, self.trip, self.start, self.last, self.step, self.avail_chunk, self.per_thread_start_point);
    }
    
}