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

type Histogram = HashMap<i64, f64>; /// maybe can be replaced by a struct
type Chunk = (i32, i32); // maybe can be replaced by a HashMap

struct chunk_dispatcher{
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
    per_thread_start_point: [i32; 4],
}

impl chunk_dispatcher {
    
}