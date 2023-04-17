use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use lazy_static::lazy_static;

type Histogram = HashMap<i64, f64>;
const THREAD_NUM: usize = 4;

lazy_static! {
    pub static ref _NoSharePRI: Arc<Mutex<[Histogram; THREAD_NUM]>> = Default::default();
    static ref _RIHist: Mutex<Histogram> = Mutex::new(Histogram::new());
    static ref _MRC: Mutex<HashMap<u64, f64>> = Mutex::new(HashMap::new());
    // static ref _NoSharePRI: Vec<Histogram> = vec![Histogram::new(); THREAD_NUM];
    static ref _SharePRI: Vec<HashMap<i32, Histogram>> = vec![HashMap::new(); THREAD_NUM];
}

pub(crate) fn pluss_cri_noshare_histogram_update(tid: usize, reuse: i64, cnt: f64, in_log_format: Option<bool>) {
    let in_log_format = in_log_format.unwrap_or(true);
    let mut local_reuse = reuse; //FIXME: this may be wrong, it's in pluss_utils.h line 680
    let mut histogram = _NoSharePRI.lock().unwrap();
    if local_reuse > 0 && in_log_format {
        local_reuse = _polybench_to_highest_power_of_two(local_reuse);
    }
    if histogram[tid].contains_key(&local_reuse) {
        let update = histogram[tid].get(&local_reuse).unwrap() + cnt;
        histogram[tid].insert(local_reuse, update);
        // *histogram[tid].get_mut(&reuse).unwrap() += cnt; // this does the same thing as above
    } else {
        histogram[tid].insert(local_reuse, cnt);
    }
}

pub(crate) fn _polybench_to_highest_power_of_two(mut x: i64) -> i64 {
    // Check for the set bits
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    // Then we remove all but the top bit by xor'ing the
    // string of 1's with that string of 1's
    // shifted one to the left, and we end up with
    // just the one top bit followed by 0's
    x ^ (x >> 1)
}

pub(crate) fn _pluss_histogram_print(title: &str, histogram: &Histogram) {
    println!("{}", title);
    let sum: f64 = histogram.values().sum(); // sum of iter->second
    for (k, v) in histogram {
        println!("{},{},{}", k, v, v / sum);
    }
}

pub(crate) fn _pluss_histogram_update(histogram: &mut Histogram, reuse: i64, cnt: f64, in_log_format: bool) {
    let reuse: i64 = if reuse > 0 && in_log_format {
        _polybench_to_highest_power_of_two(reuse)
    } else {
        reuse
    };

    let entry = histogram.entry(reuse).or_insert(0.);
    *entry += cnt;
}

pub(crate) fn pluss_cri_noshare_print_histogram() {
    let mut noshare_rih_tmp = HashMap::new();
    for i in 0..THREAD_NUM {
        let histogram = _NoSharePRI.lock().unwrap();
        for (k, v) in &histogram[i] {
            _pluss_histogram_update(&mut noshare_rih_tmp, *k, *v, false);
        }
    }
    _pluss_histogram_print("Start to dump noshare private reuse time", &noshare_rih_tmp);
}

pub(crate) fn _pluss_cri_share_print_histogram() {
    let mut share_rih_tmp = HashMap::new();
    for i in 0..THREAD_NUM {
        for (_, v) in &_SharePRI[i] {
            for (kk, vv) in v {
                _pluss_histogram_update(&mut share_rih_tmp, *kk, *vv, false);
            }
        }
    }
    _pluss_histogram_print("Start to dump share private reuse time", &share_rih_tmp);
}

pub(crate) fn pluss_print_histogram() {
    let rihist = _RIHist.lock().unwrap();
    _pluss_histogram_print("Start to dump reuse time", &rihist);
}

pub(crate) fn pluss_print_mrc() {
    //TODO: have not finished yet
    println!("miss ratio");
    let mrc = _MRC.lock().unwrap();
    let mut it1 = mrc.iter();

    while let Some((k1, v1)) = it1.next() {
        while let Some((_, v3)) = it1.clone().skip(1).next() { // not certain....
            println!("{}, {}", k1, v1);
            if v1 - v3 < 0.00001 {
                it1.next();
                println!("{}, {}", k1, v1);
            } else {
                break;
            }
        }
    }
}
