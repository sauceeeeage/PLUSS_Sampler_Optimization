use std::collections::HashMap;
use std::sync::{Arc, Mutex};

type Histogram = HashMap<i64, f64>;

// Arc<Mutex<T>> for globals
pub static mut RIHIST: Arc<Mutex<Histogram>> = Arc::new(Mutex::new(HashMap::new()));
pub static mut PAR_RI_HIST: Arc<Mutex<HashMap<i64, f64>>> = Arc::new(Mutex::new(HashMap::new()));
pub static mut MRC: Arc<Mutex<HashMap<u64, f64>>> = Arc::new(Mutex::new(HashMap::new()));

fn pluss_aet() {
    let mut p: HashMap<u64, f64> = HashMap::new();
    let mut histogram: HashMap<i64, f64> = HashMap::new();
    let mut total_num_rt = 0.0;
    let mut max_rt: i64 = 0;

    let rihist_guard = RIHIST.lock().unwrap();
    for (k, v) in rihist_guard.iter() {
        total_num_rt += *v;
        histogram.insert(*k, *v);
        if max_rt < *k {
            max_rt = *k;
        }
    }

    let mut accumulate_num_rt = 0.0;
    if histogram.contains_key(&(-1)) {
        accumulate_num_rt = *histogram.get(&(-1)).unwrap();
    }

    let mut sorted_keys: Vec<&i64> = histogram.keys().collect();
    sorted_keys.sort();
    sorted_keys.reverse();
    for k in sorted_keys {
        let v = histogram.get(k).unwrap();
        if *k == -1 {
            break;
        }
        p.insert(*k as u64, accumulate_num_rt / total_num_rt);
        accumulate_num_rt += v;
    }

    p.insert(0, 1.0);
    let mut sum_p = 0.0;
    let mut mrc_pred = -1.0;
    let mut t: u64 = 0;
    let mut prev_t: u64 = 0;
    let cs = 2560 * 1024 / std::mem::size_of::<f64>();

    for c in 0..=max_rt as u64 {
        if c > cs {
            break;
        }
        while sum_p < c as f64 && t <= max_rt as u64 {
            if p.contains_key(&t) {
                sum_p += p[&t];
                prev_t = t;
            } else {
                sum_p += p[&prev_t];
            }
            t += 1;
        }

        if mrc_pred != -1.0 {
            let mut mrc_guard = MRC.lock().unwrap();
            mrc_guard.insert(c, p[&prev_t]);
        } else if mrc_pred - p[&prev_t] < 0.0001 {
            let mut mrc_guard = MRC.lock().unwrap();
            mrc_guard.insert(c, p[&prev_t]);
            mrc_pred = p[&prev_t];
        }
    }
}
