use std::arch::asm;
use std::collections::{HashMap, HashSet};
use std::ops::Deref;
use std::sync::{Arc, Mutex, Once};
use std::borrow::BorrowMut;
use std::cell::RefCell;
use std::collections::hash_map::Entry;
use lazy_static::lazy_static;
use statrs::distribution::{Discrete, NegativeBinomial};

type Histogram = HashMap<i64, f64>;

const THREAD_NUM: usize = 4;
const POLYBENCH_CACHE_SIZE_KB: usize = 2560;

lazy_static! {
    // pub static ref _NoSharePRI: [Mutex<Histogram>; THREAD_NUM] = Default::default();
    static ref _RIHist: Mutex<Histogram> = Mutex::new(Histogram::new());
    static ref _MRC: Mutex<HashMap<u64, f64>> = Mutex::new(HashMap::new());
    // static ref _NoSharePRI: Vec<Histogram> = vec![Histogram::new(); THREAD_NUM];
    static ref _SharePRI: [Mutex<HashMap<i64, Histogram>>; THREAD_NUM] = Default::default(); //FIXME: i changed the inner hashmap to i64, Histogram instead of i32, Histogram
}
/**thread local method start here */
thread_local! {static _NoSharePRI: RefCell<Histogram> = Default::default();} //FIXME: could change this to [RefCell<Histogram>; THREAD_NUM], but not sure how
// thread_local! {static _NoSharePRI: [RefCell<Histogram>; THREAD_NUM] = Default::default();}
/**thread local method end here */

/**sync once method start here */
// static mut _NoSharePRI: Option<Mutex<[Histogram; THREAD_NUM]>> = None;
// static INIT: Once = Once::new();
/**sync once method start here */

// static mut _NoSharePRI: [HashMap<i64, f64>; THREAD_NUM] = [HashMap::new(); THREAD_NUM];
// static mut some: [Histogram; THREAD_NUM] = [Histogram::new(); THREAD_NUM];

pub(crate) fn pluss_aet() {
    let mut p: HashMap<u64, f64> = HashMap::new();
    let mut histogram: HashMap<i64, f64> = HashMap::new();
    let mut total_num_rt = 0.0;
    let mut max_rt: i64 = 0;

    // let rihist_guard = RIHIST.lock().unwrap();
    let rihist_guard = _RIHist.lock().unwrap();
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
    // sorted_keys.reverse();
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
        if c > cs as u64 {
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
            let mut mrc_guard = _MRC.lock().unwrap();
            // let mut mrc_guard = &MRC;
            mrc_guard.insert(c, p[&prev_t]);
        } else if mrc_pred - p[&prev_t] < 0.0001 {
            // let mut mrc_guard = &MRC;
            let mut mrc_guard = _MRC.lock().unwrap();
            mrc_guard.insert(c, p[&prev_t]);
            mrc_pred = p[&prev_t];
        }
    }
}

pub(crate) fn pluss_cri_share_histogram_update(tid: i32, share_ratio: i32, reuse: i64, count: f64) {
    let share_ratio = share_ratio as i64;
    let mut local_share_pri = _SharePRI[tid as usize].lock().unwrap();
    if local_share_pri.contains_key(&share_ratio) {
        let mut histogram = local_share_pri.get(&share_ratio).unwrap().clone();
        _pluss_histogram_update(&mut histogram, reuse, count, Some(false));
        local_share_pri.insert(share_ratio, histogram);
    } else {
        let mut histogram = HashMap::new();
        histogram.insert(reuse, count);
        local_share_pri.insert(share_ratio, histogram);
    }
}

pub(crate) fn pluss_cri_noshare_histogram_update(tid: usize, reuse: i64, cnt: f64, in_log_format: Option<bool>) {
    // unsafe {
    //     let in_log_format = in_log_format.unwrap_or(true);
    //     let mut local_reuse = reuse;
    //     // let mut histogram = _NoSharePRI[tid as usize];
    //     if local_reuse > 0 && in_log_format {
    //         local_reuse = _polybench_to_highest_power_of_two(local_reuse);
    //     }
    //     if _NoSharePRI[tid as usize].contains_key(&local_reuse) {
    //         let update = _NoSharePRI[tid as usize].get(&local_reuse).unwrap() + cnt;
    //         _NoSharePRI[tid as usize].insert(local_reuse, update);
    //         // *histogram[tid].get_mut(&reuse).unwrap() += cnt; // this does the same thing as above
    //     } else {
    //         _NoSharePRI[tid as usize].insert(local_reuse, cnt);
    //     }
    //     // println!("histogram: {:?}", histogram);
    // }
    let in_log_format = in_log_format.unwrap_or(true);
    let mut local_reuse = reuse;
    // let mut histogram = _NoSharePRI[tid as usize];
    if local_reuse > 0 && in_log_format {
        local_reuse = _polybench_to_highest_power_of_two(local_reuse);
    }
    _NoSharePRI.with(|histogram| {
        match histogram.borrow_mut().entry(local_reuse) {
            Entry::Occupied(mut entry) => {
                *entry.get_mut() += cnt;
            },
            Entry::Vacant(entry) => {
                entry.insert(cnt);
            }
        }
        // if histogram.borrow().get(tid).unwrap().contains_key(&local_reuse) {
        //     let update = histogram.borrow().get(tid).unwrap().get(&local_reuse).unwrap() + cnt;
        //     *histogram.borrow_mut().get_mut(tid).unwrap().get_mut(&local_reuse).unwrap() = update;
        //     // *histogram[tid].get_mut(&reuse).unwrap() += cnt; // this does the same thing as above
        // } else {
        //     histogram.borrow_mut().get_mut(tid).unwrap()
        //         .insert(local_reuse, cnt);
        // }
    });
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

pub(crate) fn _pluss_histogram_update(histogram: &mut Histogram, reuse: i64, cnt: f64, in_log_format: Option<bool>) {
    let in_log_format = in_log_format.unwrap_or(true);
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
    // for i in 0..THREAD_NUM {
    //     // let histogram = _NoSharePRI[i].lock().unwrap();
    //     _NoSharePRI.with(|histogram| {
    //         for (k, v) in histogram.borrow().iter() {
    //             _pluss_histogram_update(&mut noshare_rih_tmp, *k, *v, Some(false));
    //         }
    //     });
    //     // for (k, v) in _NoSharePRI[i].iter() {
    //     //     _pluss_histogram_update(&mut noshare_rih_tmp, *k, *v, Some(false));
    //     // }
    // }
    _pluss_histogram_print("Start to dump noshare private reuse time", &noshare_rih_tmp);
}

pub(crate) fn pluss_cri_share_print_histogram() {
    let mut share_rih_tmp = HashMap::new();
    for i in 0..THREAD_NUM {
        let local_share_pri = _SharePRI[i].lock().unwrap();
        for (_, v) in local_share_pri.iter() {
            for (kk, vv) in v {
                _pluss_histogram_update(&mut share_rih_tmp, *kk, *vv, Some(false));
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
    println!("miss ratio");
    let mrc = _MRC.lock().unwrap();
    let mut it1 = mrc.iter();
    let mut it2 = mrc.iter();
    while !it1.clone().next().is_none() {
        let (mut it1_first, mut it1_second) = it1.next().unwrap();
        let (mut it2_first, mut it2_second) = it2.next().unwrap();
        loop {
            let mut it3 = it2.clone();
            let (mut it3_first, mut it3_second) = it3.next().unwrap_or((&0, &0.));
            it3.next();
            if it3.clone().next().is_none() {
                break;
            }
            if (it1_second - it3_second) < 0.00001 {
                it2.next();
            } else {
                break;
            }
        }
        println!("{}, {}", it1_first, it1_second);
        if !it1.eq(it2.clone()) {
            println!("{}, {}", it2_first, it2_second);
        }
        it2.next();
        it1 = it2.clone();
    }
}

pub(crate) fn _pluss_cri_nbd(thread_cnt: i32, n: i64, dist: &mut Histogram) {
    let p = 1.0 / thread_cnt as f64;
    let thread_cnt = thread_cnt as i64;
    let compare = ((4000 as f64) * (thread_cnt - 1) as f64) / thread_cnt as f64;
    if n as f64 >= compare { // FIXME: very strange
        let tmp = n * thread_cnt;
        // println!("n: {}, tmp: {}, thread_cnt: {}", n, tmp, thread_cnt);
        dist.insert(tmp, 1.0);
        return;
    }
    let mut k: i64 = 0;
    let mut nbd_prob: f64 = 0.0;
    let mut prob_sum: f64 = 0.0;
    let nb = NegativeBinomial::new(n as f64, p).unwrap();
    loop {
        nbd_prob = nb.pmf(k as u64); // this should work
        prob_sum += nbd_prob;
        dist.insert(k + n, nbd_prob);
        if prob_sum > 0.9999 {
            break;
        }
        k += 1;
    }
}

pub(crate) fn _pluss_cri_racetrack(thread_cnt: Option<i32>) {
    let mut thread_cnt = thread_cnt.unwrap_or(THREAD_NUM as i32);
    let mut merged_dist: HashMap<i64, Histogram> = Default::default(); //FIXME: i also changed this i32 to i64 as well
    for i in 0..thread_cnt as usize {
        let mut hash_of_hash = _SharePRI[i].lock().unwrap();
        for (k, hash) in hash_of_hash.iter() { // share entry
            if merged_dist.contains_key(k) {
                for (k1, v1) in hash { // reuse entry
                    *merged_dist.get_mut(k).unwrap().get_mut(k1).unwrap() += v1;
                }
            } else {
                for (k1, v1) in hash { // reuse entry
                    let mut tmp_hist = Histogram::new();
                    tmp_hist.insert(*k1, *v1);
                    merged_dist.insert(*k, tmp_hist);
                    // merged_dist.entry(*k).or_insert_with(|| Histogram::new()).insert(*k1, *v1); this is the same as above
                    // merged_dist.get_mut(k).unwrap().insert(*k1, *v1);
                }
            }
        }
    }
    let mut prob: HashMap<i32, f64> = Default::default();
    let mut dist = Histogram::new();
    let mut i: i32 = 1;
    let mut prob_sum: f64 = 0.0;
    let mut n: f64;
    let mut ri_to_distribute: i64;
    let mut cnt_to_distribute: f64;
    for (k, v) in merged_dist { // share entry
        n = k as f64;
        for (k1, v1) in v { // reuse entry
            if thread_cnt > 1 {
                // println!("dist before: {:?}, k1: {}", dist, k1);
                _pluss_cri_nbd(thread_cnt, k1, &mut dist);
                // println!("dist after: {:?}", dist);
                for (k2, v2) in dist.clone() { // dist entry
                    ri_to_distribute = k2;
                    cnt_to_distribute = v1 * v2;
                    prob_sum = 0.0;
                    i = 1;
                    loop {
                        if (2.0 as f64).powf(i as f64) > ri_to_distribute as f64 { break; }
                        prob.insert(i, (1 as f64 - (2.0 as f64).powf((i - 1) as f64) / ri_to_distribute as f64).powf(n) - (1 as f64 - (2.0 as f64).powf(i as f64) / ri_to_distribute as f64).powf(n));
                        prob_sum += prob.get(&i).unwrap();
                        i += 1;
                        if prob_sum == 1.0 { break; }
                    } // end of loop
                    if prob_sum != 1.0 {
                        prob.insert(i - 1, 1.0 - prob_sum);
                    }
                    for (k3, v3) in prob.clone() { // bin
                        let new_ri = (2.0 as f64).powf(k3 as f64 - 1.0) as i64;
                        // println!("new_ri:{}, v3:{}, cnt_to_distribute:{}, v3 * cnt_to_dist: {}", new_ri, v3, cnt_to_distribute, v3 * cnt_to_distribute);
                        pluss_histogram_update(new_ri, v3 * cnt_to_distribute);
                    }
                    prob.clear();
                } // end of iterating dist
                dist.clear();
            } else {
                pluss_histogram_update(k1, v1);
            }
        } // end of iterating all reuse entries
    } // end of iterating all type of reuses
} // end of void pluss_cri_racetrack()

pub(crate) fn pluss_histogram_update(reuse: i64, cnt: f64) {
    _pluss_histogram_update(&mut _RIHist.lock().unwrap(), reuse, cnt, None);
}

pub(crate) fn _pluss_cri_noshare_distribute(thread_cnt: Option<i32>) {
    let mut thread_cnt = thread_cnt.unwrap_or(THREAD_NUM as i32);
    let mut dist = Histogram::new();
    let mut merge_dist = Histogram::new();
    for i in 0..thread_cnt as usize {
        // let mut hist = _NoSharePRI[i].lock().unwrap();
        _NoSharePRI.with(|histogram| {
            for (k, v) in histogram.borrow().iter() {
                if dist.contains_key(k) {
                    *dist.get_mut(k).unwrap() += v;
                } else {
                    dist.insert(*k, *v);
                }
            }
        })
        // for (k, v) in _NoSharePRI[i].iter() {
        //     if merge_dist.contains_key(k) {
        //         // *merge_dist.get_mut(k).unwrap() += v;
        //         let tmp = merge_dist.get(&k).unwrap() + v;// this is the same as the method above
        //         merge_dist.insert(*k, tmp);
        //     } else {
        //         merge_dist.insert(*k, *v);
        //     }
        // }
    }
    let mut dist_clone;
    for (k, v) in merge_dist.clone() {
        if k < 0 {
            pluss_histogram_update(k, v);
            continue;
        }
        if thread_cnt > 1 {
            dist_clone = dist.clone();
            _pluss_cri_nbd(thread_cnt, k, &mut dist_clone); // update dist
            dist = dist_clone.clone();
            for (kk, vv) in dist.clone() {
                let ri_to_distribute = kk;
                pluss_histogram_update(ri_to_distribute, v * vv);
            }
            dist.clear(); //FIXME: i am not sure if this is correct
        } else {
            pluss_histogram_update(k, v);
        } // end of if(thread_cnt > 1)
        // println!("pluss_cri_noshare_distribute: {:?}", merge_dist.clone());
    }
    // println!("pluss_cri_noshare_distribute: {:?}", merge_dist.clone());
} // end of void pluss_cri_noshare_distribute()

pub(crate) fn pluss_cri_distribute(thread_cnt: i32) {
    // _pluss_cri_noshare_distribute(Some(thread_cnt));
    _pluss_cri_racetrack(Some(thread_cnt));
}