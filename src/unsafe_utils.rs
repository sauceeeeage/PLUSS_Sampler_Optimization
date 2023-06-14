use std::arch::asm;
use std::ops::Deref;
use std::sync::{Arc, Mutex, Once};
use std::borrow::BorrowMut;
use std::cell::RefCell;
use std::collections::hash_map::{Entry};
use ahash::RandomState;
use std::io::Read;
use std::thread;
use lazy_static::lazy_static;
use rayon::prelude::*;
use statrs::distribution::{Discrete, NegativeBinomial};
use ahash::HashMapExt;

type HashSet<K> = std::collections::HashSet<K, RandomState>;
type HashMap<K, V> = std::collections::HashMap<K, V, RandomState>;
type Histogram = HashMap<i64, f64>;

const THREAD_NUM: usize = 4;
const POLYBENCH_CACHE_SIZE_KB: usize = 2560;

lazy_static! {
    // pub static ref _NoSharePRI: [Mutex<Histogram>; THREAD_NUM] = Default::default();
    static ref _RIHist: Mutex<Histogram> = Default::default();
    // static ref _MRC: Mutex<HashMap<u64, f64>> = Mutex::new(HashMap::new());
    // static ref _SharePRI: [Mutex<HashMap<i64, Histogram>>; THREAD_NUM] = Default::default();
    static ref _MergeSharePRI: Mutex<HashMap<i64, Histogram>> = Default::default();
    static ref _MergeNoSharePRI: Mutex<Histogram> = Default::default();
}

/**thread local method start here */
thread_local! {
    static _NoSharePRI: RefCell<Histogram> = Default::default();
    static _SharePRI: RefCell<HashMap<i64, Histogram>> = Default::default(); // share ratio, histogram
}
/**thread local method end here */

// pub(crate) fn pluss_aet() {
//     let mut p: HashMap<u64, f64> = HashMap::new();
//     let mut histogram: HashMap<i64, f64> = HashMap::new();
//     let mut total_num_rt = 0.0;
//     let mut max_rt: i64 = 0;
//
//     // let rihist_guard = RIHIST.lock().unwrap();
//     let rihist_guard = _RIHist.lock().unwrap();
//     for (k, v) in rihist_guard.iter() {
//         total_num_rt += *v;
//         histogram.insert(*k, *v);
//         if max_rt < *k {
//             max_rt = *k;
//         }
//     }
//
//     let mut accumulate_num_rt = 0.0;
//     if histogram.contains_key(&(-1)) {
//         accumulate_num_rt = *histogram.get(&(-1)).unwrap();
//     }
//
//     let mut sorted_keys: Vec<&i64> = histogram.keys().collect();
//     sorted_keys.sort();
//     // sorted_keys.reverse();
//     for k in sorted_keys {
//         let v = histogram.get(k).unwrap();
//         if *k == -1 {
//             break;
//         }
//         p.insert(*k as u64, accumulate_num_rt / total_num_rt);
//         accumulate_num_rt += v;
//     }
//
//     p.insert(0, 1.0);
//     let mut sum_p = 0.0;
//     let mut mrc_pred = -1.0;
//     let mut t: u64 = 0;
//     let mut prev_t: u64 = 0;
//     let cs = 2560 * 1024 / std::mem::size_of::<f64>();
//
//     for c in 0..=max_rt as u64 {
//         if c > cs as u64 {
//             break;
//         }
//         while sum_p < c as f64 && t <= max_rt as u64 {
//             if p.contains_key(&t) {
//                 sum_p += p[&t];
//                 prev_t = t;
//             } else {
//                 sum_p += p[&prev_t];
//             }
//             t += 1;
//         }
//
//         if mrc_pred != -1.0 {
//             let mut mrc_guard = _MRC.lock().unwrap();
//             // let mut mrc_guard = &MRC;
//             mrc_guard.insert(c, p[&prev_t]);
//         } else if mrc_pred - p[&prev_t] < 0.0001 {
//             // let mut mrc_guard = &MRC;
//             let mut mrc_guard = _MRC.lock().unwrap();
//             mrc_guard.insert(c, p[&prev_t]);
//             mrc_pred = p[&prev_t];
//         }
//     }
// }
#[inline]
pub(crate) fn sync_pri() {
    sync_noshare_pri();
    sync_share_pri();
}

fn sync_share_pri() {
    let mut local_merge_share_pri = _MergeSharePRI.lock().unwrap();
    _SharePRI.with(|local_share_pri| {
        if local_merge_share_pri.is_empty() {
            local_merge_share_pri.extend(local_share_pri.borrow().deref().clone());
        } else {
            local_share_pri.borrow().iter().for_each(|(k, v)| {
                match local_merge_share_pri.entry(*k) {
                    Entry::Occupied(mut e) => {
                        for (i, j) in v.iter() {
                            e.get_mut().entry(*i).and_modify(|x| *x += j).or_insert(*j);
                        }
                    }
                    Entry::Vacant(e) => {
                        e.insert(v.clone());
                    }
                }
            });
        }
    });
}

fn sync_noshare_pri() {
    let mut local_merge_noshare_pri = _MergeNoSharePRI.lock().unwrap();
    // println!("before adding local_merge_noshare_pri: {:?}", local_merge_noshare_pri);
    _NoSharePRI.with(|noshare_pri| {
        // println!("noshare_pri: {:?}", noshare_pri.borrow().deref().clone());
        if local_merge_noshare_pri.is_empty() {
            // println!("local_merge_noshare_pri is empty");
            local_merge_noshare_pri.extend(noshare_pri.borrow().deref().clone()); //FIXME: not sure if we need deref
        } else {
            // println!("local_merge_noshare_pri is not empty");
            // local_merge_noshare_pri.borrow_mut().into_iter().chain(&noshare_pri.borrow().iter()).for_each(|(k, v)| {
            //     *v += *noshare_pri.borrow().get(*k).unwrap();
            // }); FIXME: can we do it this way, is this more efficient?
            for (k, v) in noshare_pri.borrow().iter() {
                local_merge_noshare_pri.entry(*k).and_modify(|e| *e += *v).or_insert(*v);
            }
        }
    });
    // println!("after adding local_merge_noshare_pri: {:?}", local_merge_noshare_pri);
}

#[inline]
pub(crate) fn pluss_cri_share_histogram_update(tid: i32, share_ratio: i32, reuse: i64, count: f64) {
    let share_ratio = share_ratio as i64;
    // let mut local_share_pri = _SharePRI[tid as usize].lock().unwrap();

    // if local_share_pri.contains_key(&share_ratio) {
    //     let mut histogram = local_share_pri.get(&share_ratio).unwrap().clone();
    //     _pluss_histogram_update(&mut histogram, reuse, count, Some(false));
    //     local_share_pri.insert(share_ratio, histogram);
    // } else {
    //     let mut histogram = HashMap::new();
    //     histogram.insert(reuse, count);
    //     local_share_pri.insert(share_ratio, histogram);
    // }

    // match local_share_pri.entry(share_ratio) {
    //     Entry::Occupied(mut e) => {
    //         let mut histogram = e.get_mut();
    //         _pluss_histogram_update(&mut histogram, reuse, count, Some(false));
    //         //should take care of the insert, since this is get mut
    //     }
    //     Entry::Vacant(e) => {
    //         let mut histogram = HashMap::new();
    //         histogram.insert(reuse, count);
    //         e.insert(histogram);
    //     }
    // }

    _SharePRI.with(|share_pri| {
        match share_pri.borrow_mut().entry(share_ratio) {
            Entry::Occupied(mut e) => {
                let mut histogram = e.get_mut();
                _pluss_histogram_update(&mut histogram, reuse, count, Some(false));
                //should take care of the insert, since this is get mut
            }
            Entry::Vacant(e) => {
                let mut histogram = HashMap::new();
                histogram.insert(reuse, count);
                e.insert(histogram);
            }
        }
    });
}

#[inline]
pub(crate) fn pluss_cri_noshare_histogram_update(tid: usize, reuse: i64, cnt: f64, in_log_format: Option<bool>) {
    let in_log_format = in_log_format.unwrap_or(true);
    let mut local_reuse = reuse;
    if local_reuse > 0 && in_log_format {
        local_reuse = _polybench_to_highest_power_of_two(local_reuse);
    }
    _NoSharePRI.with(|histogram| {
        // if local_reuse == -1 {
        //     println!("before tid: {}, cnt: {}, hist: {:?}", tid, cnt, histogram);
        // }
        // println!("tid: {}, hist: {:?}", tid, histogram);
        // if histogram.borrow().get(&256) == Some(&65024.) {
        //     println!("tid: {}, hist: {:?}", tid, histogram);
        // }
        match histogram.borrow_mut().entry(local_reuse) {
            Entry::Occupied(mut entry) => {
                *entry.get_mut() += cnt;
            }
            Entry::Vacant(entry) => {
                entry.insert(cnt);
            }
        }
        // if local_reuse == -1 {
        //     println!("after tid: {}, cnt: {}, hist: {:?}", tid, cnt, histogram);
        // }
        // println!("tid: {}, hist: {:?}", tid, histogram);
    });
}

#[inline(always)]
pub(crate) fn _polybench_to_highest_power_of_two(mut x: i64) -> i64 {
    (x as u64).next_power_of_two() as i64
}

#[inline]
pub(crate) fn _pluss_histogram_print(title: &str, histogram: &Histogram) {
    println!("{}", title);
    let sum: f64 = histogram.values().sum(); // sum of iter->second
    for (k, v) in histogram {
        println!("{},{},{}", k, v, v / sum);
    }
}

#[inline]
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

#[inline]
pub(crate) fn pluss_cri_noshare_print_histogram() {
    // let mut noshare_rih_tmp = HashMap::new();
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

    // let mut handles = vec![];
    // for i in 0..THREAD_NUM {
    //     let handle = thread::spawn(move || {
    //         let histogram = _MergeNoSharePRI.lock().unwrap();
    //         _NoSharePRI.with(|local_histogram| {
    //             // match local_histogram.borrow_mut().entry() {
    //             // }
    //         });
    //     });
    //     handles.push(handle);
    // }


    _pluss_histogram_print("Start to dump noshare private reuse time", &_MergeNoSharePRI.lock().unwrap());
}

#[inline]
pub(crate) fn pluss_cri_share_print_histogram() {
    let mut share_rih_tmp = HashMap::new();
    // for i in 0..THREAD_NUM {
    //     let local_share_pri = _SharePRI[i].lock().unwrap();
    //     for (_, v) in local_share_pri.iter() {
    //         for (kk, vv) in v {
    //             _pluss_histogram_update(&mut share_rih_tmp, *kk, *vv, Some(false));
    //         }
    //     }
    // }
    let local_share_pri = _MergeSharePRI.lock().unwrap();
    local_share_pri.iter().for_each(|(_, v)| { // {(_ {_, _}), ....}
        share_rih_tmp.extend(v);
    });
    _pluss_histogram_print("Start to dump share private reuse time", &share_rih_tmp);
}

#[inline]
pub(crate) fn pluss_print_histogram() {
    let rihist = _RIHist.lock().unwrap();
    _pluss_histogram_print("Start to dump reuse time", &rihist);
}

// pub(crate) fn pluss_print_mrc() {
//     println!("miss ratio");
//     let mrc = _MRC.lock().unwrap();
//     let mut it1 = mrc.iter();
//     let mut it2 = mrc.iter();
//     while !it1.clone().next().is_none() {
//         let (mut it1_first, mut it1_second) = it1.next().unwrap();
//         let (mut it2_first, mut it2_second) = it2.next().unwrap();
//         loop {
//             let mut it3 = it2.clone();
//             let (mut it3_first, mut it3_second) = it3.next().unwrap_or((&0, &0.));
//             it3.next();
//             if it3.clone().next().is_none() {
//                 break;
//             }
//             if (it1_second - it3_second) < 0.00001 {
//                 it2.next();
//             } else {
//                 break;
//             }
//         }
//         println!("{}, {}", it1_first, it1_second);
//         if !it1.eq(it2.clone()) {
//             println!("{}, {}", it2_first, it2_second);
//         }
//         it2.next();
//         it1 = it2.clone();
//     }
// }

#[inline]
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
    let mut merged_dist = _MergeSharePRI.lock().unwrap();
    // for i in 0..thread_cnt as usize {
    //     let mut hash_of_hash = _NoSharePRI[i].lock().unwrap();
    //     println!("hash_of_hash: {:?}", hash_of_hash);
    //     for (k, hash) in hash_of_hash.iter() { // share entry
    //         if merged_dist.contains_key(k) {
    //             for (k1, v1) in hash { // reuse entry
    //                 *merged_dist.get_mut(k).unwrap().get_mut(k1).unwrap() += v1;
    //             }
    //         } else {
    //             for (k1, v1) in hash { // reuse entry
    //                 let mut tmp_hist = Histogram::new();
    //                 tmp_hist.insert(*k1, *v1);
    //                 merged_dist.insert(*k, tmp_hist);
    //                 // merged_dist.entry(*k).or_insert_with(|| Histogram::new()).insert(*k1, *v1); this is the same as above
    //                 // merged_dist.get_mut(k).unwrap().insert(*k1, *v1);
    //             }
    //         }
    //     }
    // }
    // println!("merged_dist: {:?}", merged_dist);
    let mut prob: HashMap<i32, f64> = Default::default();
    let mut dist = Histogram::new();
    let mut i: i32 = 1;
    let mut prob_sum: f64 = 0.0;
    let mut n: f64;
    let mut ri_to_distribute: i64;
    let mut cnt_to_distribute: f64;
    for (k, v) in merged_dist.iter() { // share entry
        n = *k as f64;
        for (k1, v1) in v { // reuse entry
            if thread_cnt > 1 {
                // println!("dist before: {:?}, k1: {}", dist, k1);
                _pluss_cri_nbd(thread_cnt, *k1, &mut dist);
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
                        // println!("ri_hist: {:?}", _RIHist.lock().unwrap());
                    }
                    prob.clear();
                } // end of iterating dist
                dist.clear();
            } else {
                pluss_histogram_update(*k1, *v1);
                // println!("ri_hist: {:?}", _RIHist.lock().unwrap());
            }
        } // end of iterating all reuse entries
    } // end of iterating all type of reuses
} // end of void pluss_cri_racetrack()

#[inline]
pub(crate) fn pluss_histogram_update(reuse: i64, cnt: f64) {
    _pluss_histogram_update(&mut _RIHist.lock().unwrap(), reuse, cnt, None);
}

#[inline]
pub(crate) fn _pluss_cri_noshare_distribute(thread_cnt: Option<i32>) {
    let mut thread_cnt = thread_cnt.unwrap_or(THREAD_NUM as i32);
    let mut dist = Histogram::new();
    let mut merge_dist = _MergeNoSharePRI.lock().unwrap();
    // for i in 0..thread_cnt as usize { //FIXME: this for loop is wrong b/c the thread local doesn't have thread num of thread objects
    //     // let mut hist = _NoSharePRI[i].lock().unwrap();
    //     _NoSharePRI.with(|histogram| {
    //         for (k, v) in histogram.borrow().iter() {
    //             if dist.contains_key(k) {
    //                 *dist.get_mut(k).unwrap() += v;
    //             } else {
    //                 dist.insert(*k, *v);
    //             }
    //         }
    //     })
    //     // for (k, v) in _NoSharePRI[i].iter() {
    //     //     if merge_dist.contains_key(k) {
    //     //         // *merge_dist.get_mut(k).unwrap() += v;
    //     //         let tmp = merge_dist.get(&k).unwrap() + v;// this is the same as the method above
    //     //         merge_dist.insert(*k, tmp);
    //     //     } else {
    //     //         merge_dist.insert(*k, *v);
    //     //     }
    //     // }
    // }
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

#[inline]
pub(crate) fn pluss_cri_distribute(thread_cnt: i32) {
    // println!("before distribute");
    // println!("ri_hist: {:?}", _RIHist.lock().unwrap());
    // println!("NoSharePRI: {:?}", _MergeNoSharePRI.lock().unwrap());
    // println!("SharePRI: {:?}", _MergeSharePRI.lock().unwrap());
    _pluss_cri_noshare_distribute(Some(thread_cnt));
    // println!("after distribute");
    // println!("before race track");
    // println!("ri_hist: {:?}", _RIHist.lock().unwrap());
    // println!("NoSharePRI: {:?}", _MergeNoSharePRI.lock().unwrap());
    // println!("SharePRI: {:?}", _MergeSharePRI.lock().unwrap());
    _pluss_cri_racetrack(Some(thread_cnt));
    // println!("after race track");
    // println!("ri_hist: {:?}", _RIHist.lock().unwrap());
    // println!("NoSharePRI: {:?}", _MergeNoSharePRI.lock().unwrap());
    // println!("SharePRI: {:?}", _MergeSharePRI.lock().unwrap());
}