use std::cmp::{Ordering, min};
use std::hash::{Hash, Hasher};

const CHUNK_SIZE: usize = 4;
const THREAD_NUM: usize = 4;

#[derive(Clone, Debug, Eq)]
pub struct Iteration {
    pub name: String,
    pub ivs: Vec<usize>, // iteration variables
    pub pidx: usize,
    pub priority: i32,
    pub cid: usize,
    pub tid: usize,
    pub pos: usize,
    pub is_parallel_iteration: bool,
}

impl Iteration {
    pub fn new(
        name: &str,
        ivs: Vec<usize>,
        start: usize,
        step: usize,
        is_in_parallel_region: bool,
        parallel_iter_idx: usize,
        priority: Option<i32>,
    ) -> Self {
        // By default, the priority is 1, start is 0, step is 1, parallel_iter_idx is 0, and is_in_parallel_region is false
        let priority = priority.unwrap_or(1);
        let (cid, tid, pos) = if is_in_parallel_region {
            // FIXME: this might be wrong, b/c all variables are usize, this might overflow; maybe change it to isize?
            let cid = ((ivs[parallel_iter_idx] - start) / step) / (CHUNK_SIZE * THREAD_NUM);
            let tid = ((ivs[parallel_iter_idx] - start) / step) / CHUNK_SIZE - cid * THREAD_NUM;
            let pos = ((ivs[parallel_iter_idx] - start) / step) % CHUNK_SIZE;
            (cid, tid, pos)
        } else {
            (0, 0, 0) /* dummy values */
        };

        Self {
            name: name.to_string(),
            ivs,
            pidx: parallel_iter_idx,
            priority,
            cid,
            tid,
            pos,
            is_parallel_iteration: is_in_parallel_region,
        }
    }

    pub fn to_string(&self) -> String {
        let ivs_str: Vec<String> = self.ivs.iter().map(|iv| iv.to_string()).collect();
        format!("{} ({})", self.name, ivs_str.join(","))
    }

    pub fn to_addr_string(&self) -> String {
        let ivs_str: Vec<String> = self.ivs.iter().map(|iv| iv.to_string()).collect();
        format!("{}{}", self.name, ivs_str.join("-"))
    }

    pub fn compare(&self, other: &Iteration) -> Ordering {
        if self.ivs.len() != other.ivs.len() {
            return Ordering::Equal;
        }

        if self.is_parallel_iteration {
            match self.cid.cmp(&other.cid) {
                Ordering::Less => return Ordering::Less,
                Ordering::Greater => return Ordering::Greater,
                _ => {}
            }

            match self.pos.cmp(&other.pos) {
                Ordering::Less => return Ordering::Less,
                Ordering::Greater => return Ordering::Greater,
                _ => {}
            }
        }

        let common_ivs_size = min(self.ivs.len(), other.ivs.len());
        let ivs_idx = if self.is_parallel_iteration {
            let mut ivs_idx = 0;

            while ivs_idx < common_ivs_size && ivs_idx < self.pidx {
                match self.ivs[ivs_idx].cmp(&other.ivs[ivs_idx]) {
                    Ordering::Less => return Ordering::Less,
                    Ordering::Greater => return Ordering::Greater,
                    _ => {}
                }
                ivs_idx += 1;
            }

            if ivs_idx == self.pidx {
                ivs_idx += 1; // skip the parallel iteration
            }

            // compare the iteration inside the parallel region
            while ivs_idx < common_ivs_size {
                match self.ivs[ivs_idx].cmp(&other.ivs[ivs_idx]) {
                    Ordering::Less => return Ordering::Less,
                    Ordering::Greater => return Ordering::Greater,
                    _ => {}
                }
                ivs_idx += 1;
            }

            ivs_idx // return the current index
        } else {
            0 // in a non-parallel region, start from the beginning
        };

        // compare the remaining loop induction variables
        for idx in ivs_idx..common_ivs_size {
            match self.ivs[idx].cmp(&other.ivs[idx]) {
                Ordering::Less => return Ordering::Less,
                Ordering::Greater => return Ordering::Greater,
                _ => {}
            }
        }

        // priority indicates the topology order of the references, higher means
        // this reference is accessed early
        match self.priority.cmp(&other.priority) {
            Ordering::Less => return Ordering::Greater,
            Ordering::Greater => return Ordering::Less,
            _ => {}
        }

        // all equal, these two samples are equal
        Ordering::Equal
    }

}

impl PartialEq for Iteration {
    fn eq(&self, other: &Self) -> bool {
        if self.ivs.len() != other.ivs.len() {
            return false;
        }
        for i in 0..self.ivs.len() {
            if self.ivs[i] != other.ivs[i] {
                return false;
            }
        }
        self.name == other.name
    }
}
// IterationComp
impl PartialOrd for Iteration {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let mut ordering = Ordering::Equal;

        // compare the chunk id
        if self.is_parallel_iteration {
            ordering = self.cid.cmp(&other.cid);

            // compare pos and tid only works for uniform interleaving + static
            // the same chunk. compare the thread local position.
            if ordering == Ordering::Equal {
                ordering = self.pos.cmp(&other.pos);
            }
        }

        let common_ivs_size = min(self.ivs.len(), other.ivs.len());
        let mut ivs_idx = 0;

        while ordering == Ordering::Equal && ivs_idx < common_ivs_size {
            if self.is_parallel_iteration && ivs_idx == self.pidx {
                ivs_idx += 1;
                continue;
            }

            ordering = self.ivs[ivs_idx].cmp(&other.ivs[ivs_idx]);
            ivs_idx += 1;
        }

        if ordering == Ordering::Equal {
            if self.is_parallel_iteration {
                // the same thread local position. compare the thread
                ordering = self.tid.cmp(&other.tid);
            }

            if ordering == Ordering::Equal {
                // priority indicates the topology order of the references, higher means
                // this reference is accessed early
                ordering = self.priority.cmp(&other.priority).reverse();
            }
        }

        Some(ordering)
    }
}


// IterationHasher
impl Hash for Iteration {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
        let mut bitmap: u64 = 0u64;
        let mut i = 2;
        for &iv in &self.ivs {
            bitmap |= (iv as u64) << (i * 14);
            i -= 1;
            if i < 0 {
                break;
            }
        }
        bitmap.hash(state); // not sure to be honest......
    }
}
