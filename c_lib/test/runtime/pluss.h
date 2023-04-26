#ifndef PLUSS_H
#define PLUSS_H
#include <iostream>
#include <assert.h>

#ifndef CHUNK_SIZE
#define CHUNK_SIZE	4
#endif

#ifndef THREAD_NUM
#define THREAD_NUM	4
#endif	

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
#if 0
typedef pair<int, int> Chunk;

class Iteration {
public:	
	string name;
	vector<int> ivs;
	unsigned cid; // chunk this sample locates
	unsigned tid; // thread this sample belongs
	unsigned pos; // thread local position
	bool isParallelIteration;

	Iteration() {}
	~Iteration() {}
	Iteration(string ref, vector<int> iter, int start=0, int step=1, bool isInParallelRegion=false, int parallel_iter_idx=0) {
		name = ref;
		ivs = iter;
		// works when openmp static scheduling
        if (isInParallelRegion) {
            cid = floor(((iter[parallel_iter_idx]-start)/step) / (CHUNK_SIZE * THREAD_NUM));
            tid = ((iter[parallel_iter_idx]-start)/step) / CHUNK_SIZE - floor(((iter[parallel_iter_idx]-start)/step) / (CHUNK_SIZE * THREAD_NUM)) * THREAD_NUM;
            pos = ((iter[parallel_iter_idx]-start)/step) % CHUNK_SIZE;
        }
        isParallelIteration = isInParallelRegion;
	}
	string toString() {
		string s = name + " (";
		vector<int>::iterator it;
		for (it = ivs.begin(); it != ivs.end(); ++it) {
			s += to_string(*it) + ",";
		}
		s.pop_back();
		s += ")";
		return s;
	}
	int compare(Iteration other) {
		if (ivs.size() != other.ivs.size()) { return 2; }
		/* compare the chunk id */ ;
		if (isParallelIteration) {
            if (cid < other.cid) {
                return -1;
            } else if (cid > other.cid) {
                return 1;
            }
		// // compare of pos and tid only works for uniform interleaving + static 
		// /* the same chunk. compare the thread local position. */ ;
            if (pos < other.pos) {
                return -1;
            } else if (pos > other.pos) {
                return 1;
            }
		// /* the same thread local position. compare the thread */ ;
            if (tid < other.tid) {
                return -1;
            } else if (tid > other.tid) {
                return 1;
            }
		}
		/* the same thread. compare the rest loop induction variables */ ;
		vector<int>::iterator selfit = ivs.begin();
		vector<int>::iterator otherit = other.ivs.begin();
		while (selfit != ivs.end() && otherit != other.ivs.end()) {
			if (*selfit < *otherit) {
				return -1;
			} else if (*selfit > *otherit) {
				return 1;
			}
			selfit++;
			otherit++;
		}
		/* all equal, these two samples are equal */ ;
		return 0;
	}
	bool operator==(const Iteration & other) const {
		unsigned i = 0;
		for (i = 0; i < ivs.size(); i++) {
			if (ivs[i] != other.ivs[i]) {
				return false;
			}
		}
		return name == other.name;
	}
};
// Hash function for Iteration
struct IterationHasher {
	size_t operator()(const struct Iteration &iter) const {
		using std::string;
		using std::size_t;
		using std::hash;
		size_t hash_val = hash<string>()(iter.name);
		uint64_t bitmap = 0UL;
		int i = 2;
		for (auto iv : iter.ivs) {
			bitmap |= ((uint64_t)iv << (i * 14));
			i -= 1;
			if (i < 0) { break; }
		}
		hash_val ^= hash<uint64_t>()(bitmap);
		return hash_val;
	}
};

class ChunkDispatcher {
	int lb = 0;
	int ub = 0;
	int chunk_size = 0;
	int trip = 0;
	int start = 0;
	int last = 0;
	int step = 1;
	int avail_chunk = 0;
	vector<int> per_thread_start_point; // use for static sheduling only
public:
	ChunkDispatcher() {}

	// trip : total number of iterations traversed by the parallel loop
	// start_point: the iteration the parallel loop stops
	// 
	// the last iteration it will traverse will be start_point + k*step < trip
	ChunkDispatcher(int chunk_size, int trip, int start_point=0, int step=1) {
		// it is possible that the chunk_size is greater than the trip count
		this->chunk_size = chunk_size;
		this->trip = trip;
		this->start = start_point;
		this->step = step;
		this->last = start_point + (trip - 1) * step; // the last iteration
		// the range of the parallel loop will in [start, last]
		this->avail_chunk = (trip % chunk_size) == 0 ? trip / chunk_size : trip / chunk_size + 1;

		// the lb and ub of the first chunk (dynamic scheduling only)
		if (step > 0) {
			this->lb = start_point;
			this->ub = (this->lb + (chunk_size-1)*step) <= this->last ? (this->lb + (chunk_size-1)*step) : this->last;
		} else {
			this->ub = start_point;
			this->lb = (this->ub + (chunk_size-1)*step) >= this->trip ? (this->ub + (chunk_size-1)*step) : this->last;
		}

		// assign the start chunk of each thread (static scheduling only)
		for (int t = 0; t < THREAD_NUM; t++) {
			this->per_thread_start_point.push_back(start_point + (chunk_size*step) * t);
		}
#ifdef DEBUG
		cout << "total " << avail_chunk << " chunks can be assigned" << endl;
#endif 
	}
	// for dynamic, the chunk will be assigned in a FIFO order
	// we always print out the range of the next available chunk
	string getCurrentChunkRange() {
		int nextlb = 0, nextub = 0;
		if (step > 0) {
			nextlb = this->lb;
			nextub = this->ub;
		} else {
			nextub = this->ub;
			nextlb = this->lb;
		}
		return "[" + to_string(nextlb) + ", " + to_string(nextub) + "]";
	}

	// for static, the chunk has already been preassigned
	// we print out the range of the next chunk belongs to the given tid
	string getCurrentStaticChunkRange(unsigned tid) {
		int nextlb = 0, nextub = 0;
		if (step > 0) {
			nextlb = this->per_thread_start_point[tid];
			nextub = (nextlb + (chunk_size-1)*step) < this->last ? (nextlb + (chunk_size-1)*step) : this->last;
		} else {
			nextub = this->per_thread_start_point[tid];
			nextlb = (nextub + (chunk_size-1)*step) > this->last ? (nextub + (chunk_size-1)*step) : this->last;
		}
		return "[" + to_string(nextlb) + ", " + to_string(nextub) + "]";
	}

	// for dynamic, the next chunk is available if there avail_chunk is not 0
	bool hasNextChunk(bool isStatic) {
		bool ret = false;
		if (isStatic) {
			// static
			// hasNextChunk return false if hasNextStaticChunk(tid) return false for all tid
			for (int tid = 0; tid < THREAD_NUM; tid++) {
				ret |= hasNextStaticChunk(tid);
			}
			return ret; 
		}
		// dynamic
		// the lb of the next chunk exceed the parallel loop bound
		if (step > 0)
			return this->lb <= this->last;
		else
			return this->ub >= this->last;
	}
	// for static, the next chunk is available for a given tid if its next
	// available chunk still in the loop range
	bool hasNextStaticChunk(unsigned tid) {
		if (step > 0)
			return this->per_thread_start_point[tid] <= this->last;
		else
			return this->per_thread_start_point[tid] >= this->last;
	}
	// for dynamic
	Chunk getNextChunk(unsigned tid) {
		// assign the current lb, ub to thread tid and update the next chunk
		Chunk curr = make_pair(this->lb, this->ub);
		if (step > 0) {
			this->lb = this->ub + step;
			this->ub = (this->lb + (chunk_size-1)*step) < this->last ? (this->lb + (chunk_size-1)*step) : this->last;
		} else {
			this->ub = this->lb + step;
			this->lb = (this->ub + (chunk_size-1)*step) > this->last ? (this->ub + (chunk_size-1)*step) : this->last;
		}
		this->avail_chunk -= 1;
#ifdef DEBUG
		cout << "[DYNAMIC] return chunk [" << curr.first << "," << curr.second << "]" << endl;
#endif 
		return curr;
	}
	// for static
	Chunk getNextStaticChunk(unsigned tid) {
		int retlb = 0, retub = 0;
		if (step > 0) {
			retlb = this->per_thread_start_point[tid];
			retub = (retlb + (chunk_size-1)*step) < this->last ? (retlb + (chunk_size-1)*step) : this->last;
		} else {
			retub = this->per_thread_start_point[tid];
			retlb = (retub + (chunk_size-1)*step) > this->last ? (retub + (chunk_size-1)*step) : this->last;
		}
		Chunk curr = make_pair(retlb, retub);
		this->per_thread_start_point[tid] += (chunk_size*THREAD_NUM)*step;
#ifdef DEBUG
		cout << "[STATIC] return chunk [" << curr.first << "," << curr.second << "]" << endl;
#endif 
		return curr;
	}
	// the following functions derives the iteration info
	// inside the parallel region.
	// works only for static scheduling
	unsigned getStaticTid(int i) {
		return ((i-start)/step) / chunk_size - (floor(((i-start)/step) / (chunk_size * THREAD_NUM)) * THREAD_NUM);
	}
	// return the thread local chunk id
	int getStaticChunkID(int i) {
		return floor(((i-start)/step) / (chunk_size * THREAD_NUM));
	}

	int getStaticThreadLocalPos(int i) {
		return ((i-start)/step) % chunk_size;
	}
	// given i, we need to set the lb, ub and ptlb info
	// after set, each query will get a chunk after the 
	// given iteration
	void setStartPoint(int i) {
		int start_cid = getStaticChunkID(i);
		int start_tid = getStaticTid(i);
		// set the start point for static scheduling
		for (int t = 0; t < THREAD_NUM; t++) {
#ifdef DEBUG
			cout << "chunk at thread " << t << " starts from " << (this->per_thread_start_point[t] + (start_cid*chunk_size*THREAD_NUM)*step) << endl;
#endif
			this->per_thread_start_point[t] = this->per_thread_start_point[t] + (start_cid * chunk_size*THREAD_NUM)*step;
		}
		// set the start point for dynamic scheduling
		if (step > 0) {
			this->lb = start+((i-start)/(step*chunk_size))*step*chunk_size;
			this->ub = (this->lb + (chunk_size-1)*step) < this->last ? (this->lb + (chunk_size-1)*step) : this->last;
		} else {
			this->ub = start+((i-start)/(step*chunk_size))*step*chunk_size;
			this->lb = (this->ub + (chunk_size-1)*step) > this->last ? (this->ub + (chunk_size-1)*step) : this->last;
		}

		this->avail_chunk -= (start_cid * THREAD_NUM);
		// static
		// the start point of t in range [0, THREAD_NUM-1] will be 
		// [ this->ptlb[t] + start_chunk_pos, (this->ptlb[t] + chunk_size - 1) < this->trip ? (this->ptlb[t] + chunk_size - 1) : this->trip - 1;]

		// dynamic
		// the start point of the given t will be
		// [ i, this->ub]
		// the start point of other threads need to be generated manually

	}

	Chunk getStaticStartChunk(int i, int tid) {
		int start_chunk_pos = getStaticThreadLocalPos(i);
		int retlb = 0, retub = 0;
		if (step > 0) {
			retlb = per_thread_start_point[tid]+start_chunk_pos*step;
			retub = (per_thread_start_point[tid] + (chunk_size-1)*step) < this->last ? (per_thread_start_point[tid] + (chunk_size-1)*step) : this->last;
		} else {
			retub = per_thread_start_point[tid]+start_chunk_pos*step;
			retlb = (per_thread_start_point[tid] + (chunk_size-1)*step) > this->last ? (per_thread_start_point[tid] + (chunk_size-1)*step) : this->last;
		}
		Chunk curr = make_pair(retlb, retub);
#ifdef DEBUG
		cout << "[STATIC] return start chunk [" << curr.first << "," << curr.second << "]" << endl;
#endif
		this->per_thread_start_point[tid] += (chunk_size*THREAD_NUM)*step;
		return curr;
	}

	Chunk getStartChunk(int i) {
		int retlb = 0, retub = 0;
		if (step > 0) {
			retlb = i;
			if ((i - start) % (chunk_size * step) != 0) {
				retlb = ((i - start) / (chunk_size * step)) * (chunk_size * step);
			}
			retub = this->ub;
			// this->lb = this->ub + step;
			// this->ub = (this->lb + (chunk_size-1)*step) < this->last ? (this->lb + (chunk_size-1)*step) : this->last;
		} else {
			retub = i;
			if ((i - start) % (chunk_size * step) != 0) {
				retub = ((start - i) / (chunk_size * step)) * (chunk_size * step);
			}
			retlb = this->lb;
			// this->ub = this->lb + step;
			// this->lb = (this->ub + (chunk_size-1)*step) > this->last ? (this->ub + (chunk_size-1)*step) : this->last;
		}
		Chunk curr = make_pair(retlb, retub);
#ifdef DEBUG
		cout << "[DYNAMIC] return start chunk [" << curr.first << "," << curr.second << "]" << endl;
#endif
		return curr;
	}

	void getNextKChunksFrom(int k, Chunk c, vector<Chunk> &next_chunks) {
		int orig_lb = c.first, orig_ub = c.second;
		if (step > 0) {
			if ((c.first - start) % (chunk_size * step) != 0) {
				orig_lb = ((c.first - start) / (chunk_size * step)) * (chunk_size * step);
			}
		} else {
			if ((c.second - start) % (chunk_size * step) != 0) {
				orig_ub = ((start - c.second) / (chunk_size * step)) * (chunk_size * step);
			}
		}
		int i = 1;
		int lb = 0, ub = 0;
		while (i <= k) {
			if (step > 0) {
				lb = orig_lb + (chunk_size*step*i);
				if (lb > this->last)
					break;
				ub = lb + (chunk_size-1)*step <= this->last ? lb + (chunk_size-1)*step : this->last;
				next_chunks.push_back(make_pair(lb, ub));
			} else {
				ub = orig_ub + (chunk_size*step*i);
				if (ub < this->last)
					break;
				lb = ub + (chunk_size-1)*step >= this->last ? ub + (chunk_size-1)*step : this->last;
				next_chunks.push_back(make_pair(lb, ub));

			}
#ifdef DEBUG
			cout << "[DYNAMIC] push next chunk [" << lb << "," << ub << "] into the next_chunks list" << endl;
			cout << "[DYNAMIC] current chunk [" << this->lb << "," << this->ub << "]" << endl;
#endif
			i++;
		}
	}

	void getPrevKChunksFrom(int k, Chunk c, vector<Chunk> &prev_chunks) {
		int orig_lb = c.first, orig_ub = c.second;
		if (step > 0) {
			if ((c.first - start) % (chunk_size * step) != 0) {
				orig_lb = ((c.first - start) / (chunk_size * step)) * (chunk_size * step);
			}
		} else {
			if ((c.second - start) % (chunk_size * step) != 0) {
				orig_ub = ((start - c.second) / (chunk_size * step)) * (chunk_size * step);
			}
		}
		int i = 1;
		int lb = 0, ub = 0;
		while (i <= k) {
			if (step > 0) {
				lb = orig_lb - (chunk_size*step*i);
				if (lb < this->start)
					break;
				ub = (lb + (chunk_size-1)*step);
				prev_chunks.push_back(make_pair(lb, ub));
			} else {
				ub = orig_ub - (chunk_size*step*i);
				if (ub > this->start)
					break;
				lb = (ub + (chunk_size-1)*step);
				prev_chunks.push_back(make_pair(lb, ub));
			}

#ifdef DEBUG
			cout << "[DYNAMIC] push prev chunk [" << lb << "," << ub << "] into the prev_chunks list" << endl;
#endif
			i++;
		}
	}

	void moveToNextChunk() {
		if (step > 0) {
			this->lb = this->ub + step;
			this->ub = (this->lb + (chunk_size-1)*step) < this->last ? (this->lb + (chunk_size-1)*step) : this->last;
		} else {
			this->ub = this->lb + step;
			this->lb = (this->ub + (chunk_size-1)*step) > this->last ? (this->ub + (chunk_size-1)*step) : this->last;
		}
#ifdef DEBUG
		cout << "[DYNAMIC] move to chunk [" << this->lb << "," << this->ub << "]" << endl;
#endif
	}

	void moveToNextStaticChunk(int tid) {
		if (step > 0) {
			if ((per_thread_start_point[tid] + chunk_size*step) <= this->last) {
				per_thread_start_point[tid] += chunk_size*step;
			}
		} else {
			this->ub = this->lb + step;
			this->lb = (this->ub + (chunk_size-1)*step) >= this->last ? (this->ub + (chunk_size-1)*step) : this->last;
			if ((per_thread_start_point[tid] + chunk_size*step) >= this->last) {
				per_thread_start_point[tid] += chunk_size*step;
			}
		}
#ifdef DEBUG
		cout << "[STATIC] move to chunk [" << this->lb << "," << this->ub << "]" << endl;
#endif	
	}
};

class Progress {
public:
	string ref;
	Chunk chunk;
	vector<int> iteration;
	Progress() { }
	Progress(string ref, vector<int> iteration, Chunk c) {
		this->ref = ref;
		this->iteration = iteration;
		this->chunk = c;
	}
	string getIteration() {
		string ret = "(";
		for (unsigned i = 0; i < this->iteration.size(); i++) {
			ret += to_string(this->iteration[i]);
			if (i != this->iteration.size() - 1)
				ret += ",";
			}
		ret += ")";
		return ret;
	}
	string getReference() {
		return this->ref;
	}
	void increment(string ref, vector<int> iteration) {
		this->ref = ref;
		this->iteration = iteration;
	}
	void increment(string ref) {
		this->ref = ref;
	}
	bool isInBound() {
		assert(this->iteration[0] >= chunk.first);
		return this->iteration[0] <= chunk.second;
	}
    string toString() {
        stringstream ss;
        ss << "(" << this->ref << ",";
        ss << getIteration();
        ss << ")";
        return ss.str();
    }
};
#endif
#if 0
void pluss_reset();
void pluss_init();
void pluss_per_thread_init(thread::id tid);
// void pluss_ref_access(string array, vector<int> iteration);
void pluss_ref_access(string array, unsigned addr);
void pluss_sample_access(string array, vector<int> subscripts, Iteration *access, unordered_set<Iteration, IterationHasher> &samples);
void pluss_sample_access_parallel(string, vector<int>, Iteration *, unordered_set<Iteration, IterationHasher> &, thread::id);
// void pluss_access(string array, vector<unsigned> iteration, unordered_set<Iteration> &samples);
void pluss_parallel_access(unsigned addr);
void pluss_access(unsigned addr);
void pluss_parallel_bypass(unsigned skip_cnt);
void pluss_bypass(unsigned skip_cnt);
void pluss_per_thread_bypass(thread::id, unsigned);
void pluss_histogram_update(unsigned reuse, unsigned cnt);
void pluss_AET();
void pluss_print_histogram();
void pluss_print_mrc();
#endif

void pluss_timer_start();
void pluss_timer_stop();
void pluss_timer_print();
# ifndef PLUSS_CYCLE_ACCURATE_TIMER
double pluss_timer_return();
# else
unsigned long long int pluss_timer_return();
# endif
// void pluss_terminate();

#ifdef __cplusplus
}
#endif

#endif // end of PLUSS_H
