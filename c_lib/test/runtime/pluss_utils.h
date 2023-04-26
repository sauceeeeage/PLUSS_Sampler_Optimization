#ifndef PLUSS_UTILS_H
#define PLUSS_UTILS_H
#include <assert.h>
#include <cmath>
#include <array>
#include <vector>
#include <string>
#include <thread>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
// For parallel OPT
#include <shared_mutex>
#include <atomic>
// For CRI
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h> // cdf of gaussian distribution


typedef std::unordered_map<long, double> Histogram;

extern Histogram _RIHist;
extern std::unordered_map<long, std::atomic<double>> _ParRIHist;
extern std::map<uint64_t, double> _MRC;
extern std::shared_mutex _m;

extern std::array<Histogram, THREAD_NUM> _NoSharePRI; 
extern std::array<std::unordered_map<int, Histogram>, THREAD_NUM> _SharePRI;


namespace std {
typedef pair<int, int> Chunk;
class Iteration {
public: 
	string name;
	vector<int> ivs;
	int pidx;
	int priority;
	unsigned cid; // chunk this sample locates
	unsigned tid; // thread this sample belongs
	unsigned pos; // thread local position
	bool isParallelIteration;

	Iteration() {}
	~Iteration() {}
	Iteration(string ref, vector<int> iter, int start=0, int step=1, bool isInParallelRegion=false, int parallel_iter_idx=0) {
		name = ref;
		ivs = iter;
		pidx = parallel_iter_idx;
		// works when openmp static scheduling
		if (isInParallelRegion) {
			cid = floor(((iter[parallel_iter_idx]-start)/step) / (CHUNK_SIZE * THREAD_NUM));
			tid = ((iter[parallel_iter_idx]-start)/step) / CHUNK_SIZE - floor(((iter[parallel_iter_idx]-start)/step) / (CHUNK_SIZE * THREAD_NUM)) * THREAD_NUM;
			pos = ((iter[parallel_iter_idx]-start)/step) % CHUNK_SIZE;
		}
		isParallelIteration = isInParallelRegion;
		priority = 1;
	}
	Iteration(string ref, int priority, vector<int> iter, int start=0, int step=1, bool isInParallelRegion=false, int parallel_iter_idx=0) {
		name = ref;
		ivs = iter;
		pidx = parallel_iter_idx;
		// works when openmp static scheduling
		if (isInParallelRegion) {
			cid = floor(((iter[parallel_iter_idx]-start)/step) / (CHUNK_SIZE * THREAD_NUM));
			tid = ((iter[parallel_iter_idx]-start)/step) / CHUNK_SIZE - floor(((iter[parallel_iter_idx]-start)/step) / (CHUNK_SIZE * THREAD_NUM)) * THREAD_NUM;
			pos = ((iter[parallel_iter_idx]-start)/step) % CHUNK_SIZE;
		}
		isParallelIteration = isInParallelRegion;
		priority = priority; // mark the topology order
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
	string toAddrString() {
		string s = name;
		vector<int>::iterator it;
		for (it = ivs.begin(); it != ivs.end(); ++it) {
			s += to_string(*it) + "-";
		}
		return s;
	}
	int compare(Iteration other) {
		if (ivs.size() != other.ivs.size()) { return 2; }

		// first assumes the two come from the same reference, we check their
		// iteration vector

		/* compare the chunk id */ ;
		if (isParallelIteration) {
			if (cid < other.cid) {
				return -1;
			} else if (cid > other.cid) {
				return 1;
			}
		// compare of pos and tid only works for uniform interleaving + static 
		// the same chunk. compare the thread local position.
			if (pos < other.pos) {
				return -1;
			} else if (pos > other.pos) {
				return 1;
			}
		}
		/* compare the loop induction variables */ ;
		unsigned ivs_idx = 0;
		if (isParallelIteration) {
			// we first compare those idx outside the parallel region
			while (ivs_idx < ivs.size() && ivs_idx < this->pidx) {
				if (ivs[ivs_idx] < other.ivs[ivs_idx])
					return -1;
				else if (ivs[ivs_idx] > other.ivs[ivs_idx])
					return 1;
				ivs_idx++;
			}
			if (ivs_idx == this->pidx)
				ivs_idx++; // skip the parallel iteration
			// compare the iteration inside the parallel region
			while (ivs_idx < ivs.size()) {
				if (ivs[ivs_idx] < other.ivs[ivs_idx])
					return -1;
				else if (ivs[ivs_idx] > other.ivs[ivs_idx])
					return 1;
				ivs_idx++;
			}
			/* the same thread local position. compare the thread */ ;
			if (tid < other.tid) {
				return -1;
			} else if (tid > other.tid) {
				return 1;
			}
		} else {
			// the access is in sequential region, 
			// just compare the loop induction variables
			while (ivs_idx < ivs.size()) {
				if (ivs[ivs_idx] < other.ivs[ivs_idx])
					return -1;
				else if (ivs[ivs_idx] > other.ivs[ivs_idx])
					return 1;
				ivs_idx++;
			}
		}

		// priority indicates the topology order of the references, higher means
		// this reference is accessed early
		if (priority > other.priority)
			return -1;
		else if (priority < other.priority)
			return 1;
		
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
struct  IterationComp
{
	bool operator() (Iteration const& lhs, Iteration const& rhs)
	{

		// it is possible that two accesses has different itration vector size


		// this two access from the same iteration vector size;
		// - they could come from the same loop
		//		need to compare the iteration vector first, if they both in 
		// 		the same iteration, then check the access order
		// - they come from different loop with the same loop nest size
		//		compare the priority can distinguish

		// first assumes the two come from the same reference, we check their
		// iteration vector

		/* compare the chunk id */ ;
		if (lhs.isParallelIteration) {
			if (lhs.cid < rhs.cid) {
				return false;
			} else if (lhs.cid > rhs.cid ) {
				return true;
			}
		// compare of pos and tid only works for uniform interleaving + static 
		// the same chunk. compare the thread local position.
			if (lhs.pos < rhs.pos) {
				return false;
			} else if (lhs.pos > rhs.pos) {
				return true;
			}
		}
		/* the same chunk. compare the rest loop induction variables */
		unsigned ivs_idx = 0;
		unsigned common_ivs_size = min(lhs.ivs.size(), rhs.ivs.size());
		if (lhs.isParallelIteration) {

			// we first compare those idx outside the parallel region
			while (ivs_idx < common_ivs_size && ivs_idx < lhs.pidx) {
				if (lhs.ivs[ivs_idx] < rhs.ivs[ivs_idx])
					return false;
				else if (lhs.ivs[ivs_idx] > rhs.ivs[ivs_idx])
					return true;
				ivs_idx++;
			}
			if (ivs_idx == lhs.pidx)
				ivs_idx++; // skip the parallel iteration

			// compare the iteration inside the parallel region
			while (ivs_idx < common_ivs_size) {
				if (lhs.ivs[ivs_idx] < rhs.ivs[ivs_idx])
					return false;
				else if (lhs.ivs[ivs_idx] > rhs.ivs[ivs_idx])
					return true;
				ivs_idx++;
			}
			/* the same thread local position. compare the thread */ ;
			if (lhs.tid < rhs.tid) {
				return false;
			} else if (lhs.tid > rhs.tid) {
				return true;
			}
		} else {
			while (ivs_idx != common_ivs_size) {
				if (lhs.ivs[ivs_idx] < rhs.ivs[ivs_idx])
					return false;
				else if (lhs.ivs[ivs_idx] > rhs.ivs[ivs_idx])
					return true;
				ivs_idx++;
			}
		}

		// when reaching here
		// it is possible that
		// 1) they are in the same loop and the same iteration
		// 2) they are in the same loop and their common loop iteartion are 
		// the same
		// 3) they are in different loop but same iteration vector size
		// In both three cases, we can distinguish the order by comparing their 
		// access order (priority)
		//
		// priority indicates the topology order of the references, higher means
		// this reference is accessed early
		if (lhs.priority > rhs.priority)
			return false;
		else if (lhs.priority < rhs.priority)
			return true;

		/* all equal, these two samples are equal */ ;
		return false;
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
	array<int, THREAD_NUM> per_thread_start_point; // use for static sheduling only
private:
	void init() {
		// the range of the parallel loop will in [start, last]
		this->avail_chunk = (this->trip % this->chunk_size) == 0 ? this->trip / this->chunk_size : this->trip / this->chunk_size + 1;
		// the lb and ub of the first chunk (dynamic scheduling only)
		if (this->step > 0) {
			this->lb = this->start;
			this->ub = (this->lb + (this->chunk_size-1)*this->step) <= this->last ? (this->lb + (this->chunk_size-1)*this->step) : this->last;
		} else {
			this->ub = this->start;
			this->lb = (this->ub + (this->chunk_size-1)*this->step) >= this->trip ? (this->ub + (this->chunk_size-1)*this->step) : this->last;
		}

		// assign the start chunk of each thread (static scheduling only)
		for (int t = 0; t < THREAD_NUM; t++) {
			this->per_thread_start_point[t] = (this->start + (this->chunk_size*this->step) * t);
		}
#ifdef DEBUG
		cout << "total " << this->avail_chunk << " chunks can be assigned" << endl;
#endif 
	}
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

		init();
	}
	void reset() {
		init();
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

/* util func */
inline long _polybench_to_highest_power_of_two(long x)
{
	// check for the set bits
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
	return x ^ (x >> 1);
}
inline void _pluss_histogram_update(Histogram& histogram, long reuse, double cnt, bool in_log_format=true)
{
	if (reuse > 0 && in_log_format)
		reuse = _polybench_to_highest_power_of_two(reuse);
	if (histogram.find(reuse) != histogram.end()) {
		histogram[reuse] += cnt;
	} else {
		histogram[reuse] = cnt;
	}
}
inline void _pluss_histogram_print(string title, Histogram& histogram)
{
	map<long, double> tmp;
	cout << title << endl;
	double sum = 0.;
	for (auto iter = histogram.begin(); iter != histogram.end(); ++iter) {
		sum += iter->second;
		tmp[iter->first] = iter->second;
	}
	for (auto iter = tmp.begin(); iter != tmp.end(); ++iter) {
		cout << iter->first << "," << iter->second <<  "," << iter->second/sum << endl;
	}
}
inline uint64_t distance_to(uint64_t x, uint64_t y)
{
	if (x > y)
		return x - y;
	return y - x;
}
/// Reset Histogram
inline void pluss_histogram_reset()
{
	if (!_RIHist.empty())
		_RIHist.clear();
	if (!_MRC.empty())
		_MRC.clear();
}
inline void pluss_parallel_histogram_reset(Histogram &histogram)
{
	histogram.clear();
}
/// Update Histogram
inline void pluss_histogram_update(long reuse, double cnt)
{
	_pluss_histogram_update(_RIHist, reuse, cnt);
}
inline void pluss_parallel_histogram_update(Histogram &histogram, long reuse, double cnt, bool inlog=false)
{
	_pluss_histogram_update(histogram, reuse, cnt, inlog);
}
/*
inline void pluss_parallel_histogram_update(long reuse, double cnt)
{
	if (reuse > 0)
		reuse = _polybench_to_highest_power_of_two(reuse);
	_m.lock_shared();
	if (_ParRIHist.find(reuse) != _ParRIHist.end()) {
		_m.unlock_shared();
		_ParRIHist[reuse].fetch_add(cnt);
	} else {
		_m.unlock_shared();
		_m.lock();
		_ParRIHist[reuse] += cnt;
		_m.unlock();
	}
}
*/
/// Print Histogram
inline void pluss_print_histogram()
{
	_pluss_histogram_print("Start to dump reuse time", _RIHist);
}

inline void pluss_parallel_print_histogram(Histogram &histogram)
{
	_pluss_histogram_print("Start to dump reuse time", histogram);
}
/// MRC Compute
inline void pluss_AET() {


	map<uint64_t, double> P;
	map<long, double> histogram;
	double total_num_RT = 0;
	long max_RT = 0;
	for (auto it = _RIHist.begin(); it != _RIHist.end(); ++it) {
		total_num_RT += it->second;
		histogram[it->first] = it->second;
		if (max_RT < it->first) {
			max_RT = it->first;
		}
	}
	double accumulate_num_RT = 0.;
	if (histogram.find(-1) != histogram.end())
		accumulate_num_RT = histogram[-1];
	for (auto it = histogram.rbegin(); it != histogram.rend(); ++it) {
		if (it->first == -1)
			break;
		P[it->first] = accumulate_num_RT / total_num_RT;
		accumulate_num_RT += it->second;
	}
	P[0] = 1.0;
	double sum_P = 0, MRC_pred = -1.0;
	uint64_t t = 0;
	uint64_t prev_t = 0;
	uint64_t cs = 2560 * 1024 / sizeof(double);
	for (uint64_t c = 0; c <= max_RT && c <= cs; c++) {
		while (sum_P < c && t <= max_RT) {
			if (P.find(t) != P.end()) {
				sum_P += P[t];
				prev_t = t;
			} else {
				sum_P += P[prev_t];
			}
			t++;
		}
		if (MRC_pred != -1.0) {
			_MRC[c] = P[prev_t];
		} else if (MRC_pred - P[prev_t] < 0.0001) {
			_MRC[c] = P[prev_t];
			MRC_pred = P[prev_t];
		}
	}
	return;
}
inline void pluss_parallel_AET(Histogram &histogram) {
	map<uint64_t, double> P;
	map<uint64_t, double> tmp;
	double total_num_RT = 0.;
	long max_RT = 0;
	for (auto it = histogram.begin(); it != histogram.end(); ++it) {
		tmp[it->first] = it->second;
		total_num_RT += it->second;
		if (max_RT < it->first) {
			max_RT = it->first;
		}
	}
	double accumulate_num_RT = 0.0;
	if (histogram.find(-1) != histogram.end())
		accumulate_num_RT = histogram[-1];
	for (auto it = tmp.rbegin(); it != tmp.rend(); ++it) {
		if (it->first == -1)
			break;
		P[it->first] = accumulate_num_RT / total_num_RT;
		accumulate_num_RT += it->second;
	}
	P[0] = 1.0;
	double sum_P = 0.0, MRC_pred = -1.0;
	uint64_t t = 0;
	uint64_t prev_t = 0;
	uint64_t cs = 2560 * 1024 / sizeof(double);
	for (uint64_t c = 0; c <= max_RT && c <= 327680; c++) {
		while (sum_P < c && t <= max_RT) {
			if (P.find(t) != P.end()) {
				sum_P += P[t];
				prev_t = t;
			} else {
				sum_P += P[prev_t];
			}
			t++;
		}
		if (MRC_pred != -1.0) {
			_MRC[c] = P[prev_t];
		} else if (MRC_pred - P[prev_t] < 0.0001) {
			_MRC[c] = P[prev_t];
			MRC_pred = P[prev_t];
		}
	}
	return;
}
/// Print MRC
inline void pluss_print_mrc()
{
	cout << "miss ratio" << endl;
	map<uint64_t, double>::iterator it1 = _MRC.begin();
	map<uint64_t, double>::iterator it2 = _MRC.begin();
    while(it1 != _MRC.end()) {
        while(1) {
        	map<uint64_t, double>::iterator it3 = it2;
            ++it3;
            if (it3 == _MRC.end()) {
                break;
            }
            if (it1->second - it3->second < 0.00001) {
                ++it2;
            } else {
                break;
            }
        }
        cout << it1->first << ", " << it1->second << endl;
        if (it1 != it2) {
            cout << it2->first << ", " << it2->second << endl;
        }
        it1 = ++it2;
        it2 = it1;
    }
    return;
#if 0
	if (cnt > 0) {
		cout << "max iteration traversed" << endl;
		cout << cnt << endl;
	}
#endif
}
/// Store MRC
inline void pluss_write_mrc_to_file(string file)
{
	ofstream mrc_file;
	mrc_file.open (file);
	mrc_file << "miss ratio" << endl;
	map<uint64_t, double>::iterator it1 = _MRC.begin();
	map<uint64_t, double>::iterator it2 = _MRC.begin();
    while(it1 != _MRC.end()) {
        while(1) {
        	map<uint64_t, double>::iterator it3 = it2;
            ++it3;
            if (it3 == _MRC.end()) {
                break;
            }
            if (it1->second - it3->second < 0.00001) {
                ++it2;
            } else {
                break;
            }
        }
        mrc_file << it1->first << ", " << it1->second << endl;
        if (it1 != it2) {
            mrc_file << it2->first << ", " << it2->second << endl;
        }
        it1 = ++it2;
        it2 = it1;
    }
	mrc_file.close();
}
/// CRI Model Util
inline void pluss_cri_histogram_reset()
{
	for (unsigned i = 0; i < _SharePRI.size(); i++) {
		_SharePRI[i].clear();
	}
	for (unsigned i = 0; i < _NoSharePRI.size(); i++) {
		_NoSharePRI[i].clear();
	}
}
inline void pluss_cri_noshare_histogram_update(int tid, long reuse, double count)
{
	_pluss_histogram_update(_NoSharePRI[tid], reuse, count);
}
inline void pluss_cri_share_histogram_update(int tid, int share_ratio, long reuse, double count)
{
	// if (reuse > 0)
	// 	reuse = _polybench_to_highest_power_of_two(reuse);
	if (_SharePRI[tid].find(share_ratio) != _SharePRI[tid].end()) {
		_pluss_histogram_update(_SharePRI[tid][share_ratio], reuse, count, false);
	} else {
		_SharePRI[tid][share_ratio][reuse] = count;
	}
}
inline void pluss_cri_noshare_print_histogram()
{
	Histogram noshare_rih_tmp;
	for (unsigned i = 0; i < _NoSharePRI.size(); i++) {
		for (auto iter = _NoSharePRI[i].begin(); iter != _NoSharePRI[i].end(); iter++) {
			_pluss_histogram_update(noshare_rih_tmp, iter->first, iter->second, false);
		}
	}
	_pluss_histogram_print("Start to dump noshare private reuse time", noshare_rih_tmp);

}
inline void pluss_cri_share_print_histogram()
{
	Histogram share_rih_tmp;
	for (unsigned i = 0; i < _SharePRI.size(); i++) {
		for (auto iter = _SharePRI[i].begin(); iter != _SharePRI[i].end(); ++iter) {
			for (auto ii = iter->second.begin(); ii != iter->second.end(); ++ii) {
				_pluss_histogram_update(share_rih_tmp, ii->first, ii->second, false);
			}
		}
	}
	_pluss_histogram_print("Start to dump share private reuse time", share_rih_tmp);
}
inline void pluss_pri_print_histogram()
{
	Histogram merged_dist;
	for (unsigned i = 0; i < _NoSharePRI.size(); i++) {
		for (auto entry : _NoSharePRI[i]) {
			if (merged_dist.find(entry.first) != merged_dist.end()) {
				merged_dist[entry.first] += entry.second;
			} else {
				merged_dist[entry.first] = entry.second;
			}
		}
	}
	for (unsigned i = 0; i < _SharePRI.size(); i++) {
		for (auto share_entry : _SharePRI[i]) {
			for (auto reuse_entry : share_entry.second) {
				if (merged_dist.find(reuse_entry.first) != merged_dist.end()) {
					merged_dist[reuse_entry.first] += reuse_entry.second;
				} else {
					merged_dist[reuse_entry.first] += reuse_entry.second;
				}
			}
		}
	}
	_pluss_histogram_print("Start to dump private reuse time", merged_dist);
}
/// CRI Model
inline void _pluss_cri_nbd(int thread_cnt, long n, Histogram &dist)
{
#if defined(DEBUG)
    cout << "Distribute thread-local Reuse: " << n << endl;
#endif
    double p = 1.0 / thread_cnt;
    if (n >= (4000. * (thread_cnt-1)) / thread_cnt) {
        // int i = (int)log2(n);
        // long bin = (long)(pow(2.0, i));
        dist[THREAD_NUM*n] = 1.0;
        return;
    }
    long k = 0;
    double nbd_prob = 0.0, prob_sum = 0.0;
    while (true) {
        nbd_prob = gsl_ran_negative_binomial_pdf(k, p, (double)n);
        prob_sum += nbd_prob;
        dist[k+n] = nbd_prob;
        if (prob_sum > 0.9999)
            break;
        k += 1;
    }
} // end of void pluss_cri_nbd()
inline void _pluss_cri_noshare_distribute(int thread_cnt=THREAD_NUM)
{
    Histogram dist;
	Histogram merged_dist;
	for (unsigned i = 0; i < _NoSharePRI.size(); i++) {
		for (auto entry : _NoSharePRI[i]) {
			if (merged_dist.find(entry.first) != merged_dist.end()) {
				merged_dist[entry.first] += entry.second;
			} else {
				merged_dist[entry.first] = entry.second;
			}
		}
	}
	for (auto entry : merged_dist) {
		if (entry.first < 0) {
			pluss_histogram_update(entry.first,entry.second);
			continue;
		}
		if (thread_cnt > 1) {
			_pluss_cri_nbd(thread_cnt, entry.first, dist);
			for (auto dist_entry : dist) {
				long ri_to_distribute = dist_entry.first;
				pluss_histogram_update(ri_to_distribute, entry.second * dist_entry.second);
			}
			dist.clear();
		} else {
			pluss_histogram_update(entry.first,entry.second);
		} // end of if(thread_cnt > 1)
	}
} // end of void pluss_cri_noshare_distribute()
inline void _pluss_cri_racetrack(int thread_cnt=THREAD_NUM)
{
	unordered_map<int, Histogram> merged_dist;
	for (unsigned i = 0; i < _SharePRI.size(); i++) {
		for (auto share_entry : _SharePRI[i]) {
			if (merged_dist.find(share_entry.first) != merged_dist.end()) {
				for (auto reuse_entry : share_entry.second) {
					if (merged_dist.find(reuse_entry.first) != merged_dist.end()) {
						merged_dist[share_entry.first][reuse_entry.first] += reuse_entry.second;
					} else {
						merged_dist[share_entry.first][reuse_entry.first] += reuse_entry.second;
					}
				}
			} else {
				for (auto reuse_entry : share_entry.second) {
					merged_dist[share_entry.first][reuse_entry.first] = reuse_entry.second;
				}
			}
		}
	}
	for (auto share_entry : merged_dist) {
		unordered_map<int, double> prob;
		Histogram dist;
		int i = 1;
		double prob_sum = 0.0, n = (double)share_entry.first;
		for (auto reuse_entry : share_entry.second) {
			// printf("share parameter n = %f\n", n);
			if (thread_cnt > 1) {

				_pluss_cri_nbd(thread_cnt, reuse_entry.first, dist);
				for (auto dist_entry: dist) {
					long ri_to_distribute = dist_entry.first;
					double cnt_to_distribute = reuse_entry.second * dist_entry.second;
#if defined(DEBUG)
					printf("new ri: 	%lu, %f*%f=%f\n", ri_to_distribute, reuse_entry.second, dist_entry.second, cnt_to_distribute);
#endif
					prob_sum = 0.0;
					i = 1;
					while (true) {
						if (pow(2.0, (double)i) > ri_to_distribute) { break; }
						prob[i] = pow(1-(pow(2.0, (double)i-1) / ri_to_distribute), n) - pow(1-(pow(2.0, (double)i) / ri_to_distribute), n);
						prob_sum += prob[i];
#if defined(DEBUG)
						printf("prob[2^%d <= ri < 2^%d] = %f (%f)\n", i-1, i, prob[i], prob[i]*cnt_to_distribute);
#endif
						i++;
						if (prob_sum == 1.0) { break; }
					} // end of while(true)
					if (prob_sum != 1.0) {
						prob[i-1] = 1 - prob_sum;
#if defined(DEBUG)
						printf("prob[ri >= 2^%d] = %f (%f)\n", i-1, prob[i], prob[i]*cnt_to_distribute);
#endif
					}
					for (auto bin : prob) {
						long new_ri = (long)pow(2.0, bin.first-1);
						pluss_histogram_update(new_ri, bin.second*cnt_to_distribute);
					}
					prob.clear();
					// pluss_histogram_update(ri_to_distribute, cnt_to_distribute);
				} // end of iterating dist
				dist.clear();
#if 0
				long ri_to_distribute = reuse_entry.first;
				double cnt_to_distribute = reuse_entry.second;
				printf("new ri: 	%lu, %f\n", ri_to_distribute, cnt_to_distribute);
				prob_sum = 0.0;
				i = 1;
				while (true) {
					if (pow(2.0, (double)i) > ri_to_distribute) { break; }
					prob[i] = pow(1-(pow(2.0, (double)i-1) / ri_to_distribute), n) - pow(1-(pow(2.0, (double)i) / ri_to_distribute), n);
					prob_sum += prob[i];
					printf("prob[2^%d <= ri < 2^%d] = %f (%f)\n", i-1, i, prob[i], prob[i]*cnt_to_distribute);
					i++;
					if (prob_sum == 1.0) { break; }
				} // end of while(true)
				if (prob_sum != 1.0) {
					prob[i-1] = 1 - prob_sum;
					printf("prob[ri >= 2^%d] = %f (%f)\n", i-1, prob[i], prob[i]*cnt_to_distribute);
				}
				for (auto bin : prob) {
					long new_ri = (long)pow(2.0, bin.first-1);
					pluss_histogram_update(new_ri, bin.second*cnt_to_distribute);
				}
				prob.clear();
#endif
			} else {
				pluss_histogram_update(reuse_entry.first,reuse_entry.second);
			}
		} // end of iterating all reuse entries
	} // end of iterating all type of reuses
} // end of void pluss_cri_racetrack()
#if 0
inline void _pluss_cri_racetrack(int thread_cnt=THREAD_NUM)
{
	unordered_map<int, Histogram> merged_dist;
	for (unsigned i = 0; i < _SharePRI.size(); i++) {
		for (auto share_entry : _SharePRI[i]) {
			if (merged_dist.find(share_entry.first) != merged_dist.end()) {
				for (auto reuse_entry : share_entry.second) {
					printf("merge ri: 	%lu, %f\n", reuse_entry.first, reuse_entry.second);
					if (merged_dist.find(reuse_entry.first) != merged_dist.end()) {
						merged_dist[share_entry.first][reuse_entry.first] += reuse_entry.second;
					} else {
						merged_dist[share_entry.first][reuse_entry.first] += reuse_entry.second;
					}
				}
			} else {
				for (auto reuse_entry : share_entry.second) {
					printf("merge ri: 	%lu, %f\n", reuse_entry.first, reuse_entry.second);
					merged_dist[share_entry.first][reuse_entry.first] = reuse_entry.second;
				}
			}
		}
	}
	for (auto share_entry : merged_dist) {
		unordered_map<int, double> prob;
		Histogram dist;
		int i = 1;
		double prob_sum = 0.0, n = (double)share_entry.first;
		for (auto reuse_entry : share_entry.second) {
			// printf("share parameter n = %f\n", n);
			if (thread_cnt > 1) {

				long ri_to_distribute = reuse_entry.first;
				double cnt_to_distribute = reuse_entry.second;
				
// #if defined(DEBUG)
				printf("new ri: 	%lu, %f\n", ri_to_distribute, cnt_to_distribute);
// #endif
				prob_sum = 0.0;
				double p = 1.0 / ri_to_distribute;
				i = 1;
				while (true) {
					if (pow(2.0, (double)i) > ri_to_distribute) { break; }
					unsigned int b = (unsigned int)pow(2.0, (double)i);
					unsigned int a = (unsigned int)pow(2.0, (double)(i-1));
					prob[i] = gsl_cdf_geometric_P(b, p) - gsl_cdf_geometric_P(a, p);
					prob_sum += prob[i];
// #if defined(DEBUG)
					printf("prob[2^%d <= ri < 2^%d] = %f (%f)\n", i-1, i, prob[i], prob[i]*cnt_to_distribute);
// #endif
					i++;
					if (prob_sum == 1.0) { break; }
				} // end of while(true)
				if (prob_sum != 1.0) {
					prob[i-1] = 1 - prob_sum;
// #if defined(DEBUG)
					printf("prob[ri >= 2^%d] = %f (%f)\n", i-1, prob[i], prob[i]*cnt_to_distribute);
// #endif
				}
				for (auto bin : prob) {
					long new_ri = (long)pow(2.0, bin.first-1);
					pluss_histogram_update(new_ri, bin.second*cnt_to_distribute);
				}
				prob.clear();
				// pluss_histogram_update(ri_to_distribute, cnt_to_distribute);
			} else {
				pluss_histogram_update(reuse_entry.first,reuse_entry.second);
			}
		} // end of iterating all reuse entries
	} // end of iterating all type of reuses
} // end of void pluss_cri_racetrack()
#endif
inline void pluss_cri_distribute(int thread_cnt)
{
	_pluss_cri_noshare_distribute(thread_cnt);
	_pluss_cri_racetrack(thread_cnt);
}
} // namespace std

#endif