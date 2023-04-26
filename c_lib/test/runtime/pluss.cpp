#include "pluss.h"
#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <time.h>
#include <sys/time.h>

/* Total LLC cache size. By default 32+MB.. */
#ifndef POLYBENCH_CACHE_SIZE_KB
# define POLYBENCH_CACHE_SIZE_KB 2560
#endif

#define CACHE_MASK  ~((1UL << 6) - 1)


/* Timer code (gettimeofday). */
double pluss_t_start, pluss_t_end;
/* Timer code (RDTSC). */
unsigned long long int pluss_c_start, pluss_c_end;

#if 0
/* locality analysis global variable */
unsigned cnt;
unordered_map<string, unordered_map<unsigned, unsigned>> PerArrayLAT;
unordered_map<thread::id, unordered_map<string, unordered_map<string, unsigned>>> PerThreadPerArrayLAT;
unordered_map<unsigned, unsigned> LAT;
#endif
// unordered_map<uint64_t, uint64_t> RIHist;
// unordered_map<unsigned, atomic<unsigned>> RIHistTmp;
// unordered_map<unsigned, unsigned> RIHistTmp;
// map<uint64_t, double> MRC;
#if 0
// access time -> pointer points to an Iteration object
unordered_map<unsigned, Iteration *> LATSampleIterMap;
unordered_map<thread::id, unordered_map<unsigned, Iteration *>> PerThreadLATSampleIterMap;
unordered_map<thread::id, unsigned> counter;

shared_mutex _m;
#endif
/* Locality anlaysis utils */
unsigned _polybench_to_highest_power_of_two(unsigned x);
void _polybench_flush_cache();

/* Timer */
static
double rtclock()
{
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, NULL);
    if (stat != 0)
      printf ("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}


#ifdef PLUSS_CYCLE_ACCURATE_TIMER
static
unsigned long long int rdtsc()
{
  unsigned long long int ret = 0;
  unsigned int cycles_lo;
  unsigned int cycles_hi;
  __asm__ volatile ("RDTSC" : "=a" (cycles_lo), "=d" (cycles_hi));
  ret = (unsigned long long int)cycles_hi << 32 | cycles_lo;

  return ret;
}
#endif

void _polybench_flush_cache()
{
	int cs = POLYBENCH_CACHE_SIZE_KB * 1024 / sizeof(double);
	double* flush = (double*) calloc (cs, sizeof(double));
	int i;
	double tmp = 0.0;
	for (i = 0; i < cs; i++)
	tmp += flush[i];
	assert (tmp <= 10.0);
	free (flush);
}



/* Timer */
void pluss_timer_start()
{
  _polybench_flush_cache ();
#ifndef PLUSS_CYCLE_ACCURATE_TIMER
  pluss_t_start = rtclock ();
#else
  pluss_c_start = rdtsc ();
#endif
}

void pluss_timer_stop()
{
#ifndef PLUSS_CYCLE_ACCURATE_TIMER
  pluss_t_end = rtclock ();
#else
  pluss_c_end = rdtsc ();
#endif
}

void pluss_timer_print()
{
# ifndef PLUSS_CYCLE_ACCURATE_TIMER
      printf ("%0.6f\n", pluss_t_end - pluss_t_start);
# else
      printf ("%Ld\n", pluss_c_end - pluss_c_start);
# endif
}

# ifndef PLUSS_CYCLE_ACCURATE_TIMER
double pluss_timer_return()
{
	return pluss_t_end - pluss_t_start;
}
# else
unsigned long long int pluss_timer_return()
{
	return pluss_c_end - pluss_c_start;
}
# endif

#if 0
void pluss_parallel_histogram_update(unsigned reuse, unsigned cnt)
{
    unsigned bin = _polybench_to_highest_power_of_two(reuse);
    _m.lock_shared();
    if (RIHistTmp.find(bin) != RIHistTmp.end()) {
        _m.unlock_shared();
        RIHistTmp[bin].fetch_add(cnt);
    } else {
        _m.unlock_shared();
        _m.lock();
        RIHistTmp[bin] += cnt;
        _m.unlock();
    }
}
void pluss_reset()
{
	cnt = 0;
	for (auto pair : PerArrayLAT) {
		pair.second.clear();
	}
	PerArrayLAT.clear();
}
void pluss_init()
{
	cnt = 0UL;;
}
void pluss_per_thread_init(thread::id tid)
{
    counter[tid] = 0;
}

#if 0
void pluss_ref_access(string array, vector<int> iteration)
{
	string addr = "";
	for (unsigned i = 0; i < iteration.size(); i++) {
		if (i == iteration.size()-1)
			addr += to_string(iteration[i]);
		else
			addr += to_string(iteration[i])+"_";
	}
	if (PerArrayLAT.find(array) != PerArrayLAT.end()) {
		// could have reuse
		if (PerArrayLAT[array].find(addr) != PerArrayLAT[array].end()) {
			unsigned reuse = cnt - PerArrayLAT[array][addr];
			pluss_histogram_update(reuse, 1);
		}	
	}
	PerArrayLAT[array][addr] = cnt;
	cnt++;
}
#endif
void pluss_ref_access(string array, unsigned addr)
{
	if (PerArrayLAT.find(array) != PerArrayLAT.end()) {
		// could have reuse
		if (PerArrayLAT[array].find(addr) != PerArrayLAT[array].end()) {
			unsigned reuse = cnt - PerArrayLAT[array][addr];
			pluss_histogram_update(reuse, 1);
		}	
	}
	PerArrayLAT[array][addr] = cnt;
	cnt++;
}


void pluss_sample_access(string array, vector<int> subscripts, Iteration *access, unordered_set<Iteration, IterationHasher> &samples)
{
#if defined(DEBUG)
	// if (access)
	// 	cout << access->toString() << endl;
#endif
	unsigned addr = 0;
	for (unsigned i = 0; i < subscripts.size(); i++) {
        addr += subscripts[i];
	}

	if (PerArrayLAT.find(array) != PerArrayLAT.end()) {
		if (PerArrayLAT[array].find(addr) != PerArrayLAT[array].end()) {
			unsigned reuse = cnt - PerArrayLAT[array][addr];
			pluss_histogram_update(reuse, 1);
			Iteration *src = LATSampleIterMap[PerArrayLAT[array][addr]];
#if defined(DEBUG)
			string s = array + " (";
			vector<int>::iterator it = subscripts.begin();
			for (; it != subscripts.end(); ++it) {
				s += to_string(*it) + ",";
			}
			s.pop_back();
			s += ")";
			cout << reuse << " " <<  src->toString() << " -> " << s << endl;
#endif
			PerArrayLAT[array].erase(addr);
			samples.erase(*src);
		}
	}

	if (access && (samples.find(*access) != samples.end())) {
		PerArrayLAT[array][addr] = cnt;
		LATSampleIterMap[cnt] = access;
		// cout << "Update LAT and LATSampleIterMap for sample: " <<  access->toString() << " " << cnt << endl;
	}
#if defined(DEBUG)
	if (samples.empty())
		cout << "sampler terminate" << endl;
#endif
	cnt++;
}

void pluss_sample_access_parallel(string array, vector<int> subscripts, Iteration *access, unordered_set<Iteration, IterationHasher> &samples, thread::id tid)
{
    _m.lock();
#if defined(DEBUG)
	// if (access)
	// 	cout << access->toString() << endl;
#endif
	string addr = "";
	for (unsigned i = 0; i < subscripts.size(); i++) {
		if (i == subscripts.size()-1)
			addr += to_string(subscripts[i]);
		else
			addr += to_string(subscripts[i])+"_";
	}

	if (PerThreadPerArrayLAT[tid].find(array) != PerThreadPerArrayLAT[tid].end()) {
		if (PerThreadPerArrayLAT[tid][array].find(addr) != PerThreadPerArrayLAT[tid][array].end()) {
			unsigned reuse = counter[tid] - PerThreadPerArrayLAT[tid][array][addr];
			pluss_histogram_update(reuse, 1);
			Iteration *src = PerThreadLATSampleIterMap[tid][PerThreadPerArrayLAT[tid][array][addr]];
#if defined(DEBUG)
			string s = array + " (";
			vector<int>::iterator it = subscripts.begin();
			for (; it != subscripts.end(); ++it) {
				s += to_string(*it) + ",";
			}
			s.pop_back();
			s += ")";
			cout << reuse << " " <<  src->toString() << " -> " << s << endl;
#endif
			PerThreadPerArrayLAT[tid][array].erase(addr);
			samples.erase(*src);
		}
	}

	if (access && (samples.find(*access) != samples.end())) {
		PerThreadPerArrayLAT[tid][array][addr] = counter[tid];
		PerThreadLATSampleIterMap[tid][counter[tid]] = access;
		// cout << "Update LAT and LATSampleIterMap for sample: " <<  access->toString() << " " << cnt << endl;
	}
#if defined(DEBUG)
	if (samples.empty())
		cout << "sampler terminate" << endl;
#endif
	counter[tid]++;
    _m.unlock();
}
void pluss_access(unsigned addr)
{
	unsigned tmp = addr & CACHE_MASK;
	if (LAT.find(tmp) != LAT.end()) {
		unsigned reuse = cnt - LAT[tmp];
		pluss_histogram_update(reuse, 1);
	}
	LAT[tmp] = cnt;
	cnt++;
}

void pluss_parallel_access(unsigned addr)
{
    _m.lock();
	unsigned tmp = addr & CACHE_MASK;
	if (LAT.find(tmp) != LAT.end()) {
		unsigned reuse = cnt - LAT[tmp];
		pluss_histogram_update(reuse, 1);
	}
	LAT[tmp] = cnt;
	cnt++;
    _m.unlock();

}

void pluss_bypass(unsigned skip_cnt)
{
	cnt += skip_cnt;
}

void pluss_per_thread_bypass(thread::id tid, unsigned skip_cnt)
{
	counter[tid] += skip_cnt;
}

void pluss_print_histogram()
{
    cout << "Start to dump reuse time" << endl;
	uint64_t sum = 0;
	for (auto pair : RIHist) {
		sum += pair.second;
		
	}
	for (auto pair : RIHist) {
		cout << pair.first << "," << pair.second << "," << (double)pair.second/sum << endl;
	}
#if defined(DEBUG)
	cout << "Iterate " << cnt << " memory accesses" << endl;
	if (!LAT.empty()) {
		cout << LAT.size() << " unique accesses" << endl;
	}
	else if (!PerArrayLAT.empty()) {
		unsigned unique_cnt = 0;
		for (auto pair : PerArrayLAT) {
			unique_cnt += pair.second.size();
		}
		cout << unique_cnt << " unique accesses" << endl;
	}
#endif
}

void pluss_AET() {
	map<uint64_t, double> P;
	map<uint64_t, uint64_t> histogram;
	uint64_t total_num_RT = 0;
	uint64_t max_RT = 0;
	for (auto it = RIHist.begin(); it != RIHist.end(); ++it) {
		total_num_RT += it->second;
		histogram[it->first] = it->second;
		if (max_RT < it->first) {
			max_RT = it->first;
		}
	}
	uint64_t accumulate_num_RT = 0;
	for (auto it = histogram.rbegin(); it != histogram.rend(); ++it) {
		P[it->first] = (double)accumulate_num_RT / total_num_RT;
		accumulate_num_RT += it->second;
	}
	P[0] = 1;
	double sum_P = 0;
	uint64_t t = 0;
	uint64_t prev_t = 0;
	uint64_t cs = POLYBENCH_CACHE_SIZE_KB * 1024 / sizeof(double);
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
		MRC[c] = P[prev_t];
	}
	return;
}

void pluss_print_mrc()
{
	cout << "miss ratio" << endl;
	for (auto mit = MRC.begin(); mit != MRC.end(); mit++) {
		cout << mit->first << "," << mit->second << endl;
	}
#if 0	
    if (cnt > 0) {
        cout << "max iteration traversed" << endl;
        cout << cnt << endl;
    }
#endif
}
void pluss_terminate()
{
	if (!LATSampleIterMap.empty()) {
		for (auto entry : LATSampleIterMap) {
			free(entry.second);
		}
	}
}
#endif
