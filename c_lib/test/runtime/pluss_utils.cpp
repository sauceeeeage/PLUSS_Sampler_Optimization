#include "pluss_utils.h"

using namespace std;
Histogram _RIHist;
unordered_map<long, atomic<double>> _ParRIHist;
map<uint64_t, double> _MRC;


// CRI
array<Histogram, THREAD_NUM> _NoSharePRI;
array<unordered_map<int, Histogram>, THREAD_NUM> _SharePRI;


shared_mutex _m;