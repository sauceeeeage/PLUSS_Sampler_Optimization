#include "pluss_utils_v2.h"

using namespace std;
using namespace phmap;
Histogram _RIHist;
map<uint64_t, double> _MRC;


// CRI
array<Histogram, THREAD_NUM> _NoSharePRI;
array<flat_hash_map<int, Histogram>, THREAD_NUM> _SharePRI;