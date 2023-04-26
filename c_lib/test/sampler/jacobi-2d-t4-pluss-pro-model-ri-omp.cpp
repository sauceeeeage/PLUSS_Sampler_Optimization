#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <queue>
#include <numeric>
#include <time.h>
#include "pluss.h"
#include "pluss_utils.h"
using namespace std;
long max_iteration_count = 0;
uint64_t GetAddress_B1(int idx0,int idx1) {
    uint64_t addr_B1 = (idx0*1024)+(idx1*1);
    return addr_B1*DS/CLS;
}
uint64_t GetAddress_B0(int idx0,int idx1) {
    uint64_t addr_B0 = (idx0*1024)+(idx1*1);
    return addr_B0*DS/CLS;
}
uint64_t GetAddress_A4(int idx0,int idx1) {
    uint64_t addr_A4 = (idx0*1024)+(idx1*1);
    return addr_A4*DS/CLS;
}
uint64_t GetAddress_B4(int idx0,int idx1) {
    uint64_t addr_B4 = (idx0*1024)+(idx1*1);
    return addr_B4*DS/CLS;
}
uint64_t GetAddress_A0(int idx0,int idx1) {
    uint64_t addr_A0 = (idx0*1024)+(idx1*1);
    return addr_A0*DS/CLS;
}
uint64_t GetAddress_B2(int idx0,int idx1) {
    uint64_t addr_B2 = (idx0*1024)+(idx1*1);
    return addr_B2*DS/CLS;
}
uint64_t GetAddress_B5(int idx0,int idx1) {
    uint64_t addr_B5 = (idx0*1024)+(idx1*1);
    return addr_B5*DS/CLS;
}
uint64_t GetAddress_B3(int idx0,int idx1) {
    uint64_t addr_B3 = (idx0*1024)+(idx1*1);
    return addr_B3*DS/CLS;
}
uint64_t GetAddress_A3(int idx0,int idx1) {
    uint64_t addr_A3 = (idx0*1024)+(idx1*1);
    return addr_A3*DS/CLS;
}
uint64_t GetAddress_A1(int idx0,int idx1) {
    uint64_t addr_A1 = (idx0*1024)+(idx1*1);
    return addr_A1*DS/CLS;
}
uint64_t GetAddress_A5(int idx0,int idx1) {
    uint64_t addr_A5 = (idx0*1024)+(idx1*1);
    return addr_A5*DS/CLS;
}
uint64_t GetAddress_A2(int idx0,int idx1) {
    uint64_t addr_A2 = (idx0*1024)+(idx1*1);
    return addr_A2*DS/CLS;
}
//Generate sampler without sampling
void sampler() {
    //Declare components will be used in Parallel RI search
    array<Progress *, THREAD_NUM> progress = { nullptr };
    array<int, THREAD_NUM> idle_threads;
    vector<int> subscripts;
    ChunkDispatcher dispatcher;
    int tid_to_run = 0, start_tid=0, working_threads=THREAD_NUM;
    uint64_t addr = 0, loop_cnt = 0;
    array<long, THREAD_NUM+1> count = {0};
    srand(time(NULL));
    array<unordered_map<uint64_t, long>, THREAD_NUM> LAT_B;
    array<unordered_map<uint64_t, long>, THREAD_NUM> LAT_A;
#if defined(DEBUG)
    unordered_map<int, unordered_map<long, Iteration *>> LATIterMap;
#endif
#if defined(DEBUG)
    Iteration *access = nullptr;
#endif
    //Generate parallel code for (c0,1,\<,1024,(c0 + 1))
    //Threads are scheduled using static scheduling
    //Threads are interleaved using uniform interleaving
    loop_cnt += 1;
    dispatcher = ChunkDispatcher(CHUNK_SIZE,((1024-1)),1,1);
    for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
#if defined(DEBUG)
        cout << "Move " << tid_to_run << " to idle_threads" << endl;
#endif
        idle_threads[tid_to_run] = 1;
    }
    #pragma omp parallel for num_threads(THREAD_NUM) private(addr)
    for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
        while(true) {
            if (idle_threads[tid_to_run] == 1 && dispatcher.hasNextChunk(1)) {
                Chunk c = dispatcher.getNextStaticChunk(tid_to_run);
                vector<int> parallel_iteration_vector;
                parallel_iteration_vector.push_back(c.first);
                parallel_iteration_vector.push_back(1 );
                if (progress[tid_to_run]) {
                    progress[tid_to_run]->ref = "A0";
                    progress[tid_to_run]->iteration = parallel_iteration_vector;
                    progress[tid_to_run]->chunk = c;
                } else {
                    Progress *p = new Progress("A0", parallel_iteration_vector, c);
                    progress[tid_to_run] = p;
                }
                idle_threads[tid_to_run] = 0;
#if defined(DEBUG)
                cout << "Move " << tid_to_run << " to worker_threads" << endl;
#endif
                ///* end of progress assignment */
            } /* end of chunk availability check */
            //UNIFORM INTERLEAVING
            if (!progress[tid_to_run] || !progress[tid_to_run]->isInBound()) {
                idle_threads[tid_to_run] = 1;
                break;
            }
#if defined(DEBUG)
            cout << "[" << tid_to_run << "] Access " << progress[tid_to_run]->ref << " {";
            for (unsigned i = 0; i < progress[tid_to_run]->iteration.size(); i++) {
                cout << progress[tid_to_run]->iteration[i];
                if (i < progress[tid_to_run]->iteration.size()-1) { cout << ","; }
            }
            cout << "}" << endl;
#endif
            if (progress[tid_to_run]->ref == "A0") {
                addr = GetAddress_A0(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
                access = new Iteration("A0", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_A[tid_to_run].find(addr) != LAT_A[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_A[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_A[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_A[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                progress[tid_to_run]->increment("A1");
                continue;
            } /* end of check to A0 */
            if (progress[tid_to_run]->ref == "A1") {
                addr = GetAddress_A1(progress[tid_to_run]->iteration[0],(progress[tid_to_run]->iteration[1] + -1));
#if defined(DEBUG)
                access = new Iteration("A1", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_A[tid_to_run].find(addr) != LAT_A[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_A[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_A[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_A[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                progress[tid_to_run]->increment("A2");
                continue;
            } /* end of check to A1 */
            if (progress[tid_to_run]->ref == "A2") {
                addr = GetAddress_A2(progress[tid_to_run]->iteration[0],(progress[tid_to_run]->iteration[1] + 1));
#if defined(DEBUG)
                access = new Iteration("A2", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_A[tid_to_run].find(addr) != LAT_A[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_A[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_A[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_A[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                progress[tid_to_run]->increment("A3");
                continue;
            } /* end of check to A2 */
            if (progress[tid_to_run]->ref == "A3") {
                addr = GetAddress_A3((progress[tid_to_run]->iteration[0] + 1),progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
                access = new Iteration("A3", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_A[tid_to_run].find(addr) != LAT_A[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_A[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_A[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_A[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                progress[tid_to_run]->increment("A4");
                continue;
            } /* end of check to A3 */
            if (progress[tid_to_run]->ref == "A4") {
                addr = GetAddress_A4((progress[tid_to_run]->iteration[0] + -1),progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
                access = new Iteration("A4", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_A[tid_to_run].find(addr) != LAT_A[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_A[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_A[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_A[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                progress[tid_to_run]->increment("B0");
                continue;
            } /* end of check to A4 */
            if (progress[tid_to_run]->ref == "B0") {
                addr = GetAddress_B0(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
                access = new Iteration("B0", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_B[tid_to_run].find(addr) != LAT_B[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_B[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_B[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_B[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                //CASE 3
                if ((progress[tid_to_run]->iteration[1] + 1)<1024) {
                    progress[tid_to_run]->iteration[1] = (progress[tid_to_run]->iteration[1] + 1);
                    progress[tid_to_run]->increment("A0");
                    continue;
                } /* end of check to B0 */
                //CASE 1
                progress[tid_to_run]->iteration[0] = (progress[tid_to_run]->iteration[0] + 1);
                if (progress[tid_to_run]->isInBound()) {
                    progress[tid_to_run]->iteration.pop_back();
                    progress[tid_to_run]->iteration.emplace_back(1);
                    progress[tid_to_run]->increment("A0");
                    continue;
                } /* end of check to B0 */
            } /* end of thread interleaving loop */
            if (idle_threads[tid_to_run] == 0) {
#if defined(DEBUG)
                cout << "Move " << tid_to_run << " to idle threads" << endl;
#endif
                idle_threads[tid_to_run] = 1;
            } // end of if (idle_threads[tid_to_run] == 0)
            if (idle_threads[tid_to_run] == 1 && !dispatcher.hasNextChunk(1)) {
#if defined(DEBUG)
                cout << "No chunk available, sampler terminates" << endl;
#endif
                break;
            } /* end of break condition check */
        } /* end of while(true) */
    } /* end of omp #pragma */
    //reset both lists so they can be reused for later parallel loop
    idle_threads.fill(0);
    //Addresses in A with no data reuse will be marked as -1
    for (unsigned i = 0; i < LAT_A.size(); i++) {
    pluss_cri_noshare_histogram_update(i,-1,LAT_A[i].size());
    LAT_A[i].clear();
    }
    //Addresses in B with no data reuse will be marked as -1
    for (unsigned i = 0; i < LAT_B.size(); i++) {
    pluss_cri_noshare_histogram_update(i,-1,LAT_B[i].size());
    LAT_B[i].clear();
    }
    for (unsigned i = 0; i < progress.size(); i++) {
        if (progress[i]) {
            delete progress[i];
            progress[i] = nullptr;
        }
    }
    //Generate parallel code for (c0,1,\<,1024,(c0 + 1))
    //Threads are scheduled using static scheduling
    //Threads are interleaved using uniform interleaving
    loop_cnt += 1;
    dispatcher = ChunkDispatcher(CHUNK_SIZE,((1024-1)),1,1);
    for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
#if defined(DEBUG)
        cout << "Move " << tid_to_run << " to idle_threads" << endl;
#endif
        idle_threads[tid_to_run] = 1;
    }
    #pragma omp parallel for num_threads(THREAD_NUM) private(addr)
    for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
        while(true) {
            if (idle_threads[tid_to_run] == 1 && dispatcher.hasNextChunk(1)) {
                Chunk c = dispatcher.getNextStaticChunk(tid_to_run);
                vector<int> parallel_iteration_vector;
                parallel_iteration_vector.push_back(c.first);
                parallel_iteration_vector.push_back(1 );
                if (progress[tid_to_run]) {
                    progress[tid_to_run]->ref = "B1";
                    progress[tid_to_run]->iteration = parallel_iteration_vector;
                    progress[tid_to_run]->chunk = c;
                } else {
                    Progress *p = new Progress("B1", parallel_iteration_vector, c);
                    progress[tid_to_run] = p;
                }
                idle_threads[tid_to_run] = 0;
#if defined(DEBUG)
                cout << "Move " << tid_to_run << " to worker_threads" << endl;
#endif
                ///* end of progress assignment */
            } /* end of chunk availability check */
            //UNIFORM INTERLEAVING
            if (!progress[tid_to_run] || !progress[tid_to_run]->isInBound()) {
                idle_threads[tid_to_run] = 1;
                break;
            }
#if defined(DEBUG)
            cout << "[" << tid_to_run << "] Access " << progress[tid_to_run]->ref << " {";
            for (unsigned i = 0; i < progress[tid_to_run]->iteration.size(); i++) {
                cout << progress[tid_to_run]->iteration[i];
                if (i < progress[tid_to_run]->iteration.size()-1) { cout << ","; }
            }
            cout << "}" << endl;
#endif
            if (progress[tid_to_run]->ref == "B1") {
                addr = GetAddress_B1(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
                access = new Iteration("B1", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_B[tid_to_run].find(addr) != LAT_B[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_B[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_B[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_B[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                progress[tid_to_run]->increment("B2");
                continue;
            } /* end of check to B1 */
            if (progress[tid_to_run]->ref == "B2") {
                addr = GetAddress_B2(progress[tid_to_run]->iteration[0],(progress[tid_to_run]->iteration[1] + -1));
#if defined(DEBUG)
                access = new Iteration("B2", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_B[tid_to_run].find(addr) != LAT_B[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_B[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_B[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_B[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                progress[tid_to_run]->increment("B3");
                continue;
            } /* end of check to B2 */
            if (progress[tid_to_run]->ref == "B3") {
                addr = GetAddress_B3(progress[tid_to_run]->iteration[0],(progress[tid_to_run]->iteration[1] + 1));
#if defined(DEBUG)
                access = new Iteration("B3", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_B[tid_to_run].find(addr) != LAT_B[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_B[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_B[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_B[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                progress[tid_to_run]->increment("B4");
                continue;
            } /* end of check to B3 */
            if (progress[tid_to_run]->ref == "B4") {
                addr = GetAddress_B4((progress[tid_to_run]->iteration[0] + 1),progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
                access = new Iteration("B4", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_B[tid_to_run].find(addr) != LAT_B[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_B[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_B[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_B[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                progress[tid_to_run]->increment("B5");
                continue;
            } /* end of check to B4 */
            if (progress[tid_to_run]->ref == "B5") {
                addr = GetAddress_B5((progress[tid_to_run]->iteration[0] + -1),progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
                access = new Iteration("B5", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_B[tid_to_run].find(addr) != LAT_B[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_B[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_B[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_B[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                progress[tid_to_run]->increment("A5");
                continue;
            } /* end of check to B5 */
            if (progress[tid_to_run]->ref == "A5") {
                addr = GetAddress_A5(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
                access = new Iteration("A5", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
                if (LAT_A[tid_to_run].find(addr) != LAT_A[tid_to_run].end()) {
/* (c1,1,\<,1024,(c1 + 1)), 0*/
/* (c0,1,\<,1024,(c0 + 1)), 1*/
                    long reuse = count[tid_to_run] - LAT_A[tid_to_run][addr];
#if defined(DEBUG)
                    if (reuse >= 512) {
                        Iteration *src = LATIterMap[tid_to_run][LAT_A[tid_to_run][addr]];
                        cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                    }
#endif
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
                LAT_A[tid_to_run][addr] = count[tid_to_run];
#if defined(DEBUG)
                LATIterMap[tid_to_run][count[tid_to_run]] = access;
#endif
                count[tid_to_run] += 1;
                //CASE 3
                if ((progress[tid_to_run]->iteration[1] + 1)<1024) {
                    progress[tid_to_run]->iteration[1] = (progress[tid_to_run]->iteration[1] + 1);
                    progress[tid_to_run]->increment("B1");
                    continue;
                } /* end of check to A5 */
                //CASE 1
                progress[tid_to_run]->iteration[0] = (progress[tid_to_run]->iteration[0] + 1);
                if (progress[tid_to_run]->isInBound()) {
                    progress[tid_to_run]->iteration.pop_back();
                    progress[tid_to_run]->iteration.emplace_back(1);
                    progress[tid_to_run]->increment("B1");
                    continue;
                } /* end of check to A5 */
            } /* end of thread interleaving loop */
            if (idle_threads[tid_to_run] == 0) {
#if defined(DEBUG)
                cout << "Move " << tid_to_run << " to idle threads" << endl;
#endif
                idle_threads[tid_to_run] = 1;
            } // end of if (idle_threads[tid_to_run] == 0)
            if (idle_threads[tid_to_run] == 1 && !dispatcher.hasNextChunk(1)) {
#if defined(DEBUG)
                cout << "No chunk available, sampler terminates" << endl;
#endif
                break;
            } /* end of break condition check */
        } /* end of while(true) */
    } /* end of omp #pragma */
    //reset both lists so they can be reused for later parallel loop
    idle_threads.fill(0);
    //Addresses in B with no data reuse will be marked as -1
    for (unsigned i = 0; i < LAT_B.size(); i++) {
    pluss_cri_noshare_histogram_update(i,-1,LAT_B[i].size());
    LAT_B[i].clear();
    }
    //Addresses in A with no data reuse will be marked as -1
    for (unsigned i = 0; i < LAT_A.size(); i++) {
    pluss_cri_noshare_histogram_update(i,-1,LAT_A[i].size());
    LAT_A[i].clear();
    }
    for (unsigned i = 0; i < progress.size(); i++) {
        if (progress[i]) {
            delete progress[i];
            progress[i] = nullptr;
        }
    }
#if defined(DEBUG)
    for (int t = 0; t < THREAD_NUM; t++) {
        for (auto entry : LATIterMap[t]) { delete entry.second; }
    }
    LATIterMap.clear();
#endif
    max_iteration_count = accumulate(count.begin(), count.end(), 0);
}
int main() {
    pluss_timer_start();
    sampler();
    pluss_cri_distribute(THREAD_NUM);
    pluss_AET();
    pluss_timer_stop();
    pluss_timer_print();
    pluss_cri_noshare_print_histogram();
    pluss_cri_share_print_histogram();
    pluss_print_histogram();
    pluss_print_mrc();
    cout << "max iteration traversed" << endl;
    cout << max_iteration_count << endl;
    return 0;
}
