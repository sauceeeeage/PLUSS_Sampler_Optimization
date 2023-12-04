#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <queue>
#include <numeric>
#include <time.h>
#include "pluss.h"
#include "pluss_utils_v2.h"
using namespace std;
long max_iteration_count = 0;
uint64_t GetAddress_B0(int idx0,int idx1) {
    uint64_t addr_B0 = (idx0*128)+(idx1*1);
    return addr_B0*DS/CLS;
}
uint64_t GetAddress_A0(int idx0,int idx1) {
    uint64_t addr_A0 = (idx0*128)+(idx1*1);
    return addr_A0*DS/CLS;
}
uint64_t GetAddress_C0(int idx0,int idx1) {
    uint64_t addr_C0 = (idx0*128)+(idx1*1);
    return addr_C0*DS/CLS;
}
uint64_t GetAddress_C3(int idx0,int idx1) {
    uint64_t addr_C3 = (idx0*128)+(idx1*1);
    return addr_C3*DS/CLS;
}
uint64_t GetAddress_C1(int idx0,int idx1) {
    uint64_t addr_C1 = (idx0*128)+(idx1*1);
    return addr_C1*DS/CLS;
}
uint64_t GetAddress_C2(int idx0,int idx1) {
    uint64_t addr_C2 = (idx0*128)+(idx1*1);
    return addr_C2*DS/CLS;
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
    srand(time(NULL));
#if defined(DEBUG)
    unordered_map<int, unordered_map<long, Iteration *>> LATIterMap;
#endif
#if defined(DEBUG)
    Iteration *access = nullptr;
#endif
    //Generate parallel code for (c0,0,\<,128,(c0 + 1))
    //Threads are scheduled using static scheduling
    //Threads are interleaved using uniform interleaving
    loop_cnt += 1;
    dispatcher = ChunkDispatcher(CHUNK_SIZE,((128-0)),0,1);
    for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
#if defined(DEBUG)
        cout << "Move " << tid_to_run << " to idle_threads" << endl;
#endif
        idle_threads[tid_to_run] = 1;
    }
    #pragma omp parallel for num_threads(1) private(addr)
    for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
        long count = 0;
        unordered_map<uint64_t, long> LAT_B;
        unordered_map<uint64_t, long> LAT_A;
        unordered_map<uint64_t, long> LAT_C;
        while(true) {
            if (idle_threads[tid_to_run] == 1 && dispatcher.hasNextChunk(1)) {
                Chunk c = dispatcher.getNextStaticChunk(tid_to_run);
                vector<int> parallel_iteration_vector;
                parallel_iteration_vector.push_back(c.first);
                parallel_iteration_vector.push_back(0 );
                if (progress[tid_to_run]) {
                    // progress[tid_to_run]->ref = "C0";
                    progress[tid_to_run]->iteration = parallel_iteration_vector;
                    progress[tid_to_run]->chunk = c;
                } else {
                    Progress *p = new Progress("C0", parallel_iteration_vector, c);
                    progress[tid_to_run] = p;
                }
                idle_threads[tid_to_run] = 0;
#if defined(DEBUG)
                cout << "Move " << tid_to_run << " to worker_threads" << endl;
#endif
                ///* end of progress assignment */
            } /* end of chunk availability check */
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
LOOP_J:
            addr = GetAddress_C0(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
            access = new Iteration("C0", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
            if (LAT_C.find(addr) != LAT_C.end()) {
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                long reuse = count- LAT_C[addr];
#if defined(DEBUG)
                if (reuse >= 512) {
                    Iteration *src = LATIterMap[tid_to_run][LAT_C[addr]];
                    cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                }
#endif
                pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
            }
            LAT_C[addr] = count;
#if defined(DEBUG)
            LATIterMap[tid_to_run][count] = access;
#endif
            count += 1;
            /* start check to C1 */
            addr = GetAddress_C1(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
            access = new Iteration("C1", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
#endif
            if (LAT_C.find(addr) != LAT_C.end()) {
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                long reuse = count - LAT_C[addr];
#if defined(DEBUG)
                if (reuse >= 512) {
                    Iteration *src = LATIterMap[tid_to_run][LAT_C[addr]];
                    cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                }
#endif
                pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
            }
            LAT_C[addr] = count;
#if defined(DEBUG)
            LATIterMap[tid_to_run][count] = access;
#endif
            count += 1;
            //CASE 2
            progress[tid_to_run]->iteration.emplace_back(0);
            /* end of check to C0 and C1 */
LOOP_K:
			/* check to A0 */
            addr = GetAddress_A0(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[2]);
#if defined(DEBUG)
            access = new Iteration("A0", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1],progress[tid_to_run]->iteration[2]});
#endif
            if (LAT_A.find(addr) != LAT_A.end()) {
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                long reuse = count - LAT_A[addr];
#if defined(DEBUG)
                if (reuse >= 512) {
                    Iteration *src = LATIterMap[tid_to_run][LAT_A[addr]];
                    cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                }
#endif
                pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
            }
            LAT_A[addr] = count;
#if defined(DEBUG)
            LATIterMap[tid_to_run][count] = access;
#endif
            count += 1;
            /* check to B0 */
            addr = GetAddress_B0(progress[tid_to_run]->iteration[2],progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
            access = new Iteration("B0", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1],progress[tid_to_run]->iteration[2]});
#endif
            if (LAT_B.find(addr) != LAT_B.end()) {
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                long reuse = count - LAT_B[addr];
#if defined(DEBUG)
                if (reuse >= 512) {
                    Iteration *src = LATIterMap[tid_to_run][LAT_B[addr]];
                    cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                }
#endif
/* Compare c2*/
/* With c2*/
/* With c1*/
/* Compare c1*/
/* With c2*/
/* With c1*/
                //B[c2][c1] is carried by (c1,0,\<,128,(c1 + 1))
                if (distance_to(reuse,0) > distance_to(reuse,(((1)*((128-0)/1)+1)*((128-0)/1)+1))) {
                    pluss_cri_share_histogram_update(tid_to_run,THREAD_NUM-1,reuse,1);
                } else {
                    pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
                }
            }
            LAT_B[addr] = count;
#if defined(DEBUG)
            LATIterMap[tid_to_run][count] = access;
#endif
            count += 1;
            /* check to C2 */
            addr = GetAddress_C2(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
            access = new Iteration("C2", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1],progress[tid_to_run]->iteration[2]});
#endif
            if (LAT_C.find(addr) != LAT_C.end()) {
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                long reuse = count - LAT_C[addr];
#if defined(DEBUG)
                if (reuse >= 512) {
                    Iteration *src = LATIterMap[tid_to_run][LAT_C[addr]];
                    cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                }
#endif
                pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
            }
            LAT_C[addr] = count;
#if defined(DEBUG)
            LATIterMap[tid_to_run][count] = access;
#endif
            count += 1;
            /* start check to C3 */
            addr = GetAddress_C3(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
#if defined(DEBUG)
            access = new Iteration("C3", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1],progress[tid_to_run]->iteration[2]});
#endif
            if (LAT_C.find(addr) != LAT_C.end()) {
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                long reuse = count - LAT_C[addr];
#if defined(DEBUG)
                if (reuse >= 512) {
                    Iteration *src = LATIterMap[tid_to_run][LAT_C[addr]];
                    cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
                }
#endif
                pluss_cri_noshare_histogram_update(tid_to_run,reuse,1);
            }
            LAT_C[addr] = count;
#if defined(DEBUG)
            LATIterMap[tid_to_run][count] = access;
#endif
            count += 1;
            //CASE 3
            if ((progress[tid_to_run]->iteration[2] + 1)<128) {
                progress[tid_to_run]->iteration[2] = (progress[tid_to_run]->iteration[2] + 1);
                goto LOOP_K;
            } /* end of check to C3 */
            //CASE 1
            if ((progress[tid_to_run]->iteration[1] + 1)<128) {
                progress[tid_to_run]->iteration[1] = (progress[tid_to_run]->iteration[1] + 1);
                progress[tid_to_run]->iteration.pop_back();
                goto LOOP_J;
            } /* end of check to C3 */
            //CASE 1
            progress[tid_to_run]->iteration[0] = (progress[tid_to_run]->iteration[0] + 1);
            if (progress[tid_to_run]->isInBound()) {
                progress[tid_to_run]->iteration.pop_back();
                progress[tid_to_run]->iteration.pop_back();
                progress[tid_to_run]->iteration.emplace_back(0);
                progress[tid_to_run]->increment("C0");
                continue;
            } /* end of check to C3 */
            /* end of check to A0 to C3 */
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
                #pragma omp critical
                {
                max_iteration_count += count;
                }
                //Addresses in A, B and C with no data reuse will be marked as -1
                pluss_cri_noshare_histogram_update(tid_to_run,-1,LAT_A.size());
                pluss_cri_noshare_histogram_update(tid_to_run,-1,LAT_B.size());
                pluss_cri_noshare_histogram_update(tid_to_run,-1,LAT_C.size());
                LAT_A.clear();
                LAT_B.clear();
                LAT_C.clear();
                break;
            } /* end of break condition check */
        } /* end of while(true) */
        // max_iteration_count += count;
    } /* end of omp #pragma */
    //reset both lists so they can be reused for later parallel loop
    idle_threads.fill(0);
    
#if 0
    //Addresses in C with no data reuse will be marked as -1
    for (unsigned i = 0; i < LAT_C.size(); i++) {
    pluss_cri_noshare_histogram_update(i,-1,LAT_C[i].size());
    LAT_C[i].clear();
    }
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
#endif
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
}
int main(int argc, char** argv) {
//    pluss_timer_start();
//    sampler();
//    // pluss_cri_distribute(THREAD_NUM);
//    pluss_AET();
//    pluss_timer_stop();
//    pluss_timer_print();
//    // pluss_cri_noshare_print_histogram();
//    // pluss_cri_share_print_histogram();
//    // pluss_print_histogram();
//    // pluss_print_mrc();
//    // cout << "max iteration traversed" << endl;
//    // cout << max_iteration_count << endl;
//    return 0;
    string method = argv[1];
    if (method == "acc") {
        pluss_timer_start();
        sampler();
        pluss_cri_distribute(THREAD_NUM);
        pluss_timer_stop();
        cout << "OPENMP OPT C++: ";
        pluss_timer_print();
        pluss_cri_noshare_print_histogram();
        pluss_cri_share_print_histogram();
        pluss_print_histogram();
        cout << "max iteration traversed" << endl;
        cout << max_iteration_count << endl;
        printf("\n");
    } else if (method == "speed") {
        cout << "OPENMP OPT C++:\n";
        for (int i = 0; i < 10; i++) {
            pluss_timer_start();
            sampler();
            pluss_cri_distribute(THREAD_NUM);
            pluss_timer_stop();
            pluss_timer_print();
        }
        printf("\n");
    }

    return 0;
}
