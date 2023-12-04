#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <queue>
#include <thread>
#include <sched.h>
#include <pthread.h>
#include <numeric>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h> // cdf of gaussian distribution
#include <time.h>
#include "pluss.h"
#include "pluss_utils.h"
using namespace std;
unordered_map<string, long> iteration_traversed_map;
uint64_t GetAddress_C3(int idx0,int idx1) {
    uint64_t addr_C3 = (idx0*128)+(idx1*1);
    return addr_C3*DS/CLS;
}
uint64_t GetAddress_C2(int idx0,int idx1) {
    uint64_t addr_C2 = (idx0*128)+(idx1*1);
    return addr_C2*DS/CLS;
}
uint64_t GetAddress_A0(int idx0,int idx1) {
    uint64_t addr_A0 = (idx0*128)+(idx1*1);
    return addr_A0*DS/CLS;
}
uint64_t GetAddress_C0(int idx0,int idx1) {
    uint64_t addr_C0 = (idx0*128)+(idx1*1);
    return addr_C0*DS/CLS;
}
uint64_t GetAddress_B0(int idx0,int idx1) {
    uint64_t addr_B0 = (idx0*128)+(idx1*1);
    return addr_B0*DS/CLS;
}
uint64_t GetAddress_C1(int idx0,int idx1) {
    uint64_t addr_C1 = (idx0*128)+(idx1*1);
    return addr_C1*DS/CLS;
}
void simulate_negative_binomial(int thread_cnt, long n, unordered_map<long, double> &dist)
{
#if defined(DEBUG)
    cout << "Distribute thread-local Reuse: " << n << endl;
#endif
    double p = 1.0 / thread_cnt;
    if (n >= (4000. * (thread_cnt-1)) / thread_cnt) {
        int i = (int)log2(n);
        long bin = (long)(pow(2.0, i));
        dist[THREAD_NUM*bin] = 1.0;
        return;
    }
    uint64_t k = 0;
    double nbd_prob = 0.0, prob_sum = 0.0;
    while (true) {
        nbd_prob = gsl_ran_negative_binomial_pdf(k, p, n);
        prob_sum += nbd_prob;
        dist[k+n] = nbd_prob;
        if (prob_sum > 0.999)
            break;
        k += 1;
    }
} // end of void simulate_negative_binomial()
void no_share_distribute(unordered_map<long, double>& histogram_to_distribute, unordered_map<long, double>& target_histogram, int thread_cnt=THREAD_NUM)
{
    unordered_map<long, double> dist;
    for (auto entry : histogram_to_distribute) {
        if (entry.first < 0) {
            pluss_parallel_histogram_update(target_histogram,entry.first,entry.second);
            continue;
        };
        if (thread_cnt > 1) {
            simulate_negative_binomial(thread_cnt, entry.first, dist);
            for (auto dist_entry : dist) {
                long ri_to_distribute = dist_entry.first;
                pluss_parallel_histogram_update(target_histogram,ri_to_distribute,entry.second * dist_entry.second);
            }
            dist.clear();
        } else {
            pluss_parallel_histogram_update(target_histogram,entry.first,entry.second);
        } // end of if(thread_cnt > 1)
    }
} // end of void no_share_distribute()
void share_distribute(unordered_map<int, unordered_map<long, double>>&histogram_to_distribute, unordered_map<long, double> &target_histogram, int thread_cnt=THREAD_NUM)
{
    unordered_map<int, double> prob;
    unordered_map<long, double> dist;
    for (auto share_entry : histogram_to_distribute) {
        int i = 1;
        double prob_sum = 0.0, n = (double)share_entry.first;
        for (auto reuse_entry : share_entry.second) {
            if (thread_cnt > 1) {
                simulate_negative_binomial(1.0/THREAD_NUM, reuse_entry.first, dist);
                for (auto dist_entry: dist) {
                    long ri_to_distribute = dist_entry.first;
                    double cnt_to_distribute = reuse_entry.second * dist_entry.second;
#if defined(DEBUG)
                    printf("new ri: 	%lu, %lu*%f=%f\n", ri_to_distribute, reuse_entry.second, dist_entry.second, cnt_to_distribute);
#endif
                    prob_sum = 0.0;
                    i = 1;
                    while (true) {
                        if (pow(2.0, (double)i) > ri_to_distribute) { break; }
                        prob[i] = pow(1-(pow(2.0, (double)i-1) / ri_to_distribute), n-1) - pow(1-(pow(2.0, (double)i) / ri_to_distribute), n-1);
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
                        pluss_parallel_histogram_update(target_histogram,new_ri,bin.second*cnt_to_distribute);
                    }
                    prob.clear();
                } // end of iterating dist
                dist.clear();
            } else {
                pluss_parallel_histogram_update(target_histogram,reuse_entry.first,reuse_entry.second);
            }
        } // end of iterating all reuse entries
    } // end of iterating all type of reuses
} // end of void share_distribute()
//Generate sampler with sampling
//Use random start sampling with ratio 10%
//C[c0][c1]
void sampler_C3(unordered_map<long, double> &histogram) {
    //Declare components will be used in Parallel RI search
    array<Progress *, THREAD_NUM> progress = { nullptr };
    vector<int> idle_threads, worker_threads, subscripts;
    ChunkDispatcher dispatcher;
    int tid_to_run = 0, start_tid = 0, working_threads = THREAD_NUM;
    uint64_t addr = 0;
    //global counter
    array<long, THREAD_NUM> count = {0};
    unordered_map<int, unordered_map<uint64_t, long>> LAT;
    unordered_map<int, unordered_map<long, Iteration *>> LATSampleIterMap;
    unordered_map<long, double> no_share_histogram;
    unordered_map<int, unordered_map<long, double>> share_histogram;
    unsigned s=0, traversing_samples=0;
    bool start_reuse_search=false;
    Iteration start_iteration;
    priority_queue<Iteration, vector<Iteration>, IterationComp> samples;
    unordered_set<string> sample_names, samples_meet;
    hash<std::thread::id> hasher;
    static thread_local mt19937 generator(time(NULL)+hasher(this_thread::get_id()));
    //generate all samples and sort them in order
    int num_samples = 2098;
    while (s < num_samples) {
        string sample_label = "C3";
        int sample_c0 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c0: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c0 << endl;
#endif
        if (((128-0)/1-((128-0)%1==0))< 0) { continue; }
        int sample_c1 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0); //this just give whatever rand() is
#if defined(DEBUG)
        cout << "sample_c1: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c1 << endl;
#endif
        if (((128-0)/1-((128-0)%1==0))< 0) { continue; }
        int sample_c2 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c2: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c2 << endl;
#endif
        sample_label += (to_string(sample_c0) + "-");
        sample_label += (to_string(sample_c1) + "-");
        sample_label += (to_string(sample_c2) + "-");
        //avoid duplicate sampling
        if (sample_names.find(sample_label)!=sample_names.end()) { continue; }
        sample_names.insert(sample_label);
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
        Iteration sample("C3", {sample_c0,sample_c1,sample_c2},0,1,true,0);
        samples.push(sample);
        s++;
    } // the end of while(s < NUM_SAMPLES)
    while (!samples.empty()) {
START_SAMPLE_C3:
        start_reuse_search = false;
        start_iteration = samples.top();
        traversing_samples++;
        samples.pop();
        if (!idle_threads.empty()) { idle_threads.clear(); };
        if (!worker_threads.empty()) { worker_threads.clear(); };
        if (!LAT.empty()) { ;
            for (unsigned i = 0; i < LAT.size(); i++) {
                pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
                LAT.clear();;
            }
        }
        if (!LATSampleIterMap.empty()) {
            for (auto LATSampleIterMapEntry : LATSampleIterMap) {
                for (auto entry : LATSampleIterMapEntry.second) {
                    delete entry.second;
                }
                LATSampleIterMapEntry.second.clear();
            }
            LATSampleIterMap.clear();
        }
#if defined(DEBUG)
        cout << "Start tracking sample " << start_iteration.toString() << endl;
#endif
        int sample_c0_start=start_iteration.ivs[0];
        int c0_start=0;
        int sample_c1_start=start_iteration.ivs[1];
        int c1_start=0;
        int sample_c2_start=start_iteration.ivs[2];
        int c2_start=0;
        //Generate parallel code for (c0,0,\<,128,(c0 + 1))
        //Threads are scheduled using static scheduling
        //Threads are interleaved using uniform interleaving
        dispatcher = ChunkDispatcher(CHUNK_SIZE,((128-0)),0,1);
        dispatcher.reset();
        idle_threads.clear();
        worker_threads.clear();
        for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
#if defined(DEBUG)
            cout << "Move " << tid_to_run << " to idle_threads" << endl;
#endif
            idle_threads.emplace_back(tid_to_run);
        }
        bool start_parallel_sample = true;
        if (start_parallel_sample) {
            dispatcher.setStartPoint(start_iteration.ivs[0]);
            start_tid = dispatcher.getStaticTid(start_iteration.ivs[0]);
            for (auto idle_threads_iterator = idle_threads.begin();idle_threads_iterator != idle_threads.end();) {
                int idle_tid = *idle_threads_iterator;
                if (dispatcher.hasNextStaticChunk(idle_tid)) {
                    Chunk start_chunk = dispatcher.getStaticStartChunk(start_iteration.ivs[0], idle_tid);
                    if (start_chunk.first > start_chunk.second) { 
                        //though the idle_tid has an available chunk, but at the given sampling point, there is no iteration to execute for this idle_tid
#if defined(DEBUG)
                        cout << "[" << start_chunk.first << ">" << start_chunk.second << "] No available iteration for " << idle_tid <<  endl;
#endif
                        idle_threads_iterator++;
                        continue;
                    }
                    if (progress[idle_tid]) {
                        progress[idle_tid]->ref = "C3";
                        progress[idle_tid]->iteration = {start_chunk.first, start_iteration.ivs[1], start_iteration.ivs[2], };
                        progress[idle_tid]->chunk = start_chunk;
                    } else {
                        Progress *p = new Progress("C3", {start_chunk.first, start_iteration.ivs[1], start_iteration.ivs[2], }, start_chunk);
                        progress[idle_tid] = p;
                    }
#if defined(DEBUG)
                    cout << "[" << idle_tid << "] starts from " << progress[idle_tid]->toString() << endl;
#endif
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
#if defined(DEBUG)
                    cout << "Move " << idle_tid << " to worker_threads" << endl;
#endif
                    worker_threads.push_back(idle_tid);
                    continue;
                }
                idle_threads_iterator++;
            }
            start_parallel_sample = false;
            //assign the start sample to each thread
            if (!worker_threads.empty()) {
                goto INTERLEAVING_LOOP_0;
            } else {
                goto END_SAMPLE;
            }
        } // end of if (start_parallel_sample)
        while(true) {
            if (!idle_threads.empty() && dispatcher.hasNextChunk(1)) {
                auto idle_threads_iterator = idle_threads.begin();
                while(idle_threads_iterator != idle_threads.end() && dispatcher.hasNextChunk(1)) {
                    tid_to_run = *idle_threads_iterator;
                    if (!dispatcher.hasNextStaticChunk(tid_to_run)) {
                        idle_threads_iterator++;
                        continue;
                    }
                    Chunk c = dispatcher.getNextStaticChunk(tid_to_run);
                    vector<int> parallel_iteration_vector;
                    parallel_iteration_vector.push_back(c.first);
                    parallel_iteration_vector.push_back(0);
                    if (progress[tid_to_run]) {
                        progress[tid_to_run]->ref = "C0";
                        progress[tid_to_run]->iteration = parallel_iteration_vector;
                        progress[tid_to_run]->chunk = c;
                    } else {
                        Progress *p = new Progress("C0", parallel_iteration_vector, c);
                        progress[tid_to_run] = p;
                    }
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
                    worker_threads.emplace_back(tid_to_run);
#if defined(DEBUG)
                    cout << "Move " << tid_to_run << " to worker_threads" << endl;
#endif
                } /* end of progress assignment */
            } /* end of chunk availability check */
INTERLEAVING_LOOP_0:
            //UNIFORM INTERLEAVING
            working_threads = (int)worker_threads.size();
            sort(worker_threads.begin(), worker_threads.end());
            auto worker_thread_iterator = worker_threads.begin();
            while (worker_thread_iterator != worker_threads.end()) {
                tid_to_run = *worker_thread_iterator;
                if (!progress[tid_to_run]->isInBound()) {
                    worker_thread_iterator++;
#if defined(DEBUG)
                    cout << "[" << tid_to_run << "] Out of bound" << endl;;
#endif
                    continue;
                }
#if defined(DEBUG)
                if (start_reuse_search) {
                    cout << "[" << tid_to_run << "] Access " << progress[tid_to_run]->ref << " {";
                    for (unsigned i = 0; i < progress[tid_to_run]->iteration.size(); i++) {
                        cout << progress[tid_to_run]->iteration[i];
                        if (i < progress[tid_to_run]->iteration.size()-1) { cout << ","; }
                    }
                    cout << "} " << start_reuse_search << endl;
                }
#endif
                if (progress[tid_to_run]->ref == "C0") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C0(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C3;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C1");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C0 */
                if (progress[tid_to_run]->ref == "C1") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C1(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C3;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    //CASE 2
                    progress[tid_to_run]->iteration.emplace_back(0);
                    progress[tid_to_run]->increment("A0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C1 */
                if (progress[tid_to_run]->ref == "A0") {
                    if (start_reuse_search) {
                        //A[c0][c2] will not access C
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("B0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to A0 */
                if (progress[tid_to_run]->ref == "B0") {
                    if (start_reuse_search) {
                        //B[c2][c1] will not access C
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C2");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to B0 */
                if (progress[tid_to_run]->ref == "C2") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C2(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C3;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C3");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C2 */
                if (progress[tid_to_run]->ref == "C3") {
                    if (!start_reuse_search) { start_reuse_search= (start_tid == tid_to_run); }
                    if (start_reuse_search) {
                        Iteration *access = new Iteration("C3", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1],progress[tid_to_run]->iteration[2]});
                        string array = "C";
                        addr = GetAddress_C3(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
#if defined(DEBUG)
                        cout << access->toString() << " @ " << addr << endl;
#endif
                        if (access->compare(start_iteration) == 0) {
#if defined(DEBUG)
                            cout << "Meet the start sample " << access->toString() << endl;
#endif
                            isSample = true;
                        } else if ((!samples.empty() && access->compare(samples.top()) == 0) || (sample_names.find(access->toAddrString()) != sample_names.end())) {
#if defined(DEBUG)
                            cout << "Meet a new sample " << access->toString() << " while searching reuses" << endl;
#endif
                            traversing_samples++;
                            if (!samples.empty() && access->compare(samples.top()) == 0) {
                                samples.pop();
                            }
                            isSample = true;
                            samples_meet.insert(access->toAddrString());
                        }
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
#if defined(DEBUG)
                            cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
#endif
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C3;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        if (isSample) {
                            LAT[tid_to_run][addr] = count[tid_to_run];
                            LATSampleIterMap[tid_to_run][count[tid_to_run]] = access;
                        } else { delete access; }
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    //CASE 3
                    if ((progress[tid_to_run]->iteration[2] + 1)<128) {
                        progress[tid_to_run]->iteration[2] = (progress[tid_to_run]->iteration[2] + 1);
                        progress[tid_to_run]->increment("A0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    if ((progress[tid_to_run]->iteration[1] + 1)<128) {
                        progress[tid_to_run]->iteration[1] = (progress[tid_to_run]->iteration[1] + 1);
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    progress[tid_to_run]->iteration[0] = (progress[tid_to_run]->iteration[0] + 1);
                    if (progress[tid_to_run]->isInBound()) {
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.emplace_back(0);
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    if (find(idle_threads.begin(), idle_threads.end(), tid_to_run) == idle_threads.end()) {
#if defined(DEBUG)
                        cout << "Move " << tid_to_run << " to idle threads" << endl;
#endif
                        idle_threads.emplace_back(tid_to_run);
                        worker_thread_iterator = worker_threads.erase(find(worker_threads.begin(), worker_threads.end(), tid_to_run));
                    }
                }
            } /* end of thread interleaving loop */
            if (idle_threads.size() == THREAD_NUM && !dispatcher.hasNextChunk(1)) {
                break;
            } /* end of break condition check */
        } /* end of while(true) */
        //reset both lists so they can be reused for later parallel loop
        idle_threads.clear();
        worker_threads.clear();
        for (unsigned i = 0; i < progress.size(); i++) {
            if (progress[i]) {
                delete progress[i];
                progress[i] = nullptr;
            }
        } // end of progress traversal
        goto END_SAMPLE;
    } // end of while (!samples.empty())
END_SAMPLE:
    if (!idle_threads.empty()) { idle_threads.clear(); };
    if (!worker_threads.empty()) { worker_threads.clear(); };
    if (!LAT.empty()) { ;
        for (unsigned i = 0; i < LAT.size(); i++) {
            pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
            LAT.clear();;
        }
    }
    for (unsigned i = 0; i < progress.size(); i++) {
        if (progress[i]) {
            delete progress[i];
            progress[i] = nullptr;
        }
    }
    if (!LATSampleIterMap.empty()) {
        for (auto LATSampleIterMapEntry : LATSampleIterMap) {
            for (auto entry : LATSampleIterMapEntry.second) {
                delete entry.second;
            }
            LATSampleIterMapEntry.second.clear();
        }
        LATSampleIterMap.clear();
    }
    no_share_distribute(no_share_histogram,histogram);
    share_distribute(share_histogram,histogram);
    no_share_histogram.clear();
    share_histogram.clear();
    iteration_traversed_map["C3"] = accumulate(count.begin(), count.end(), 0);
    return;
} // end of the sampler function
//C[c0][c1]
void sampler_C2(unordered_map<long, double> &histogram) {
    //Declare components will be used in Parallel RI search
    array<Progress *, THREAD_NUM> progress = { nullptr };
    vector<int> idle_threads, worker_threads, subscripts;
    ChunkDispatcher dispatcher;
    int tid_to_run = 0, start_tid = 0, working_threads = THREAD_NUM;
    uint64_t addr = 0;
    //global counter
    array<long, THREAD_NUM> count = {0};
    unordered_map<int, unordered_map<uint64_t, long>> LAT;
    unordered_map<int, unordered_map<long, Iteration *>> LATSampleIterMap;
    unordered_map<long, double> no_share_histogram;
    unordered_map<int, unordered_map<long, double>> share_histogram;
    unsigned s=0, traversing_samples=0;
    bool start_reuse_search=false;
    Iteration start_iteration;
    priority_queue<Iteration, vector<Iteration>, IterationComp> samples;
    unordered_set<string> sample_names, samples_meet;
    hash<std::thread::id> hasher;
    static thread_local mt19937 generator(time(NULL)+hasher(this_thread::get_id()));
    //generate all samples and sort them in order
    int num_samples = 2098;
    while (s < num_samples) {
        string sample_label = "C2";
        int sample_c0 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c0: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c0 << endl;
#endif
        if (((128-0)/1-((128-0)%1==0))< 0) { continue; }
        int sample_c1 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c1: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c1 << endl;
#endif
        if (((128-0)/1-((128-0)%1==0))< 0) { continue; }
        int sample_c2 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c2: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c2 << endl;
#endif
        sample_label += (to_string(sample_c0) + "-");
        sample_label += (to_string(sample_c1) + "-");
        sample_label += (to_string(sample_c2) + "-");
        //avoid duplicate sampling
        if (sample_names.find(sample_label)!=sample_names.end()) { continue; }
        sample_names.insert(sample_label);
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
        Iteration sample("C2", {sample_c0,sample_c1,sample_c2},0,1,true,0);
        samples.push(sample);
        s++;
    } // the end of while(s < NUM_SAMPLES)
    while (!samples.empty()) {
START_SAMPLE_C2:
        start_reuse_search = false;
        start_iteration = samples.top();
        traversing_samples++;
        samples.pop();
        if (!idle_threads.empty()) { idle_threads.clear(); };
        if (!worker_threads.empty()) { worker_threads.clear(); };
        if (!LAT.empty()) { ;
            for (unsigned i = 0; i < LAT.size(); i++) {
                pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
                LAT.clear();;
            }
        }
        if (!LATSampleIterMap.empty()) {
            for (auto LATSampleIterMapEntry : LATSampleIterMap) {
                for (auto entry : LATSampleIterMapEntry.second) {
                    delete entry.second;
                }
                LATSampleIterMapEntry.second.clear();
            }
            LATSampleIterMap.clear();
        }
#if defined(DEBUG)
        cout << "Start tracking sample " << start_iteration.toString() << endl;
#endif
        int sample_c0_start=start_iteration.ivs[0];
        int c0_start=0;
        int sample_c1_start=start_iteration.ivs[1];
        int c1_start=0;
        int sample_c2_start=start_iteration.ivs[2];
        int c2_start=0;
        //Generate parallel code for (c0,0,\<,128,(c0 + 1))
        //Threads are scheduled using static scheduling
        //Threads are interleaved using uniform interleaving
        dispatcher = ChunkDispatcher(CHUNK_SIZE,((128-0)),0,1);
        dispatcher.reset();
        idle_threads.clear();
        worker_threads.clear();
        for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
#if defined(DEBUG)
            cout << "Move " << tid_to_run << " to idle_threads" << endl;
#endif
            idle_threads.emplace_back(tid_to_run);
        }
        bool start_parallel_sample = true;
        if (start_parallel_sample) {
            dispatcher.setStartPoint(start_iteration.ivs[0]);
            start_tid = dispatcher.getStaticTid(start_iteration.ivs[0]);
            for (auto idle_threads_iterator = idle_threads.begin();idle_threads_iterator != idle_threads.end();) {
                int idle_tid = *idle_threads_iterator;
                if (dispatcher.hasNextStaticChunk(idle_tid)) {
                    Chunk start_chunk = dispatcher.getStaticStartChunk(start_iteration.ivs[0], idle_tid);
                    if (start_chunk.first > start_chunk.second) { 
                        //though the idle_tid has an available chunk, but at the given sampling point, there is no iteration to execute for this idle_tid
#if defined(DEBUG)
                        cout << "[" << start_chunk.first << ">" << start_chunk.second << "] No available iteration for " << idle_tid <<  endl;
#endif
                        idle_threads_iterator++;
                        continue;
                    }
                    if (progress[idle_tid]) {
                        progress[idle_tid]->ref = "C2";
                        progress[idle_tid]->iteration = {start_chunk.first, start_iteration.ivs[1], start_iteration.ivs[2], };
                        progress[idle_tid]->chunk = start_chunk;
                    } else {
                        Progress *p = new Progress("C2", {start_chunk.first, start_iteration.ivs[1], start_iteration.ivs[2], }, start_chunk);
                        progress[idle_tid] = p;
                    }
#if defined(DEBUG)
                    cout << "[" << idle_tid << "] starts from " << progress[idle_tid]->toString() << endl;
#endif
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
#if defined(DEBUG)
                    cout << "Move " << idle_tid << " to worker_threads" << endl;
#endif
                    worker_threads.push_back(idle_tid);
                    continue;
                }
                idle_threads_iterator++;
            }
            start_parallel_sample = false;
            //assign the start sample to each thread
            if (!worker_threads.empty()) {
                goto INTERLEAVING_LOOP_0;
            } else {
                goto END_SAMPLE;
            }
        } // end of if (start_parallel_sample)
        while(true) {
            if (!idle_threads.empty() && dispatcher.hasNextChunk(1)) {
                auto idle_threads_iterator = idle_threads.begin();
                while(idle_threads_iterator != idle_threads.end() && dispatcher.hasNextChunk(1)) {
                    tid_to_run = *idle_threads_iterator;
                    if (!dispatcher.hasNextStaticChunk(tid_to_run)) {
                        idle_threads_iterator++;
                        continue;
                    }
                    Chunk c = dispatcher.getNextStaticChunk(tid_to_run);
                    vector<int> parallel_iteration_vector;
                    parallel_iteration_vector.push_back(c.first);
                    parallel_iteration_vector.push_back(0);
                    if (progress[tid_to_run]) {
                        progress[tid_to_run]->ref = "C0";
                        progress[tid_to_run]->iteration = parallel_iteration_vector;
                        progress[tid_to_run]->chunk = c;
                    } else {
                        Progress *p = new Progress("C0", parallel_iteration_vector, c);
                        progress[tid_to_run] = p;
                    }
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
                    worker_threads.emplace_back(tid_to_run);
#if defined(DEBUG)
                    cout << "Move " << tid_to_run << " to worker_threads" << endl;
#endif
                } /* end of progress assignment */
            } /* end of chunk availability check */
INTERLEAVING_LOOP_0:
            //UNIFORM INTERLEAVING
            working_threads = (int)worker_threads.size();
            sort(worker_threads.begin(), worker_threads.end());
            auto worker_thread_iterator = worker_threads.begin();
            while (worker_thread_iterator != worker_threads.end()) {
                tid_to_run = *worker_thread_iterator;
                if (!progress[tid_to_run]->isInBound()) {
                    worker_thread_iterator++;
#if defined(DEBUG)
                    cout << "[" << tid_to_run << "] Out of bound" << endl;;
#endif
                    continue;
                }
#if defined(DEBUG)
                if (start_reuse_search) {
                    cout << "[" << tid_to_run << "] Access " << progress[tid_to_run]->ref << " {";
                    for (unsigned i = 0; i < progress[tid_to_run]->iteration.size(); i++) {
                        cout << progress[tid_to_run]->iteration[i];
                        if (i < progress[tid_to_run]->iteration.size()-1) { cout << ","; }
                    }
                    cout << "} " << start_reuse_search << endl;
                }
#endif
                if (progress[tid_to_run]->ref == "C0") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C0(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C2;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C1");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C0 */
                if (progress[tid_to_run]->ref == "C1") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C1(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C2;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    //CASE 2
                    progress[tid_to_run]->iteration.emplace_back(0);
                    progress[tid_to_run]->increment("A0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C1 */
                if (progress[tid_to_run]->ref == "A0") {
                    if (start_reuse_search) {
                        //A[c0][c2] will not access C
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("B0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to A0 */
                if (progress[tid_to_run]->ref == "B0") {
                    if (start_reuse_search) {
                        //B[c2][c1] will not access C
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C2");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to B0 */
                if (progress[tid_to_run]->ref == "C2") {
                    if (!start_reuse_search) { start_reuse_search= (start_tid == tid_to_run); }
                    if (start_reuse_search) {
                        Iteration *access = new Iteration("C2", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1],progress[tid_to_run]->iteration[2]});
                        string array = "C";
                        addr = GetAddress_C2(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
#if defined(DEBUG)
                        cout << access->toString() << " @ " << addr << endl;
#endif
                        if (access->compare(start_iteration) == 0) {
#if defined(DEBUG)
                            cout << "Meet the start sample " << access->toString() << endl;
#endif
                            isSample = true;
                        } else if ((!samples.empty() && access->compare(samples.top()) == 0) || (sample_names.find(access->toAddrString()) != sample_names.end())) {
#if defined(DEBUG)
                            cout << "Meet a new sample " << access->toString() << " while searching reuses" << endl;
#endif
                            traversing_samples++;
                            if (!samples.empty() && access->compare(samples.top()) == 0) {
                                samples.pop();
                            }
                            isSample = true;
                            samples_meet.insert(access->toAddrString());
                        }
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
#if defined(DEBUG)
                            cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
#endif
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C2;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        if (isSample) {
                            LAT[tid_to_run][addr] = count[tid_to_run];
                            LATSampleIterMap[tid_to_run][count[tid_to_run]] = access;
                        } else { delete access; }
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C3");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C2 */
                if (progress[tid_to_run]->ref == "C3") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C3(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C2;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    //CASE 3
                    if ((progress[tid_to_run]->iteration[2] + 1)<128) {
                        progress[tid_to_run]->iteration[2] = (progress[tid_to_run]->iteration[2] + 1);
                        progress[tid_to_run]->increment("A0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    if ((progress[tid_to_run]->iteration[1] + 1)<128) {
                        progress[tid_to_run]->iteration[1] = (progress[tid_to_run]->iteration[1] + 1);
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    progress[tid_to_run]->iteration[0] = (progress[tid_to_run]->iteration[0] + 1);
                    if (progress[tid_to_run]->isInBound()) {
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.emplace_back(0);
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    if (find(idle_threads.begin(), idle_threads.end(), tid_to_run) == idle_threads.end()) {
#if defined(DEBUG)
                        cout << "Move " << tid_to_run << " to idle threads" << endl;
#endif
                        idle_threads.emplace_back(tid_to_run);
                        worker_thread_iterator = worker_threads.erase(find(worker_threads.begin(), worker_threads.end(), tid_to_run));
                    }
                }
            } /* end of thread interleaving loop */
            if (idle_threads.size() == THREAD_NUM && !dispatcher.hasNextChunk(1)) {
                break;
            } /* end of break condition check */
        } /* end of while(true) */
        //reset both lists so they can be reused for later parallel loop
        idle_threads.clear();
        worker_threads.clear();
        for (unsigned i = 0; i < progress.size(); i++) {
            if (progress[i]) {
                delete progress[i];
                progress[i] = nullptr;
            }
        } // end of progress traversal
        goto END_SAMPLE;
    } // end of while (!samples.empty())
END_SAMPLE:
    if (!idle_threads.empty()) { idle_threads.clear(); };
    if (!worker_threads.empty()) { worker_threads.clear(); };
    if (!LAT.empty()) { ;
        for (unsigned i = 0; i < LAT.size(); i++) {
            pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
            LAT.clear();;
        }
    }
    for (unsigned i = 0; i < progress.size(); i++) {
        if (progress[i]) {
            delete progress[i];
            progress[i] = nullptr;
        }
    }
    if (!LATSampleIterMap.empty()) {
        for (auto LATSampleIterMapEntry : LATSampleIterMap) {
            for (auto entry : LATSampleIterMapEntry.second) {
                delete entry.second;
            }
            LATSampleIterMapEntry.second.clear();
        }
        LATSampleIterMap.clear();
    }
    no_share_distribute(no_share_histogram,histogram);
    share_distribute(share_histogram,histogram);
    no_share_histogram.clear();
    share_histogram.clear();
    iteration_traversed_map["C2"] = accumulate(count.begin(), count.end(), 0);
    return;
} // end of the sampler function
//A[c0][c2]
void sampler_A0(unordered_map<long, double> &histogram) {
    //Declare components will be used in Parallel RI search
    array<Progress *, THREAD_NUM> progress = { nullptr };
    vector<int> idle_threads, worker_threads, subscripts;
    ChunkDispatcher dispatcher;
    int tid_to_run = 0, start_tid = 0, working_threads = THREAD_NUM;
    uint64_t addr = 0;
    //global counter
    array<long, THREAD_NUM> count = {0};
    unordered_map<int, unordered_map<uint64_t, long>> LAT;
    unordered_map<int, unordered_map<long, Iteration *>> LATSampleIterMap;
    unordered_map<long, double> no_share_histogram;
    unordered_map<int, unordered_map<long, double>> share_histogram;
    unsigned s=0, traversing_samples=0;
    bool start_reuse_search=false;
    Iteration start_iteration;
    priority_queue<Iteration, vector<Iteration>, IterationComp> samples;
    unordered_set<string> sample_names, samples_meet;
    hash<std::thread::id> hasher;
    static thread_local mt19937 generator(time(NULL)+hasher(this_thread::get_id()));
    //generate all samples and sort them in order
    int num_samples = 2098;
    while (s < num_samples) {
        string sample_label = "A0";
        int sample_c0 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c0: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c0 << endl;
#endif
        if (((128-0)/1-((128-0)%1==0))< 0) { continue; }
        int sample_c1 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c1: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c1 << endl;
#endif
        if (((128-0)/1-((128-0)%1==0))< 0) { continue; }
        int sample_c2 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c2: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c2 << endl;
#endif
        sample_label += (to_string(sample_c0) + "-");
        sample_label += (to_string(sample_c1) + "-");
        sample_label += (to_string(sample_c2) + "-");
        //avoid duplicate sampling
        if (sample_names.find(sample_label)!=sample_names.end()) { continue; }
        sample_names.insert(sample_label);
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
        Iteration sample("A0", {sample_c0,sample_c1,sample_c2},0,1,true,0);
        samples.push(sample);
        s++;
    } // the end of while(s < NUM_SAMPLES)
    while (!samples.empty()) {
START_SAMPLE_A0:
        start_reuse_search = false;
        start_iteration = samples.top();
        traversing_samples++;
        samples.pop();
        if (!idle_threads.empty()) { idle_threads.clear(); };
        if (!worker_threads.empty()) { worker_threads.clear(); };
        if (!LAT.empty()) { ;
            for (unsigned i = 0; i < LAT.size(); i++) {
                pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
                LAT.clear();;
            }
        }
        if (!LATSampleIterMap.empty()) {
            for (auto LATSampleIterMapEntry : LATSampleIterMap) {
                for (auto entry : LATSampleIterMapEntry.second) {
                    delete entry.second;
                }
                LATSampleIterMapEntry.second.clear();
            }
            LATSampleIterMap.clear();
        }
#if defined(DEBUG)
        cout << "Start tracking sample " << start_iteration.toString() << endl;
#endif
        int sample_c0_start=start_iteration.ivs[0];
        int c0_start=0;
        int sample_c1_start=start_iteration.ivs[1];
        int c1_start=0;
        int sample_c2_start=start_iteration.ivs[2];
        int c2_start=0;
        //Generate parallel code for (c0,0,\<,128,(c0 + 1))
        //Threads are scheduled using static scheduling
        //Threads are interleaved using uniform interleaving
        dispatcher = ChunkDispatcher(CHUNK_SIZE,((128-0)),0,1);
        dispatcher.reset();
        idle_threads.clear();
        worker_threads.clear();
        for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
#if defined(DEBUG)
            cout << "Move " << tid_to_run << " to idle_threads" << endl;
#endif
            idle_threads.emplace_back(tid_to_run);
        }
        bool start_parallel_sample = true;
        if (start_parallel_sample) {
            dispatcher.setStartPoint(start_iteration.ivs[0]);
            start_tid = dispatcher.getStaticTid(start_iteration.ivs[0]);
            for (auto idle_threads_iterator = idle_threads.begin();idle_threads_iterator != idle_threads.end();) {
                int idle_tid = *idle_threads_iterator;
                if (dispatcher.hasNextStaticChunk(idle_tid)) {
                    Chunk start_chunk = dispatcher.getStaticStartChunk(start_iteration.ivs[0], idle_tid);
                    if (start_chunk.first > start_chunk.second) { 
                        //though the idle_tid has an available chunk, but at the given sampling point, there is no iteration to execute for this idle_tid
#if defined(DEBUG)
                        cout << "[" << start_chunk.first << ">" << start_chunk.second << "] No available iteration for " << idle_tid <<  endl;
#endif
                        idle_threads_iterator++;
                        continue;
                    }
                    if (progress[idle_tid]) {
                        progress[idle_tid]->ref = "A0";
                        progress[idle_tid]->iteration = {start_chunk.first, start_iteration.ivs[1], start_iteration.ivs[2], };
                        progress[idle_tid]->chunk = start_chunk;
                    } else {
                        Progress *p = new Progress("A0", {start_chunk.first, start_iteration.ivs[1], start_iteration.ivs[2], }, start_chunk);
                        progress[idle_tid] = p;
                    }
#if defined(DEBUG)
                    cout << "[" << idle_tid << "] starts from " << progress[idle_tid]->toString() << endl;
#endif
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
#if defined(DEBUG)
                    cout << "Move " << idle_tid << " to worker_threads" << endl;
#endif
                    worker_threads.push_back(idle_tid);
                    continue;
                }
                idle_threads_iterator++;
            }
            start_parallel_sample = false;
            //assign the start sample to each thread
            if (!worker_threads.empty()) {
                goto INTERLEAVING_LOOP_0;
            } else {
                goto END_SAMPLE;
            }
        } // end of if (start_parallel_sample)
        while(true) {
            if (!idle_threads.empty() && dispatcher.hasNextChunk(1)) {
                auto idle_threads_iterator = idle_threads.begin();
                while(idle_threads_iterator != idle_threads.end() && dispatcher.hasNextChunk(1)) {
                    tid_to_run = *idle_threads_iterator;
                    if (!dispatcher.hasNextStaticChunk(tid_to_run)) {
                        idle_threads_iterator++;
                        continue;
                    }
                    Chunk c = dispatcher.getNextStaticChunk(tid_to_run);
                    vector<int> parallel_iteration_vector;
                    parallel_iteration_vector.push_back(c.first);
                    parallel_iteration_vector.push_back(0);
                    if (progress[tid_to_run]) {
                        progress[tid_to_run]->ref = "C0";
                        progress[tid_to_run]->iteration = parallel_iteration_vector;
                        progress[tid_to_run]->chunk = c;
                    } else {
                        Progress *p = new Progress("C0", parallel_iteration_vector, c);
                        progress[tid_to_run] = p;
                    }
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
                    worker_threads.emplace_back(tid_to_run);
#if defined(DEBUG)
                    cout << "Move " << tid_to_run << " to worker_threads" << endl;
#endif
                } /* end of progress assignment */
            } /* end of chunk availability check */
INTERLEAVING_LOOP_0:
            //UNIFORM INTERLEAVING
            working_threads = (int)worker_threads.size();
            sort(worker_threads.begin(), worker_threads.end());
            auto worker_thread_iterator = worker_threads.begin();
            while (worker_thread_iterator != worker_threads.end()) {
                tid_to_run = *worker_thread_iterator;
                if (!progress[tid_to_run]->isInBound()) {
                    worker_thread_iterator++;
#if defined(DEBUG)
                    cout << "[" << tid_to_run << "] Out of bound" << endl;;
#endif
                    continue;
                }
#if defined(DEBUG)
                if (start_reuse_search) {
                    cout << "[" << tid_to_run << "] Access " << progress[tid_to_run]->ref << " {";
                    for (unsigned i = 0; i < progress[tid_to_run]->iteration.size(); i++) {
                        cout << progress[tid_to_run]->iteration[i];
                        if (i < progress[tid_to_run]->iteration.size()-1) { cout << ","; }
                    }
                    cout << "} " << start_reuse_search << endl;
                }
#endif
                if (progress[tid_to_run]->ref == "C0") {
                    if (start_reuse_search) {
                        //C[c0][c1] will not access A
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C1");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C0 */
                if (progress[tid_to_run]->ref == "C1") {
                    if (start_reuse_search) {
                        //C[c0][c1] will not access A
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    //CASE 2
                    progress[tid_to_run]->iteration.emplace_back(0);
                    progress[tid_to_run]->increment("A0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C1 */
                if (progress[tid_to_run]->ref == "A0") {
                    if (!start_reuse_search) { start_reuse_search= (start_tid == tid_to_run); }
                    if (start_reuse_search) {
                        Iteration *access = new Iteration("A0", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1],progress[tid_to_run]->iteration[2]});
                        string array = "A";
                        addr = GetAddress_A0(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[2]);
                        bool isSample = false;
#if defined(DEBUG)
                        cout << access->toString() << " @ " << addr << endl;
#endif
                        if (access->compare(start_iteration) == 0) {
#if defined(DEBUG)
                            cout << "Meet the start sample " << access->toString() << endl;
#endif
                            isSample = true;
                        } else if ((!samples.empty() && access->compare(samples.top()) == 0) || (sample_names.find(access->toAddrString()) != sample_names.end())) {
#if defined(DEBUG)
                            cout << "Meet a new sample " << access->toString() << " while searching reuses" << endl;
#endif
                            traversing_samples++;
                            if (!samples.empty() && access->compare(samples.top()) == 0) {
                                samples.pop();
                            }
                            isSample = true;
                            samples_meet.insert(access->toAddrString());
                        }
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
#if defined(DEBUG)
                            cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
#endif
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_A0;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        if (isSample) {
                            LAT[tid_to_run][addr] = count[tid_to_run];
                            LATSampleIterMap[tid_to_run][count[tid_to_run]] = access;
                        } else { delete access; }
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("B0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to A0 */
                if (progress[tid_to_run]->ref == "B0") {
                    if (start_reuse_search) {
                        //B[c2][c1] will not access A
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C2");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to B0 */
                if (progress[tid_to_run]->ref == "C2") {
                    if (start_reuse_search) {
                        //C[c0][c1] will not access A
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C3");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C2 */
                if (progress[tid_to_run]->ref == "C3") {
                    if (start_reuse_search) {
                        //C[c0][c1] will not access A
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    //CASE 3
                    if ((progress[tid_to_run]->iteration[2] + 1)<128) {
                        progress[tid_to_run]->iteration[2] = (progress[tid_to_run]->iteration[2] + 1);
                        progress[tid_to_run]->increment("A0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    if ((progress[tid_to_run]->iteration[1] + 1)<128) {
                        progress[tid_to_run]->iteration[1] = (progress[tid_to_run]->iteration[1] + 1);
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    progress[tid_to_run]->iteration[0] = (progress[tid_to_run]->iteration[0] + 1);
                    if (progress[tid_to_run]->isInBound()) {
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.emplace_back(0);
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    if (find(idle_threads.begin(), idle_threads.end(), tid_to_run) == idle_threads.end()) {
#if defined(DEBUG)
                        cout << "Move " << tid_to_run << " to idle threads" << endl;
#endif
                        idle_threads.emplace_back(tid_to_run);
                        worker_thread_iterator = worker_threads.erase(find(worker_threads.begin(), worker_threads.end(), tid_to_run));
                    }
                }
            } /* end of thread interleaving loop */
            if (idle_threads.size() == THREAD_NUM && !dispatcher.hasNextChunk(1)) {
                break;
            } /* end of break condition check */
        } /* end of while(true) */
        //reset both lists so they can be reused for later parallel loop
        idle_threads.clear();
        worker_threads.clear();
        for (unsigned i = 0; i < progress.size(); i++) {
            if (progress[i]) {
                delete progress[i];
                progress[i] = nullptr;
            }
        } // end of progress traversal
        goto END_SAMPLE; // TODO: do not understand why do we have a while (!samples.empty()) loop when we just jump to END_SAMPLE in the first iteration???
    } // end of while (!samples.empty())
END_SAMPLE:
    if (!idle_threads.empty()) { idle_threads.clear(); };
    if (!worker_threads.empty()) { worker_threads.clear(); };
    if (!LAT.empty()) { ;
        for (unsigned i = 0; i < LAT.size(); i++) {
            pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
            LAT.clear();;
        }
    }
    for (unsigned i = 0; i < progress.size(); i++) {
        if (progress[i]) {
            delete progress[i];
            progress[i] = nullptr;
        }
    }
    if (!LATSampleIterMap.empty()) {
        for (auto LATSampleIterMapEntry : LATSampleIterMap) {
            for (auto entry : LATSampleIterMapEntry.second) {
                delete entry.second;
            }
            LATSampleIterMapEntry.second.clear();
        }
        LATSampleIterMap.clear();
    }
    no_share_distribute(no_share_histogram,histogram);
    share_distribute(share_histogram,histogram);
    no_share_histogram.clear();
    share_histogram.clear();
    iteration_traversed_map["A0"] = accumulate(count.begin(), count.end(), 0);
    return;
} // end of the sampler function
//C[c0][c1]
void sampler_C0(unordered_map<long, double> &histogram) {
    //Declare components will be used in Parallel RI search
    array<Progress *, THREAD_NUM> progress = { nullptr };
    vector<int> idle_threads, worker_threads, subscripts;
    ChunkDispatcher dispatcher;
    int tid_to_run = 0, start_tid = 0, working_threads = THREAD_NUM;
    uint64_t addr = 0;
    //global counter
    array<long, THREAD_NUM> count = {0};
    unordered_map<int, unordered_map<uint64_t, long>> LAT;
    unordered_map<int, unordered_map<long, Iteration *>> LATSampleIterMap;
    unordered_map<long, double> no_share_histogram;
    unordered_map<int, unordered_map<long, double>> share_histogram;
    unsigned s=0, traversing_samples=0;
    bool start_reuse_search=false;
    Iteration start_iteration;
    priority_queue<Iteration, vector<Iteration>, IterationComp> samples;
    unordered_set<string> sample_names, samples_meet;
    hash<std::thread::id> hasher;
    static thread_local mt19937 generator(time(NULL)+hasher(this_thread::get_id()));
    //generate all samples and sort them in order
    int num_samples = 164;
    while (s < num_samples) {
        string sample_label = "C0";
        int sample_c0 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c0: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c0 << endl;
#endif
        if (((128-0)/1-((128-0)%1==0))< 0) { continue; }
        int sample_c1 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c1: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c1 << endl;
#endif
        sample_label += (to_string(sample_c0) + "-");
        sample_label += (to_string(sample_c1) + "-");
        //avoid duplicate sampling
        if (sample_names.find(sample_label)!=sample_names.end()) { continue; }
        sample_names.insert(sample_label);
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
        Iteration sample("C0", {sample_c0,sample_c1},0,1,true,0);
        samples.push(sample);
        s++;
    } // the end of while(s < NUM_SAMPLES)
    while (!samples.empty()) {
START_SAMPLE_C0:
        start_reuse_search = false;
        start_iteration = samples.top();
        traversing_samples++;
        samples.pop();
        if (!idle_threads.empty()) { idle_threads.clear(); };
        if (!worker_threads.empty()) { worker_threads.clear(); };
        if (!LAT.empty()) { ;
            for (unsigned i = 0; i < LAT.size(); i++) {
                pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
                LAT.clear();;
            }
        }
        if (!LATSampleIterMap.empty()) {
            for (auto LATSampleIterMapEntry : LATSampleIterMap) {
                for (auto entry : LATSampleIterMapEntry.second) {
                    delete entry.second;
                }
                LATSampleIterMapEntry.second.clear();
            }
            LATSampleIterMap.clear();
        }
#if defined(DEBUG)
        cout << "Start tracking sample " << start_iteration.toString() << endl;
#endif
        int sample_c0_start=start_iteration.ivs[0];
        int c0_start=0;
        int sample_c1_start=start_iteration.ivs[1];
        int c1_start=0;
        //Generate parallel code for (c0,0,\<,128,(c0 + 1))
        //Threads are scheduled using static scheduling
        //Threads are interleaved using uniform interleaving
        dispatcher = ChunkDispatcher(CHUNK_SIZE,((128-0)),0,1);
        dispatcher.reset();
        idle_threads.clear();
        worker_threads.clear();
        for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
#if defined(DEBUG)
            cout << "Move " << tid_to_run << " to idle_threads" << endl;
#endif
            idle_threads.emplace_back(tid_to_run);
        }
        bool start_parallel_sample = true;
        if (start_parallel_sample) {
            dispatcher.setStartPoint(start_iteration.ivs[0]);
            start_tid = dispatcher.getStaticTid(start_iteration.ivs[0]);
            for (auto idle_threads_iterator = idle_threads.begin();idle_threads_iterator != idle_threads.end();) {
                int idle_tid = *idle_threads_iterator;
                if (dispatcher.hasNextStaticChunk(idle_tid)) {
                    Chunk start_chunk = dispatcher.getStaticStartChunk(start_iteration.ivs[0], idle_tid);
                    if (start_chunk.first > start_chunk.second) { 
                        //though the idle_tid has an available chunk, but at the given sampling point, there is no iteration to execute for this idle_tid
#if defined(DEBUG)
                        cout << "[" << start_chunk.first << ">" << start_chunk.second << "] No available iteration for " << idle_tid <<  endl;
#endif
                        idle_threads_iterator++;
                        continue;
                    }
                    if (progress[idle_tid]) {
                        progress[idle_tid]->ref = "C0";
                        progress[idle_tid]->iteration = {start_chunk.first, start_iteration.ivs[1], };
                        progress[idle_tid]->chunk = start_chunk;
                    } else {
                        Progress *p = new Progress("C0", {start_chunk.first, start_iteration.ivs[1], }, start_chunk);
                        progress[idle_tid] = p;
                    }
#if defined(DEBUG)
                    cout << "[" << idle_tid << "] starts from " << progress[idle_tid]->toString() << endl;
#endif
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
#if defined(DEBUG)
                    cout << "Move " << idle_tid << " to worker_threads" << endl;
#endif
                    worker_threads.push_back(idle_tid);
                    continue;
                }
                idle_threads_iterator++;
            }
            start_parallel_sample = false;
            //assign the start sample to each thread
            if (!worker_threads.empty()) {
                goto INTERLEAVING_LOOP_0;
            } else {
                goto END_SAMPLE;
            }
        } // end of if (start_parallel_sample)
        while(true) {
            if (!idle_threads.empty() && dispatcher.hasNextChunk(1)) {
                auto idle_threads_iterator = idle_threads.begin();
                while(idle_threads_iterator != idle_threads.end() && dispatcher.hasNextChunk(1)) {
                    tid_to_run = *idle_threads_iterator;
                    if (!dispatcher.hasNextStaticChunk(tid_to_run)) {
                        idle_threads_iterator++;
                        continue;
                    }
                    Chunk c = dispatcher.getNextStaticChunk(tid_to_run);
                    vector<int> parallel_iteration_vector;
                    parallel_iteration_vector.push_back(c.first);
                    parallel_iteration_vector.push_back(0);
                    if (progress[tid_to_run]) {
                        progress[tid_to_run]->ref = "C0";
                        progress[tid_to_run]->iteration = parallel_iteration_vector;
                        progress[tid_to_run]->chunk = c;
                    } else {
                        Progress *p = new Progress("C0", parallel_iteration_vector, c);
                        progress[tid_to_run] = p;
                    }
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
                    worker_threads.emplace_back(tid_to_run);
#if defined(DEBUG)
                    cout << "Move " << tid_to_run << " to worker_threads" << endl;
#endif
                } /* end of progress assignment */
            } /* end of chunk availability check */
INTERLEAVING_LOOP_0:
            //UNIFORM INTERLEAVING
            working_threads = (int)worker_threads.size();
            sort(worker_threads.begin(), worker_threads.end());
            auto worker_thread_iterator = worker_threads.begin();
            while (worker_thread_iterator != worker_threads.end()) {
                tid_to_run = *worker_thread_iterator;
                if (!progress[tid_to_run]->isInBound()) {
                    worker_thread_iterator++;
#if defined(DEBUG)
                    cout << "[" << tid_to_run << "] Out of bound" << endl;;
#endif
                    continue;
                }
#if defined(DEBUG)
                if (start_reuse_search) {
                    cout << "[" << tid_to_run << "] Access " << progress[tid_to_run]->ref << " {";
                    for (unsigned i = 0; i < progress[tid_to_run]->iteration.size(); i++) {
                        cout << progress[tid_to_run]->iteration[i];
                        if (i < progress[tid_to_run]->iteration.size()-1) { cout << ","; }
                    }
                    cout << "} " << start_reuse_search << endl;
                }
#endif
                if (progress[tid_to_run]->ref == "C0") {
                    if (!start_reuse_search) { start_reuse_search= (start_tid == tid_to_run); }
                    if (start_reuse_search) {
                        Iteration *access = new Iteration("C0", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
                        string array = "C";
                        addr = GetAddress_C0(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
#if defined(DEBUG)
                        cout << access->toString() << " @ " << addr << endl;
#endif
                        if (access->compare(start_iteration) == 0) {
#if defined(DEBUG)
                            cout << "Meet the start sample " << access->toString() << endl;
#endif
                            isSample = true;
                        } else if ((!samples.empty() && access->compare(samples.top()) == 0) || (sample_names.find(access->toAddrString()) != sample_names.end())) {
#if defined(DEBUG)
                            cout << "Meet a new sample " << access->toString() << " while searching reuses" << endl;
#endif
                            traversing_samples++;
                            if (!samples.empty() && access->compare(samples.top()) == 0) {
                                samples.pop();
                            }
                            isSample = true;
                            samples_meet.insert(access->toAddrString());
                        }
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
#if defined(DEBUG)
                            cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
#endif
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C0;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        if (isSample) {
                            LAT[tid_to_run][addr] = count[tid_to_run];
                            LATSampleIterMap[tid_to_run][count[tid_to_run]] = access;
                        } else { delete access; }
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C1");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C0 */
                if (progress[tid_to_run]->ref == "C1") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C1(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C0;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    //CASE 2
                    progress[tid_to_run]->iteration.emplace_back(0);
                    progress[tid_to_run]->increment("A0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C1 */
                if (progress[tid_to_run]->ref == "A0") {
                    if (start_reuse_search) {
                        //A[c0][c2] will not access C
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("B0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to A0 */
                if (progress[tid_to_run]->ref == "B0") {
                    if (start_reuse_search) {
                        //B[c2][c1] will not access C
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C2");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to B0 */
                if (progress[tid_to_run]->ref == "C2") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C2(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C0;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C3");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C2 */
                if (progress[tid_to_run]->ref == "C3") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C3(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C0;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    //CASE 3
                    if ((progress[tid_to_run]->iteration[2] + 1)<128) {
                        progress[tid_to_run]->iteration[2] = (progress[tid_to_run]->iteration[2] + 1);
                        progress[tid_to_run]->increment("A0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    if ((progress[tid_to_run]->iteration[1] + 1)<128) {
                        progress[tid_to_run]->iteration[1] = (progress[tid_to_run]->iteration[1] + 1);
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    progress[tid_to_run]->iteration[0] = (progress[tid_to_run]->iteration[0] + 1);
                    if (progress[tid_to_run]->isInBound()) {
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.emplace_back(0);
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    if (find(idle_threads.begin(), idle_threads.end(), tid_to_run) == idle_threads.end()) {
#if defined(DEBUG)
                        cout << "Move " << tid_to_run << " to idle threads" << endl;
#endif
                        idle_threads.emplace_back(tid_to_run);
                        worker_thread_iterator = worker_threads.erase(find(worker_threads.begin(), worker_threads.end(), tid_to_run));
                    }
                }
            } /* end of thread interleaving loop */
            if (idle_threads.size() == THREAD_NUM && !dispatcher.hasNextChunk(1)) {
                break;
            } /* end of break condition check */
        } /* end of while(true) */
        //reset both lists so they can be reused for later parallel loop
        idle_threads.clear();
        worker_threads.clear();
        for (unsigned i = 0; i < progress.size(); i++) {
            if (progress[i]) {
                delete progress[i];
                progress[i] = nullptr;
            }
        } // end of progress traversal
        goto END_SAMPLE;
    } // end of while (!samples.empty())
END_SAMPLE:
    if (!idle_threads.empty()) { idle_threads.clear(); };
    if (!worker_threads.empty()) { worker_threads.clear(); };
    if (!LAT.empty()) { ;
        for (unsigned i = 0; i < LAT.size(); i++) {
            pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
            LAT.clear();;
        }
    }
    for (unsigned i = 0; i < progress.size(); i++) {
        if (progress[i]) {
            delete progress[i];
            progress[i] = nullptr;
        }
    }
    if (!LATSampleIterMap.empty()) {
        for (auto LATSampleIterMapEntry : LATSampleIterMap) {
            for (auto entry : LATSampleIterMapEntry.second) {
                delete entry.second;
            }
            LATSampleIterMapEntry.second.clear();
        }
        LATSampleIterMap.clear();
    }
    no_share_distribute(no_share_histogram,histogram);
    share_distribute(share_histogram,histogram);
    no_share_histogram.clear();
    share_histogram.clear();
    iteration_traversed_map["C0"] = accumulate(count.begin(), count.end(), 0);
    return;
} // end of the sampler function
//B[c2][c1]
void sampler_B0(unordered_map<long, double> &histogram) {
    //Declare components will be used in Parallel RI search
    array<Progress *, THREAD_NUM> progress = { nullptr };
    vector<int> idle_threads, worker_threads, subscripts;
    ChunkDispatcher dispatcher;
    int tid_to_run = 0, start_tid = 0, working_threads = THREAD_NUM;
    uint64_t addr = 0;
    //global counter
    array<long, THREAD_NUM> count = {0};
    unordered_map<int, unordered_map<uint64_t, long>> LAT;
    unordered_map<int, unordered_map<long, Iteration *>> LATSampleIterMap;
    unordered_map<long, double> no_share_histogram;
    unordered_map<int, unordered_map<long, double>> share_histogram;
    unsigned s=0, traversing_samples=0;
    bool start_reuse_search=false;
    Iteration start_iteration;
    priority_queue<Iteration, vector<Iteration>, IterationComp> samples;
    unordered_set<string> sample_names, samples_meet;
    hash<std::thread::id> hasher;
    static thread_local mt19937 generator(time(NULL)+hasher(this_thread::get_id()));
    //generate all samples and sort them in order
    int num_samples = 2098; // TODO: do not understand why this is 2098 and others are 164
    while (s < num_samples) {
        string sample_label = "B0";
        int sample_c0 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c0: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c0 << endl;
#endif
        if (((128-0)/1-((128-0)%1==0))< 0) { continue; }
        int sample_c1 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c1: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c1 << endl;
#endif
        if (((128-0)/1-((128-0)%1==0))< 0) { continue; }
        int sample_c2 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c2: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c2 << endl;
#endif
        sample_label += (to_string(sample_c0) + "-");
        sample_label += (to_string(sample_c1) + "-");
        sample_label += (to_string(sample_c2) + "-"); //TODO: do not understand why B0 never concatenate A0 but can concatenate C2
        //avoid duplicate sampling
        if (sample_names.find(sample_label)!=sample_names.end()) { continue; }
        sample_names.insert(sample_label);
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
        Iteration sample("B0", {sample_c0,sample_c1,sample_c2},0,1,true,0);
        samples.push(sample);
        s++;
    } // the end of while(s < NUM_SAMPLES)
    while (!samples.empty()) {
START_SAMPLE_B0:
        start_reuse_search = false;
        start_iteration = samples.top();
        traversing_samples++;
        samples.pop();
        if (!idle_threads.empty()) { idle_threads.clear(); };
        if (!worker_threads.empty()) { worker_threads.clear(); };
        if (!LAT.empty()) { ;
            for (unsigned i = 0; i < LAT.size(); i++) {
                pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
                LAT.clear();;
            }
        }
        if (!LATSampleIterMap.empty()) {
            for (auto LATSampleIterMapEntry : LATSampleIterMap) {
                for (auto entry : LATSampleIterMapEntry.second) {
                    delete entry.second;
                }
                LATSampleIterMapEntry.second.clear();
            }
            LATSampleIterMap.clear();
        }
#if defined(DEBUG)
        cout << "Start tracking sample " << start_iteration.toString() << endl;
#endif
        int sample_c0_start=start_iteration.ivs[0];
        int c0_start=0;
        int sample_c1_start=start_iteration.ivs[1];
        int c1_start=0;
        int sample_c2_start=start_iteration.ivs[2];
        int c2_start=0;
        //Generate parallel code for (c0,0,\<,128,(c0 + 1))
        //Threads are scheduled using static scheduling
        //Threads are interleaved using uniform interleaving
        dispatcher = ChunkDispatcher(CHUNK_SIZE,((128-0)),0,1);
        dispatcher.reset();
        idle_threads.clear();
        worker_threads.clear();
        for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
#if defined(DEBUG)
            cout << "Move " << tid_to_run << " to idle_threads" << endl;
#endif
            idle_threads.emplace_back(tid_to_run);
        }
        bool start_parallel_sample = true;
        if (start_parallel_sample) {
            dispatcher.setStartPoint(start_iteration.ivs[0]);
            start_tid = dispatcher.getStaticTid(start_iteration.ivs[0]);
            for (auto idle_threads_iterator = idle_threads.begin();idle_threads_iterator != idle_threads.end();) {
                int idle_tid = *idle_threads_iterator;
                if (dispatcher.hasNextStaticChunk(idle_tid)) {
                    Chunk start_chunk = dispatcher.getStaticStartChunk(start_iteration.ivs[0], idle_tid);
                    if (start_chunk.first > start_chunk.second) { 
                        //though the idle_tid has an available chunk, but at the given sampling point, there is no iteration to execute for this idle_tid
#if defined(DEBUG)
                        cout << "[" << start_chunk.first << ">" << start_chunk.second << "] No available iteration for " << idle_tid <<  endl;
#endif
                        idle_threads_iterator++;
                        continue;
                    }
                    if (progress[idle_tid]) {
                        progress[idle_tid]->ref = "B0";
                        progress[idle_tid]->iteration = {start_chunk.first, start_iteration.ivs[1], start_iteration.ivs[2], };
                        progress[idle_tid]->chunk = start_chunk;
                    } else {
                        Progress *p = new Progress("B0", {start_chunk.first, start_iteration.ivs[1], start_iteration.ivs[2], }, start_chunk);
                        progress[idle_tid] = p;
                    }
#if defined(DEBUG)
                    cout << "[" << idle_tid << "] starts from " << progress[idle_tid]->toString() << endl;
#endif
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
#if defined(DEBUG)
                    cout << "Move " << idle_tid << " to worker_threads" << endl;
#endif
                    worker_threads.push_back(idle_tid);
                    continue;
                }
                idle_threads_iterator++;
            }
            start_parallel_sample = false;
            //assign the start sample to each thread
            if (!worker_threads.empty()) {
                goto INTERLEAVING_LOOP_0;
            } else {
                goto END_SAMPLE;
            }
        } // end of if (start_parallel_sample)
        while(true) {
            if (!idle_threads.empty() && dispatcher.hasNextChunk(1)) {
                auto idle_threads_iterator = idle_threads.begin();
                while(idle_threads_iterator != idle_threads.end() && dispatcher.hasNextChunk(1)) {
                    tid_to_run = *idle_threads_iterator;
                    if (!dispatcher.hasNextStaticChunk(tid_to_run)) {
                        idle_threads_iterator++;
                        continue;
                    }
                    Chunk c = dispatcher.getNextStaticChunk(tid_to_run);
                    vector<int> parallel_iteration_vector;
                    parallel_iteration_vector.push_back(c.first);
                    parallel_iteration_vector.push_back(0);
                    if (progress[tid_to_run]) {
                        progress[tid_to_run]->ref = "C0";
                        progress[tid_to_run]->iteration = parallel_iteration_vector;
                        progress[tid_to_run]->chunk = c;
                    } else {
                        Progress *p = new Progress("C0", parallel_iteration_vector, c);
                        progress[tid_to_run] = p;
                    }
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
                    worker_threads.emplace_back(tid_to_run);
#if defined(DEBUG)
                    cout << "Move " << tid_to_run << " to worker_threads" << endl;
#endif
                } /* end of progress assignment */
            } /* end of chunk availability check */
INTERLEAVING_LOOP_0:
            //UNIFORM INTERLEAVING
            working_threads = (int)worker_threads.size();
            sort(worker_threads.begin(), worker_threads.end());
            auto worker_thread_iterator = worker_threads.begin();
            while (worker_thread_iterator != worker_threads.end()) {
                tid_to_run = *worker_thread_iterator;
                if (!progress[tid_to_run]->isInBound()) {
                    worker_thread_iterator++;
#if defined(DEBUG)
                    cout << "[" << tid_to_run << "] Out of bound" << endl;;
#endif
                    continue;
                }
#if defined(DEBUG)
                if (start_reuse_search) {
                    cout << "[" << tid_to_run << "] Access " << progress[tid_to_run]->ref << " {";
                    for (unsigned i = 0; i < progress[tid_to_run]->iteration.size(); i++) {
                        cout << progress[tid_to_run]->iteration[i];
                        if (i < progress[tid_to_run]->iteration.size()-1) { cout << ","; }
                    }
                    cout << "} " << start_reuse_search << endl;
                }
#endif
                if (progress[tid_to_run]->ref == "C0") {
                    if (start_reuse_search) {
                        //C[c0][c1] will not access B
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C1");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C0 */
                if (progress[tid_to_run]->ref == "C1") {
                    if (start_reuse_search) {
                        //C[c0][c1] will not access B
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    //CASE 2
                    progress[tid_to_run]->iteration.emplace_back(0);
                    progress[tid_to_run]->increment("A0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C1 */
                if (progress[tid_to_run]->ref == "A0") {
                    if (start_reuse_search) {
                        //A[c0][c2] will not access B
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("B0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to A0 */
                if (progress[tid_to_run]->ref == "B0") {
                    if (!start_reuse_search) { start_reuse_search= (start_tid == tid_to_run); }
                    if (start_reuse_search) {
                        Iteration *access = new Iteration("B0", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1],progress[tid_to_run]->iteration[2]});
                        string array = "B";
                        addr = GetAddress_B0(progress[tid_to_run]->iteration[2],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
#if defined(DEBUG)
                        cout << access->toString() << " @ " << addr << endl;
#endif
                        if (access->compare(start_iteration) == 0) {
#if defined(DEBUG)
                            cout << "Meet the start sample " << access->toString() << endl;
#endif
                            isSample = true;
                        } else if ((!samples.empty() && access->compare(samples.top()) == 0) || (sample_names.find(access->toAddrString()) != sample_names.end())) {
#if defined(DEBUG)
                            cout << "Meet a new sample " << access->toString() << " while searching reuses" << endl;
#endif
                            traversing_samples++;
                            if (!samples.empty() && access->compare(samples.top()) == 0) {
                                samples.pop();
                            }
                            isSample = true;
                            samples_meet.insert(access->toAddrString());
                        }
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 0
                            //Sink has Parallel Induction : 0
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
/* Compare c2*/
/* With c2*/
/* With c1*/
/* Compare c1*/
/* With c2*/
/* With c1*/
                            //B[c2][c1] is carried by (c1,0,\<,128,(c1 + 1))
                            if (distance_to(reuse,0) > distance_to(reuse,(((4)*((128-0)/1)+2)*((128-0)/1)+0))) {
                                pluss_parallel_histogram_update(share_histogram[THREAD_NUM-1],reuse,1.);
                            } else {
                                pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            }
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
#if defined(DEBUG)
                            cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
#endif
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_B0;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        if (isSample) {
                            LAT[tid_to_run][addr] = count[tid_to_run];
                            LATSampleIterMap[tid_to_run][count[tid_to_run]] = access;
                        } else { delete access; }
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C2");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to B0 */
                if (progress[tid_to_run]->ref == "C2") {
                    if (start_reuse_search) {
                        //C[c0][c1] will not access B
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C3");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C2 */
                if (progress[tid_to_run]->ref == "C3") {
                    if (start_reuse_search) {
                        //C[c0][c1] will not access B
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    //CASE 3
                    if ((progress[tid_to_run]->iteration[2] + 1)<128) {
                        progress[tid_to_run]->iteration[2] = (progress[tid_to_run]->iteration[2] + 1);
                        progress[tid_to_run]->increment("A0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    if ((progress[tid_to_run]->iteration[1] + 1)<128) {
                        progress[tid_to_run]->iteration[1] = (progress[tid_to_run]->iteration[1] + 1);
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    progress[tid_to_run]->iteration[0] = (progress[tid_to_run]->iteration[0] + 1);
                    if (progress[tid_to_run]->isInBound()) {
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.emplace_back(0);
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    if (find(idle_threads.begin(), idle_threads.end(), tid_to_run) == idle_threads.end()) {
#if defined(DEBUG)
                        cout << "Move " << tid_to_run << " to idle threads" << endl;
#endif
                        idle_threads.emplace_back(tid_to_run);
                        worker_thread_iterator = worker_threads.erase(find(worker_threads.begin(), worker_threads.end(), tid_to_run));
                    }
                }
            } /* end of thread interleaving loop */
            if (idle_threads.size() == THREAD_NUM && !dispatcher.hasNextChunk(1)) {
                break;
            } /* end of break condition check */
        } /* end of while(true) */
        //reset both lists so they can be reused for later parallel loop
        idle_threads.clear();
        worker_threads.clear();
        for (unsigned i = 0; i < progress.size(); i++) {
            if (progress[i]) {
                delete progress[i];
                progress[i] = nullptr;
            }
        } // end of progress traversal
        goto END_SAMPLE;
    } // end of while (!samples.empty())
END_SAMPLE:
    if (!idle_threads.empty()) { idle_threads.clear(); };
    if (!worker_threads.empty()) { worker_threads.clear(); };
    if (!LAT.empty()) { ;
        for (unsigned i = 0; i < LAT.size(); i++) {
            pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
            LAT.clear();;
        }
    }
    for (unsigned i = 0; i < progress.size(); i++) {
        if (progress[i]) {
            delete progress[i];
            progress[i] = nullptr;
        }
    }
    if (!LATSampleIterMap.empty()) {
        for (auto LATSampleIterMapEntry : LATSampleIterMap) {
            for (auto entry : LATSampleIterMapEntry.second) {
                delete entry.second;
            }
            LATSampleIterMapEntry.second.clear();
        }
        LATSampleIterMap.clear();
    }
    no_share_distribute(no_share_histogram,histogram);
    share_distribute(share_histogram,histogram);
    no_share_histogram.clear();
    share_histogram.clear();
    iteration_traversed_map["B0"] = accumulate(count.begin(), count.end(), 0);
    return;
} // end of the sampler function
//C[c0][c1]
void sampler_C1(unordered_map<long, double> &histogram) {
    //Declare components will be used in Parallel RI search
    array<Progress *, THREAD_NUM> progress = { nullptr };
    vector<int> idle_threads, worker_threads, subscripts;
    ChunkDispatcher dispatcher;
    int tid_to_run = 0, start_tid = 0, working_threads = THREAD_NUM;
    uint64_t addr = 0;
    //global counter
    array<long, THREAD_NUM> count = {0};
    unordered_map<int, unordered_map<uint64_t, long>> LAT;
    unordered_map<int, unordered_map<long, Iteration *>> LATSampleIterMap;
    unordered_map<long, double> no_share_histogram;
    unordered_map<int, unordered_map<long, double>> share_histogram;
    unsigned s=0, traversing_samples=0;
    bool start_reuse_search=false;
    Iteration start_iteration;
    priority_queue<Iteration, vector<Iteration>, IterationComp> samples;
    unordered_set<string> sample_names, samples_meet;
    hash<std::thread::id> hasher;
    static thread_local mt19937 generator(time(NULL)+hasher(this_thread::get_id()));
    //generate all samples and sort them in order
    int num_samples = 164;
    while (s < num_samples) {
        string sample_label = "C1";
        int sample_c0 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c0: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c0 << endl;
#endif
        if (((128-0)/1-((128-0)%1==0))< 0) { continue; }
        int sample_c1 = (((128-0)/1-((128-0)%1==0)) == 0)?(0):(rand()%(((128-0)/1-((128-0)%1==0))))*(1)+(0);
#if defined(DEBUG)
        cout << "sample_c1: [" << (0) << "," << (((128-0)/1-((128-0)%1==0))*(1)+0) << "] samples " << sample_c1 << endl;
#endif
        sample_label += (to_string(sample_c0) + "-");
        sample_label += (to_string(sample_c1) + "-");
        //avoid duplicate sampling
        if (sample_names.find(sample_label)!=sample_names.end()) { continue; }
        sample_names.insert(sample_label);
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
        Iteration sample("C1", {sample_c0,sample_c1},0,1,true,0);
        samples.push(sample);
        s++;
    } // the end of while(s < NUM_SAMPLES)
    while (!samples.empty()) {
START_SAMPLE_C1:
        start_reuse_search = false;
        start_iteration = samples.top();
        traversing_samples++;
        samples.pop();
        if (!idle_threads.empty()) { idle_threads.clear(); };
        if (!worker_threads.empty()) { worker_threads.clear(); };
        if (!LAT.empty()) { ;
            for (unsigned i = 0; i < LAT.size(); i++) {
                pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
                LAT.clear();;
            }
        }
        if (!LATSampleIterMap.empty()) {
            for (auto LATSampleIterMapEntry : LATSampleIterMap) {
                for (auto entry : LATSampleIterMapEntry.second) {
                    delete entry.second;
                }
                LATSampleIterMapEntry.second.clear();
            }
            LATSampleIterMap.clear();
        }
#if defined(DEBUG)
        cout << "Start tracking sample " << start_iteration.toString() << endl;
#endif
        int sample_c0_start=start_iteration.ivs[0];
        int c0_start=0;
        int sample_c1_start=start_iteration.ivs[1];
        int c1_start=0;
        //Generate parallel code for (c0,0,\<,128,(c0 + 1))
        //Threads are scheduled using static scheduling
        //Threads are interleaved using uniform interleaving
        dispatcher = ChunkDispatcher(CHUNK_SIZE,((128-0)),0,1);
        dispatcher.reset();
        idle_threads.clear();
        worker_threads.clear();
        for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
#if defined(DEBUG)
            cout << "Move " << tid_to_run << " to idle_threads" << endl;
#endif
            idle_threads.emplace_back(tid_to_run);
        }
        bool start_parallel_sample = true;
        if (start_parallel_sample) {
            dispatcher.setStartPoint(start_iteration.ivs[0]);
            start_tid = dispatcher.getStaticTid(start_iteration.ivs[0]);
            for (auto idle_threads_iterator = idle_threads.begin();idle_threads_iterator != idle_threads.end();) {
                int idle_tid = *idle_threads_iterator;
                if (dispatcher.hasNextStaticChunk(idle_tid)) {
                    Chunk start_chunk = dispatcher.getStaticStartChunk(start_iteration.ivs[0], idle_tid);
                    if (start_chunk.first > start_chunk.second) { 
                        //though the idle_tid has an available chunk, but at the given sampling point, there is no iteration to execute for this idle_tid
#if defined(DEBUG)
                        cout << "[" << start_chunk.first << ">" << start_chunk.second << "] No available iteration for " << idle_tid <<  endl;
#endif
                        idle_threads_iterator++;
                        continue;
                    }
                    if (progress[idle_tid]) {
                        progress[idle_tid]->ref = "C1";
                        progress[idle_tid]->iteration = {start_chunk.first, start_iteration.ivs[1], };
                        progress[idle_tid]->chunk = start_chunk;
                    } else {
                        Progress *p = new Progress("C1", {start_chunk.first, start_iteration.ivs[1], }, start_chunk);
                        progress[idle_tid] = p;
                    }
#if defined(DEBUG)
                    cout << "[" << idle_tid << "] starts from " << progress[idle_tid]->toString() << endl;
#endif
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
#if defined(DEBUG)
                    cout << "Move " << idle_tid << " to worker_threads" << endl;
#endif
                    worker_threads.push_back(idle_tid);
                    continue;
                }
                idle_threads_iterator++;
            }
            start_parallel_sample = false;
            //assign the start sample to each thread
            if (!worker_threads.empty()) {
                goto INTERLEAVING_LOOP_0;
            } else {
                goto END_SAMPLE;
            }
        } // end of if (start_parallel_sample)
        while(true) {
            if (!idle_threads.empty() && dispatcher.hasNextChunk(1)) {
                auto idle_threads_iterator = idle_threads.begin();
                while(idle_threads_iterator != idle_threads.end() && dispatcher.hasNextChunk(1)) {
                    tid_to_run = *idle_threads_iterator;
                    if (!dispatcher.hasNextStaticChunk(tid_to_run)) {
                        idle_threads_iterator++;
                        continue;
                    }
                    Chunk c = dispatcher.getNextStaticChunk(tid_to_run);
                    vector<int> parallel_iteration_vector;
                    parallel_iteration_vector.push_back(c.first);
                    parallel_iteration_vector.push_back(0);
                    if (progress[tid_to_run]) {
                        progress[tid_to_run]->ref = "C0";
                        progress[tid_to_run]->iteration = parallel_iteration_vector;
                        progress[tid_to_run]->chunk = c;
                    } else {
                        Progress *p = new Progress("C0", parallel_iteration_vector, c);
                        progress[tid_to_run] = p;
                    }
                    idle_threads_iterator = idle_threads.erase(idle_threads_iterator);
                    worker_threads.emplace_back(tid_to_run);
#if defined(DEBUG)
                    cout << "Move " << tid_to_run << " to worker_threads" << endl;
#endif
                } /* end of progress assignment */
            } /* end of chunk availability check */
INTERLEAVING_LOOP_0:
            //UNIFORM INTERLEAVING
            working_threads = (int)worker_threads.size();
            sort(worker_threads.begin(), worker_threads.end());
            auto worker_thread_iterator = worker_threads.begin();
            while (worker_thread_iterator != worker_threads.end()) {
                tid_to_run = *worker_thread_iterator;
                if (!progress[tid_to_run]->isInBound()) {
                    worker_thread_iterator++;
#if defined(DEBUG)
                    cout << "[" << tid_to_run << "] Out of bound" << endl;;
#endif
                    continue;
                }
#if defined(DEBUG)
                if (start_reuse_search) {
                    cout << "[" << tid_to_run << "] Access " << progress[tid_to_run]->ref << " {";
                    for (unsigned i = 0; i < progress[tid_to_run]->iteration.size(); i++) {
                        cout << progress[tid_to_run]->iteration[i];
                        if (i < progress[tid_to_run]->iteration.size()-1) { cout << ","; }
                    }
                    cout << "} " << start_reuse_search << endl;
                }
#endif
                if (progress[tid_to_run]->ref == "C0") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C0(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samples are found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C1;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C1");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C0 */
                if (progress[tid_to_run]->ref == "C1") {
                    if (!start_reuse_search) { start_reuse_search= (start_tid == tid_to_run); }
                    if (start_reuse_search) {
                        Iteration *access = new Iteration("C1", {progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]});
                        string array = "C";
                        addr = GetAddress_C1(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
#if defined(DEBUG)
                        cout << access->toString() << " @ " << addr << endl;
#endif
                        if (access->compare(start_iteration) == 0) {
#if defined(DEBUG)
                            cout << "Meet the start sample " << access->toString() << endl;
#endif
                            isSample = true;
                        } else if ((!samples.empty() && access->compare(samples.top()) == 0) || (sample_names.find(access->toAddrString()) != sample_names.end())) {
#if defined(DEBUG)
                            cout << "Meet a new sample " << access->toString() << " while searching reuses" << endl;
#endif
                            traversing_samples++;
                            if (!samples.empty() && access->compare(samples.top()) == 0) {
                                samples.pop();
                            }
                            isSample = true;
                            samples_meet.insert(access->toAddrString());
                        }
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
#if defined(DEBUG)
                            cout << "[" << reuse << "] " << src->toString() << " -> " << access->toString() << endl;
#endif
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C1;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        if (isSample) {
                            LAT[tid_to_run][addr] = count[tid_to_run];
                            LATSampleIterMap[tid_to_run][count[tid_to_run]] = access;
                        } else { delete access; }
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    //CASE 2
                    progress[tid_to_run]->iteration.emplace_back(0);
                    progress[tid_to_run]->increment("A0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C1 */
                if (progress[tid_to_run]->ref == "A0") {
                    if (start_reuse_search) {
                        //A[c0][c2] will not access C
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("B0");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to A0 */
                if (progress[tid_to_run]->ref == "B0") {
                    if (start_reuse_search) {
                        //B[c2][c1] will not access C
                        count[tid_to_run] += 1;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C2");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to B0 */
                if (progress[tid_to_run]->ref == "C2") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C2(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C1;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    progress[tid_to_run]->increment("C3");
                    worker_thread_iterator++;
                    continue;
                } /* end of check to C2 */
                if (progress[tid_to_run]->ref == "C3") {
                    if (start_reuse_search) {
                        string array = "C";
                        addr = GetAddress_C3(progress[tid_to_run]->iteration[0],progress[tid_to_run]->iteration[1]);
                        bool isSample = false;
                        if (LAT[tid_to_run].find(addr) != LAT[tid_to_run].end()) {
                            long reuse = count[tid_to_run] - LAT[tid_to_run][addr];
                            //Src has Parallel Induction : 1
                            //Sink has Parallel Induction : 1
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
                            pluss_parallel_histogram_update(no_share_histogram,reuse,1.);
                            Iteration *src = LATSampleIterMap[tid_to_run][LAT[tid_to_run][addr]];
                            traversing_samples--;
                            //stop traversing if reuse of all samplesare found
#if defined(DEBUG)
                            if (samples.empty() && traversing_samples == 0) { cout << "[" << reuse << "] for last sample " << src->toString() << endl; };
#endif
                            if (samples.empty() && traversing_samples == 0) { goto END_SAMPLE; }
                            if (traversing_samples == 0) { 
#if defined(DEBUG)
                                cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                                delete src;
                                LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                                LAT[tid_to_run].erase(addr);
                                //Here we examine if there is an out-of-order effect.
                                //if the next sample we should jump has been traversed before, we will pop this sample directly.
                                //It is safe to call samples.top() once, since when entering here, 'samples' queue is not empty()
                                if (samples_meet.size() >= samples.size()) { goto END_SAMPLE; }
                                Iteration next = samples.top();
                                while(samples_meet.find(next.toAddrString()) != samples_meet.end()) {
#if defined(DEBUG)
                                    cout << "Skip " << next.toString() << " because we met this sample already " << endl;
#endif
                                    samples.pop();
                                    //All samples has been traversed, no need to jump
                                    if (samples.empty()) { break; }
                                    next = samples.top();
                                } // end of out-of-order check
                                if (!samples.empty()) {
#if defined(DEBUG)
                                    next = samples.top();
                                    cout << "Jump to next sample " << next.toString() << endl;
#endif
                                    goto START_SAMPLE_C1;
                                } else { goto END_SAMPLE; }
                            } // end of if (traversing_samples == 0)
#if defined(DEBUG)
                            cout << "delete sample " << src->toString() << ", active:" << traversing_samples << ", remain:" << samples.size() << endl;
#endif
                            delete src;
                            LATSampleIterMap[tid_to_run].erase(LAT[tid_to_run][addr]);
                            LAT[tid_to_run].erase(addr);
                        } // end of if (LAT.find(addr) != LAT.end())
                        count[tid_to_run]++;
                    } // end of if (start_reuse_search)
                    //CASE 3
                    if ((progress[tid_to_run]->iteration[2] + 1)<128) {
                        progress[tid_to_run]->iteration[2] = (progress[tid_to_run]->iteration[2] + 1);
                        progress[tid_to_run]->increment("A0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    if ((progress[tid_to_run]->iteration[1] + 1)<128) {
                        progress[tid_to_run]->iteration[1] = (progress[tid_to_run]->iteration[1] + 1);
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    //CASE 1
                    progress[tid_to_run]->iteration[0] = (progress[tid_to_run]->iteration[0] + 1);
                    if (progress[tid_to_run]->isInBound()) {
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.pop_back();
                        progress[tid_to_run]->iteration.emplace_back(0);
                        progress[tid_to_run]->increment("C0");
                        worker_thread_iterator++;
                        continue;
                    } /* end of check to C3 */
                    if (find(idle_threads.begin(), idle_threads.end(), tid_to_run) == idle_threads.end()) {
#if defined(DEBUG)
                        cout << "Move " << tid_to_run << " to idle threads" << endl;
#endif
                        idle_threads.emplace_back(tid_to_run);
                        worker_thread_iterator = worker_threads.erase(find(worker_threads.begin(), worker_threads.end(), tid_to_run));
                    }
                }
            } /* end of thread interleaving loop */
            if (idle_threads.size() == THREAD_NUM && !dispatcher.hasNextChunk(1)) {
                break;
            } /* end of break condition check */
        } /* end of while(true) */
        //reset both lists so they can be reused for later parallel loop
        idle_threads.clear();
        worker_threads.clear();
        for (unsigned i = 0; i < progress.size(); i++) {
            if (progress[i]) {
                delete progress[i];
                progress[i] = nullptr;
            }
        } // end of progress traversal
        goto END_SAMPLE;
    } // end of while (!samples.empty())
END_SAMPLE:
    if (!idle_threads.empty()) { idle_threads.clear(); };
    if (!worker_threads.empty()) { worker_threads.clear(); };
    if (!LAT.empty()) { ;
        for (unsigned i = 0; i < LAT.size(); i++) {
            pluss_parallel_histogram_update(no_share_histogram,-1,LAT[i].size());
            LAT.clear();;
        }
    }
    for (unsigned i = 0; i < progress.size(); i++) {
        if (progress[i]) {
            delete progress[i];
            progress[i] = nullptr;
        }
    }
    if (!LATSampleIterMap.empty()) {
        for (auto LATSampleIterMapEntry : LATSampleIterMap) {
            for (auto entry : LATSampleIterMapEntry.second) {
                delete entry.second;
            }
            LATSampleIterMapEntry.second.clear();
        }
        LATSampleIterMap.clear();
    }
    no_share_distribute(no_share_histogram,histogram);
    share_distribute(share_histogram,histogram);
    no_share_histogram.clear();
    share_histogram.clear();
    iteration_traversed_map["C1"] = accumulate(count.begin(), count.end(), 0);
    return;
} // end of the sampler function
int main() {
    long max_iteration_count = 0;
    unordered_map<long, double> histogram_C3;
    unordered_map<long, double> histogram_C2;
    unordered_map<long, double> histogram_A0;
    unordered_map<long, double> histogram_C0;
    unordered_map<long, double> histogram_B0;
    unordered_map<long, double> histogram_C1;
    pluss_timer_start();
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
    thread t_sampler_C3(sampler_C3,ref(histogram_C3));
    //set the thread affinity
    cpu_set_t cpuset_C3;
    CPU_ZERO(&cpuset_C3);
    CPU_SET(0,&cpuset_C3);
    pthread_setaffinity_np(t_sampler_C3.native_handle(),sizeof(cpu_set_t),&cpuset_C3);
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
    thread t_sampler_C2(sampler_C2,ref(histogram_C2));
    //set the thread affinity
    cpu_set_t cpuset_C2;
    CPU_ZERO(&cpuset_C2);
    CPU_SET(2,&cpuset_C2);
    pthread_setaffinity_np(t_sampler_C2.native_handle(),sizeof(cpu_set_t),&cpuset_C2);
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
    thread t_sampler_A0(sampler_A0,ref(histogram_A0));
    //set the thread affinity
    cpu_set_t cpuset_A0;
    CPU_ZERO(&cpuset_A0);
    CPU_SET(4,&cpuset_A0);
    pthread_setaffinity_np(t_sampler_A0.native_handle(),sizeof(cpu_set_t),&cpuset_A0);
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
    thread t_sampler_C0(sampler_C0,ref(histogram_C0));
    //set the thread affinity
    cpu_set_t cpuset_C0;
    CPU_ZERO(&cpuset_C0);
    CPU_SET(6,&cpuset_C0); // TODO: do not understand why we space the threads out like (4, 6, 8, 10)
    pthread_setaffinity_np(t_sampler_C0.native_handle(),sizeof(cpu_set_t),&cpuset_C0);
/* (c2,0,\<,128,(c2 + 1)), 0*/
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
    thread t_sampler_B0(sampler_B0,ref(histogram_B0));
    //set the thread affinity
    cpu_set_t cpuset_B0;
    CPU_ZERO(&cpuset_B0);
    CPU_SET(8,&cpuset_B0);
    pthread_setaffinity_np(t_sampler_B0.native_handle(),sizeof(cpu_set_t),&cpuset_B0);
/* (c1,0,\<,128,(c1 + 1)), 0*/
/* (c0,0,\<,128,(c0 + 1)), 1*/
    thread t_sampler_C1(sampler_C1,ref(histogram_C1));
    //set the thread affinity
    cpu_set_t cpuset_C1;
    CPU_ZERO(&cpuset_C1);
    CPU_SET(10,&cpuset_C1);
    pthread_setaffinity_np(t_sampler_C1.native_handle(),sizeof(cpu_set_t),&cpuset_C1);
    t_sampler_C3.join();
    t_sampler_C2.join();
    t_sampler_A0.join();
    t_sampler_C0.join();
    t_sampler_B0.join();
    t_sampler_C1.join();
    //Merge the histogram found by each reference
    for (auto pair : histogram_C3) {
        pluss_histogram_update(pair.first,pair.second);
    }
    for (auto pair : histogram_C2) {
        pluss_histogram_update(pair.first,pair.second);
    }
    for (auto pair : histogram_A0) {
        pluss_histogram_update(pair.first,pair.second);
    }
    for (auto pair : histogram_C0) {
        pluss_histogram_update(pair.first,pair.second);
    }
    for (auto pair : histogram_B0) {
        pluss_histogram_update(pair.first,pair.second);
    }
    for (auto pair : histogram_C1) {
        pluss_histogram_update(pair.first,pair.second);
    }
    pluss_AET();
    pluss_timer_stop();
    pluss_timer_print();
    _pluss_histogram_print("C3", histogram_C3);
    _pluss_histogram_print("C2", histogram_C2);
    _pluss_histogram_print("A0", histogram_A0);
    _pluss_histogram_print("C0", histogram_C0);
    _pluss_histogram_print("B0", histogram_B0);
    _pluss_histogram_print("C1", histogram_C1);

    pluss_print_histogram();
    pluss_print_mrc();
    for (auto entry: iteration_traversed_map) {
        max_iteration_count = max(max_iteration_count, entry.second);
    }
    cout << "max iteration traversed" << endl;
    cout << max_iteration_count << endl;
    return 0;
}
