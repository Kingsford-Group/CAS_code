#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <random>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <iterator>
#include <algorithm>

using namespace std;

struct ETS_entry {
    seqan::String<seqan::Dna5> ets_seed;
    uint8_t *distance;
};

struct SP_entry{
    seqan::String<seqan::Dna5> sign_post;
    vector<unsigned> pop_count;
    uint8_t biggest_dist;
};

#define THREAD_ALO 200

ETS_entry *entry_ary;
ETS_entry **merge_workspace;
ETS_entry **swapping_space;
unsigned merge_entry_num;
vector<unsigned> *entry_off_ary;
SP_entry *sp_ary;
unsigned * sp_pop_count;

unsigned sp_size;
unsigned sp_idx;
unsigned entry_size;
unsigned total_threads;

seqan::AlignConfig<false, false, true, true> alignConfig;
seqan::Score<int, seqan::Simple> scoreConfig(0, -1, -1);

struct ETS_cmp {
    bool operator() (const int& a, const int& b) {
        for (unsigned i = 0; i <= sp_idx; ++i) {
            if (entry_ary[a].distance[i] > entry_ary[b].distance[i])
                return false;
            else if (entry_ary[a].distance[i] < entry_ary[b].distance[i])
                return true;
        }

        return true;
    }
} ETS_cmp_obj;

int ETS_compare (const void * a, const void * b) {
    ETS_entry * a_cast = (ETS_entry*) a;
    ETS_entry * b_cast = (ETS_entry*) b;
    
    for (unsigned i = 0; i <= sp_idx; ++i) {
        if (a_cast->distance[i] > b_cast->distance[i])
            return 1;
        else if (a_cast->distance[i] < b_cast->distance[i])
            return -1;
    }

    return 0;
}

// Returns the pointer to the first entry that is equal or greater than target_distance, 
// Returns -1 when further search is guaranteed to be ver
int binarySearch (unsigned int dist_idx, int head, int tail, unsigned target_distance) {
    // Insure invariance: head < target and tail >= target
    if (entry_ary[entry_off_ary->at(head)].distance[dist_idx] >= target_distance)
        return head;
    else if (entry_ary[entry_off_ary->at(tail)].distance[dist_idx] < target_distance)
        return tail + 1;

    unsigned long length = tail - head;

    while (length > 1) {
        unsigned long mid = ((unsigned long) tail + head) / 2;

        if (entry_ary[entry_off_ary->at(mid)].distance[dist_idx] < target_distance)
            head = mid;
        else
            tail = mid;
    
        length = tail - head;
    }

    //printEntry(ets_spdb[head]);
    //printEntry(ets_spdb[tail]);

    return tail;
}

int find_entry_sp_idx_helper(int head, int tail, unsigned int depth_idx, int compensation_range) {
    //cout << "Entering at -- depth_idx: " << depth_idx << " head: " << head << " tail: " << tail << endl;
    // return if the depth is too much
    if ((unsigned) depth_idx > sp_idx)
        return head;

    int comp_idx = -compensation_range;
    int result_entry_idx = 0;
    int new_head, new_tail;

    do {
        int target_dist = sp_ary[depth_idx].biggest_dist + comp_idx;
        
        if (target_dist < 0) {
            target_dist = 0;
            comp_idx = -sp_ary[depth_idx].biggest_dist;
        }

        new_head = binarySearch(depth_idx, head, tail, target_dist);
        new_tail = binarySearch(depth_idx, head, tail, target_dist + 1) - 1;

        if (new_head > new_tail) {
            //cout << "Dead end at: comp_idx: " << comp_idx << " target_dist: " << target_dist  << " depth_idx: " << depth_idx << " head: " << head << " tail: " << tail << endl;
            ++comp_idx;
            continue;
        }

        result_entry_idx = find_entry_sp_idx_helper(new_head, new_tail, depth_idx + 1, compensation_range);

        if (result_entry_idx != -1)
            return result_entry_idx;

        ++comp_idx;
    }
    while (comp_idx <= compensation_range);

    return -1;
}

int find_entry_sp_idx() {
    int head = 0;
    int tail = entry_size - 1;

    /*
    for (unsigned i = 0; i <= sp_idx; ++i) {
        head = binarySearch(i, head, tail, sp_ary[i].biggest_dist);
        tail = binarySearch(i, head, tail, sp_ary[i].biggest_dist + 1) - 1;

        assert(head > 0 && tail > 0);
    }
    */
    int comp_range = 0;
  
    while (true) {
        cout << "comp_range: " << comp_range << endl;

        int result = find_entry_sp_idx_helper(head, tail, 0, comp_range);
        cout << "result: " << result << endl;

        if (result >= 0)
            return result;

        ++comp_range;
    }
}

void ETS_merge(unsigned beg_pos, unsigned mid_pos, unsigned end_pos) {
    unsigned t_id = omp_get_thread_num();

    // A and B are the two halfs of the bigger to-be-merged array
    unsigned total_a = mid_pos - beg_pos;
    unsigned total_b = end_pos - mid_pos;

    unsigned progress_a = 0;
    unsigned progress_b = 0;

    unsigned scan_pos = beg_pos;

    while (progress_a < total_a && progress_b < total_b) {
        unsigned scan_count = 0;
        unsigned ori_scan_pos = scan_pos;
        unsigned count_b = 0;
        unsigned count_a = 0;
        unsigned last_b = 0;

        while (scan_count < merge_entry_num && scan_pos < end_pos) {

            if (entry_off_ary->at(scan_pos) >= mid_pos) {
                ++count_b;
                last_b = scan_count;
            }
            ++scan_count;
            ++scan_pos;
        }

        count_a = scan_count - count_b;

        if (count_b != 0) {
            unsigned b_idx = 0;
            unsigned merge_idx = ori_scan_pos;
            while (b_idx != count_b) {
                unsigned ws_offset = merge_idx - ori_scan_pos;
                unsigned entry_idx_offset = entry_off_ary->at(merge_idx);

                // from b
                if (entry_idx_offset >= mid_pos) {
                    ++b_idx;
                }
                // from a
                else {
                    entry_idx_offset += progress_b;
                }
                
                memcpy(&(merge_workspace[t_id][ws_offset]), &entry_ary[entry_idx_offset], sizeof(ETS_entry));

                ++merge_idx;
            }

            ++last_b;
            unsigned skip_a = last_b - count_b;

            unsigned total_swap = total_a - progress_a - skip_a;
            // move a down
            while (total_swap > 0) {
                unsigned swap_size = (total_swap < merge_entry_num) ? total_swap : merge_entry_num;
                total_swap -= swap_size;
                unsigned src_off = ori_scan_pos + skip_a + total_swap;
                unsigned dest_off = ori_scan_pos + last_b + total_swap;

                memcpy(swapping_space[t_id], &entry_ary[src_off], swap_size * sizeof(ETS_entry));
                memcpy(&entry_ary[dest_off], swapping_space[t_id], swap_size * sizeof(ETS_entry));
                //memcpy(&entry_ary[ori_scan_pos + last_b], &entry_ary[ori_scan_pos + skip_a], (total_a - progress_a - skip_a) * sizeof(ETS_entry) );
            }
            memcpy(&entry_ary[ori_scan_pos], merge_workspace[t_id], last_b * sizeof(ETS_entry) );
        }

        progress_a += count_a;
        progress_b += count_b;
    }
}

#define SRT_BATCH 10000000

void parallel_sort() {
    unsigned total_batches = (entry_size - 1) / SRT_BATCH + 1;

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < total_batches; ++i) {
        auto iter_begin = entry_off_ary->begin() + i * SRT_BATCH;
        auto iter_end = (i == total_batches - 1) ? entry_off_ary->end() : entry_off_ary->begin() + i * SRT_BATCH + SRT_BATCH;

        unsigned ary_begin = i * SRT_BATCH;
        unsigned ary_size = (i == total_batches - 1) ? entry_size - ary_begin : SRT_BATCH;

        cout << "begin sorting " << i << endl;

        //sort(iter_begin, iter_end, ETS_cmp_obj);
        qsort(entry_ary + ary_begin, ary_size, sizeof(ETS_entry), ETS_compare);

        cout << "finished sorting " << i << endl;
    }

    #pragma omp parallel for schedule(dynamic)
    for (unsigned x = 0; x < entry_size; ++x) {
        entry_off_ary->at(x) = x;
    }

    cout << "finished parallel sort" << endl;

    unsigned long merge_size = SRT_BATCH * 2;

    while (true) {
       total_batches = (entry_size - 1) / merge_size + 1;
       
        #pragma omp parallel for schedule(dynamic)
        for (unsigned i = 0; i < total_batches; ++i) {
            //cout << "starting " << i << endl;
            if (i != total_batches - 1 || entry_size > i * merge_size + merge_size / 2) {
                auto beg_iter = entry_off_ary->begin() + i * merge_size;
                auto mid_iter = beg_iter + merge_size / 2; auto end_iter = (i == total_batches - 1) ? entry_off_ary->end() : beg_iter + merge_size; //cout << "running " << i << endl;
                inplace_merge(beg_iter, mid_iter, end_iter, ETS_cmp_obj);

                unsigned beg_pos = i * merge_size;
                unsigned mid_pos = beg_pos + merge_size / 2;
                unsigned end_pos = (i == total_batches - 1) ? entry_size : beg_pos + merge_size;

                if (ETS_compare(&entry_ary[entry_off_ary->at(mid_pos - 1)], &entry_ary[entry_off_ary->at(mid_pos)]) <= 0)
                    ETS_merge(beg_pos, mid_pos, end_pos);
            }
        }

        cout << "finished parallel merge of size: " << merge_size << endl;

        #pragma omp parallel for schedule(dynamic)
        for (unsigned x = 0; x < entry_size; ++x) {
            entry_off_ary->at(x) = x;
        }

        //if (merge_size >= SRT_BATCH * 32)
        if (merge_size > entry_size)
            break;

        merge_size *= 2;
    }
}

#define TRD_BATCH 10000000

unsigned single_iteration (unsigned entry_sp_idx) {

    sp_ary[sp_idx].sign_post = entry_ary[entry_off_ary->at(entry_sp_idx)].ets_seed;

    for (unsigned i = 0; i < total_threads * THREAD_ALO; ++i)
        sp_pop_count[i] = 0;

    #pragma omp parallel for schedule(dynamic, TRD_BATCH)
    for (unsigned i = 0; i < entry_size; ++i) {

        int dist = -seqan::globalAlignmentScore(entry_ary[i].ets_seed, sp_ary[sp_idx].sign_post, scoreConfig, alignConfig);

        if (dist > 63)
            dist = 63;

        entry_ary[i].distance[sp_idx] = dist;
        unsigned t_id = omp_get_thread_num();

        ++sp_pop_count[t_id * THREAD_ALO + dist];

        //seqan::String<char> temp0 = road_signs[omp_get_thread_num()][omp_get_thread_num()];
        //string temp = seqan::toCString(temp0);
        //str_ary[i] += temp;
        //str_ary[i] += to_string(omp_get_thread_num() );
    }

    //qsort((*entry_ary), entry_size, sizeof(ETS_entry), ETS_compare);
    cout << "started sorting" << endl;
    parallel_sort();

    sp_ary[sp_idx].pop_count.resize(THREAD_ALO);

    for (unsigned i = 0; i < THREAD_ALO; ++i)
        sp_ary[sp_idx].pop_count[i] = 0;

    for (unsigned i = 0; i < total_threads; ++i) {
        for (unsigned j = 0; j < THREAD_ALO; ++j)
            sp_ary[sp_idx].pop_count[j] += sp_pop_count[i * THREAD_ALO + j];
    }
    
    unsigned max_sp_idx = 0;
    unsigned max_sp_count = 0;

    for (unsigned i = 0; i < THREAD_ALO; ++i) {
        cout << "pop_idx: " << i << endl;
        cout << "pop_count: " << sp_ary[sp_idx].pop_count[i] << endl;

        if (sp_ary[sp_idx].pop_count[i] > max_sp_count) {
            max_sp_idx = i;
            max_sp_count = sp_ary[sp_idx].pop_count[i];
        }
    }

    cout << "max_sp_idx: " << max_sp_idx << endl;

    sp_ary[sp_idx].biggest_dist = max_sp_idx;

    entry_sp_idx = find_entry_sp_idx();

    cout << "Done with sp: " << sp_idx << endl;
    cout << "*entry_off_ary[" << entry_sp_idx << "]: " << (unsigned) entry_ary[entry_off_ary->at(entry_sp_idx)].distance[sp_idx] << endl;

    return entry_sp_idx;
}

void generate_sp(char* argv[]) {

    ofstream sp_file;
    ofstream output_file;

    sp_file.open(argv[4]);
    output_file.open(argv[5]);

    unsigned entry_sp_idx;
   
    if (sp_idx == 0)
        entry_sp_idx = entry_size / 2;
    else
        entry_sp_idx = find_entry_sp_idx();

    cout << "starting with entry_sp_idx: " << entry_sp_idx << endl;
    //++sp_idx;

    for (; sp_idx < sp_size; ++sp_idx) {

        cout << "running sp: " << sp_idx;

        entry_sp_idx = single_iteration(entry_sp_idx);

/*        
        sp_ary[sp_idx].sign_post = (*entry_ary)[entry_sp_idx].ets_seed;

        for (unsigned i = 0; i < total_threads * THREAD_ALO; ++i)
            sp_pop_count[i] = 0;

        #pragma omp parallel for schedule(dynamic, TRD_BATCH)
        for (unsigned i = 0; i < entry_size; ++i) {

            int dist = -seqan::globalAlignmentScore((*entry_ary)[i].ets_seed, sp_ary[sp_idx].sign_post, scoreConfig, alignConfig);

            (*entry_ary)[i].distance[sp_idx] = dist;
            unsigned t_id = omp_get_thread_num();

            ++sp_pop_count[t_id * THREAD_ALO + dist];

            //seqan::String<char> temp0 = road_signs[omp_get_thread_num()][omp_get_thread_num()];
            //string temp = seqan::toCString(temp0);
            //str_ary[i] += temp;
            //str_ary[i] += to_string(omp_get_thread_num() );
        }

        qsort((*entry_ary), entry_size, sizeof(ETS_entry), ETS_compare);

        sp_ary[sp_idx].pop_count.resize(THREAD_ALO);

        for (unsigned i = 0; i < THREAD_ALO; ++i)
            sp_ary[sp_idx].pop_count[i] = 0;

        for (unsigned i = 0; i < total_threads; ++i) {
            for (unsigned j = 0; j < THREAD_ALO; ++j)
                sp_ary[sp_idx].pop_count[j] += sp_pop_count[i * THREAD_ALO + j];
        }
        
        unsigned max_sp_idx = 0;
        unsigned max_sp_count = 0;

        for (unsigned i = 0; i < THREAD_ALO; ++i) {
            cout << "pop_idx: " << i << endl;
            cout << "pop_count: " << sp_ary[sp_idx].pop_count[i] << endl;

            if (sp_ary[sp_idx].pop_count[i] > max_sp_count) {
                max_sp_idx = i;
                max_sp_count = sp_ary[sp_idx].pop_count[i];
            }
        }

        cout << "max_sp_idx: " << max_sp_idx << endl;

        sp_ary[sp_idx].biggest_dist = max_sp_idx;

        entry_sp_idx = find_entry_sp_idx();

        cout << "Done with sp: " << sp_idx << endl;
        cout << "sp_entry[" << entry_sp_idx << "]: " << (*entry_ary)[entry_sp_idx].distance[sp_idx] << endl;

        */
    }

    for (unsigned i = 0; i < entry_size; ++i) {
        output_file << entry_ary[i].ets_seed << " ";

        for (unsigned j = 0; j < sp_size; ++j)
            output_file << (unsigned) entry_ary[i].distance[j] << " ";

        output_file << endl;
    }

    for (unsigned i = 0; i < sp_size; ++i)
        sp_file << sp_ary[i].sign_post << " " << (unsigned) sp_ary[i].biggest_dist << endl;

    sp_file.close();
    output_file.close();
}


int main (int argc, char *argv[]) {
    merge_entry_num = 1000000;

    if (argc != 6 && argc != 8) {
        cerr << "Not enough arguments."
             << "./sp_gen input_file input_size sp_size sp_file output_file: " << endl;
        cerr << "./sp_gen input_file input_size sp_out_size sp_out_file output_file sp_in_file sp_in_size: " << endl;
        return -1;
    }

    entry_size = atoi(argv[2]);
    sp_size = atoi(argv[3]);

    total_threads = omp_get_max_threads();
    cout << "Available threads: " << total_threads << endl;

    entry_ary = new ETS_entry[entry_size];
    entry_off_ary = new vector<unsigned>;
    entry_off_ary->resize(entry_size);
    sp_ary = new SP_entry [sp_size];

    sp_pop_count = new unsigned [total_threads * THREAD_ALO];

    for (unsigned i = 0; i < entry_size; ++i) {
        entry_ary[i].distance = new uint8_t [sp_size];
    }

    merge_workspace = new ETS_entry* [total_threads];
    swapping_space = new ETS_entry* [total_threads];

    for (unsigned i = 0; i < total_threads; ++i) {
        merge_workspace[i] = new ETS_entry[merge_entry_num];
        swapping_space[i] = new ETS_entry[merge_entry_num];
    }

    // In debug mode!
    if (argc == 8) {
        string seed;
        ifstream input_file;
        input_file.open(argv[1]);

        sp_idx = atoi(argv[7]);
        unsigned temp_num;

        for (unsigned i = 0; i < entry_size; ++i) {
            input_file >> seed;
            entry_ary[i].ets_seed = seed;
            entry_off_ary->at(i) = i;

            for (unsigned j = 0; j < sp_idx; ++j) {
                input_file >> temp_num;
                entry_ary[i].distance[j] = temp_num;
            }
        }

        input_file.close();

        input_file.open(argv[6]);

        for (unsigned i = 0; i < sp_idx; ++i) {
            input_file >> seed;
            sp_ary[i].sign_post = seed;
            input_file >> temp_num;
            sp_ary[i].biggest_dist = temp_num;
        }

        input_file.close();

        //++sp_idx;

        cout << "finished loading" << endl;

        generate_sp(argv);
    }
    else {
        ifstream input_file;
        input_file.open(argv[1]);
        string seed;

        for (unsigned i = 0; i < entry_size; ++i) {
            input_file >> seed;

            entry_ary[i].ets_seed = seed;
        }

        input_file.close();
        sp_idx = 0;
        generate_sp(argv);
    }

    return 0;
}

