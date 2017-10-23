#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <stdexcept>

using namespace std;

seqan::AlignConfig<false, false, true, true> alignConfig;
seqan::Score<int, seqan::Simple> scoreConfig(0, -1, -1);

struct ETS_entry {
    string ets_seed;
    unsigned *distance;
};

struct ETS_result {
    string ets_seed;
    unsigned cloest_distance;
};

unsigned sp_size = 0;
unsigned spdb_size = 0;
ETS_entry* ets_spdb;
ETS_entry* ets_task;
ETS_result* ets_result;
unsigned max_dist = 8;

#define STR_BATCH 1000
#define TRD_BATCH 10

void printEntry(const ETS_entry& entry) {
    cout << entry.ets_seed << " ";
    for (unsigned i = 0; i < sp_size; ++i)
        cout << entry.distance[i] << " ";

    cout << endl;
}

// Returns the pointer to the first entry that is equal or greater than target_distance, 
// Returns -1 when further search is guaranteed to be ver
int binarySearch (unsigned int dist_idx, int head, int tail, unsigned target_distance) {
    // Insure invariance: head < target and tail >= target
    if (ets_spdb[head].distance[dist_idx] >= target_distance)
        return head;
    else if ( ets_spdb[tail].distance[dist_idx] < target_distance)
        return -1;

    int length = tail - head;

    while (length > 1) {
        int mid = (tail + head) / 2;

        if (ets_spdb[mid].distance[dist_idx] < target_distance)
            head = mid;
        else
            tail = mid;
    
        length = tail - head;
    }

    //printEntry(ets_spdb[head]);
    //printEntry(ets_spdb[tail]);

    return tail;
}

int closestDistanceHelper(const ETS_entry &entry, unsigned dist_idx, int head, int tail, int distance) {
    //cout << "dist_idx: " << dist_idx << " head: " << head << " tail: " << tail << endl;
    
    // compute all the distances
    if (dist_idx == sp_size) {
        seqan::String<seqan::Dna5> pattern = entry.ets_seed;

        for (int i = head; i <= tail; ++i) {

            if (ets_spdb[i].ets_seed.find('N') != string::npos)
                continue;

            seqan::String<seqan::Dna5> db_pattern = ets_spdb[i].ets_seed;

            int score = -seqan::globalAlignmentScore(pattern, db_pattern, scoreConfig, alignConfig);

            int substr_len = entry.ets_seed.length() - score;

            if (score < distance && entry.ets_seed.substr(0, substr_len) != ets_spdb[i].ets_seed.substr(score, substr_len) && entry.ets_seed.substr(score, substr_len) != ets_spdb[i].ets_seed.substr(0, substr_len)) {
                //printEntry(ets_spdb[i]);

                if (score < distance)
                    distance = score;

                //cout << score << endl;
            }

        }
    }
    // go down one more level
    else {
        int begin_dist = entry.distance[dist_idx] - distance;
        int end_dist = entry.distance[dist_idx] + distance;

        if (begin_dist < 0)
            begin_dist = 0;

        int next_head = -1;
        int next_tail = -1;

        next_head = binarySearch(dist_idx, head, tail, begin_dist);
        
        int dist;
       
        if (next_head != -1) 
            dist = ets_spdb[next_head].distance[dist_idx];

        while (next_head != -1 && dist <= end_dist) {
            next_tail = binarySearch(dist_idx, head, tail, dist + 1);

            int corrected_tail;
            if (next_tail == -1)
                corrected_tail = tail;
            else
                corrected_tail = next_tail - 1;
            
            //cout << "dist: " << dist << " next_head: " << next_head << " next_tail: " << next_tail << endl;
            // tail - 1 generates problem here.
            distance = closestDistanceHelper(entry, dist_idx + 1, next_head, corrected_tail, distance);

            next_head = next_tail;
            if (next_head != -1)
                dist = ets_spdb[next_head].distance[dist_idx];
        }
    }

    return distance;
}

int closestDistance(const ETS_entry &entry) {
    int distance = max_dist;

    unsigned dist_idx = 0;

    int head = 0;
    int tail = spdb_size - 1;

    distance = closestDistanceHelper(entry, dist_idx, head, tail, distance);

    /*
    printEntry(entry);
    cout << distance << endl;
    cout.flush();
    */

    return distance;
}

int main (int argc, char *argv[]) {
    unsigned total_threads = omp_get_max_threads();
    cout << "Available threads: " << total_threads << endl;

    if (argc != 7) {
        cerr << "Not enough arguments. ./ets_db_construct_relief input_file sp_db_file sp_file sp_db_size max_range_distance output_file" << endl;
        return -1;
    }

    ifstream input_file;
    ifstream spdb_file;
    ifstream sp_file;
    ofstream output_file;

    input_file.open(argv[1]);
    spdb_file.open(argv[2]);
    sp_file.open(argv[3]);
    spdb_size = atoi(argv[4]);
    max_dist = atoi(argv[5]);
    output_file.open(argv[6]);

    string sp;

    sp_file >> sp;

    while (!sp_file.eof() ) {
        ++sp_size;
        sp_file >> sp;
    }

    sp_file.close();

    ets_spdb = new ETS_entry[spdb_size];

    for (unsigned i = 0; i < spdb_size; ++i) {
        ets_spdb[i].distance = new unsigned [sp_size];
        spdb_file >> ets_spdb[i].ets_seed;

        for (unsigned j = 0; j < sp_size; ++j) {
            spdb_file >> ets_spdb[i].distance[j];
        }
    }

    cout << "the full database is loaded" << endl;

    ets_task = new ETS_entry[STR_BATCH];

    for (unsigned i = 0; i < STR_BATCH; ++i)
        ets_task[i].distance = new unsigned[sp_size];

    ets_result = new ETS_result[STR_BATCH];

    bool still_reading = true;
    unsigned ets_num;
    unsigned batch_num = 0;

    while (still_reading) {

        ets_num = 0;

        input_file >> ets_task[ets_num].ets_seed;

        for (unsigned i = 0; i < sp_size; ++i)
            input_file >> ets_task[ets_num].distance[i];
        
        still_reading = !input_file.eof();

        if (still_reading)
            ++ets_num;

        while (still_reading && ets_num < STR_BATCH) {
         
            input_file >> ets_task[ets_num].ets_seed;

            for (unsigned i = 0; i < sp_size; ++i)
                input_file >> ets_task[ets_num].distance[i];

            still_reading = !input_file.eof();

            if (still_reading && ets_task[ets_num].ets_seed.find('N') == string::npos)
                ++ets_num;
        }
   
        #pragma omp parallel for schedule(dynamic, TRD_BATCH)
        for (unsigned i = 0; i < ets_num; ++i) {
            ets_result[i].ets_seed = ets_task[i].ets_seed;

            //seqan::String<seqan::Dna5> pattern = str_ary[i];

            ets_result[i].cloest_distance = closestDistance(ets_task[i]);

            //cout << "i: " << i << endl;

            //seqan::String<char> temp0 = road_signs[omp_get_thread_num()][omp_get_thread_num()];
            //string temp = seqan::toCString(temp0);
            //str_ary[i] += temp;
            //str_ary[i] += to_string(omp_get_thread_num() );
        }


        for (unsigned i = 0; i < ets_num; ++i) {
            output_file << ets_result[i].ets_seed << " " << ets_result[i].cloest_distance << endl;
        }

        cout << "finished batch: " << batch_num << endl;
        ++batch_num;
    }

    input_file.close();
    sp_file.close();
    output_file.close();

    return 0;
}

