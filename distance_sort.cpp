#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <queue>
#include <fstream>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/align.h>

using namespace std;

struct ETS_entry {
    string ets_seed;
    unsigned *distance;
};

unsigned sp_size = 0;
ETS_entry **unsorted_ary;
ETS_entry *sorted_ary;
unsigned *progress_ary;
unsigned *length_ary;
ifstream *file_ary;
ETS_entry *file_ets_ary;

#define STR_BATCH 100000000
#define TRD_BATCH 1000000
#define BATCH_NUM STR_BATCH / TRD_BATCH

int ETS_compare (const void * a, const void * b) {
    ETS_entry * a_cast = (ETS_entry*) a;
    ETS_entry * b_cast = (ETS_entry*) b;
    
    for (unsigned i = 0; i < sp_size; ++i) {
        if (a_cast->distance[i] > b_cast->distance[i])
            return 1;
        else if (a_cast->distance[i] < b_cast->distance[i])
            return -1;
    }

    return 0;
}


int main (int argc, char *argv[]) {

    if (argc != 4) {
        cerr << "Not enough arguments. ./disaster unsorted_file sp_file sorted_file" << endl;
        return -1;
    }

    unsigned total_threads = omp_get_max_threads();
    cout << "Available threads: " << total_threads << endl;

    ifstream input_file;
    ifstream sp_file;

    input_file.open(argv[1]);
    sp_file.open(argv[2]);
    string tmp_dir = argv[3];

    if (tmp_dir[tmp_dir.length() - 1] != '/')
        tmp_dir += '/';

    string sp;

    sp_file >> sp;

    while (!sp_file.eof() ) {
        ++sp_size;
        sp_file >> sp;
    }

    unsorted_ary = new ETS_entry* [BATCH_NUM];

    for (unsigned i = 0; i < BATCH_NUM; ++i) {
        unsorted_ary[i] = new ETS_entry[TRD_BATCH];
    }

    sorted_ary = new ETS_entry [STR_BATCH];

    for (unsigned i = 0; i < BATCH_NUM; ++i) {
        for (unsigned j = 0; j < TRD_BATCH; ++j)
            unsorted_ary[i][j].distance = new unsigned [sp_size];
    }

    for (unsigned i = 0; i < STR_BATCH; ++i)
        sorted_ary[i].distance = new unsigned [sp_size];

    progress_ary = new unsigned [BATCH_NUM];
    length_ary = new unsigned [BATCH_NUM];

    bool still_reading = true;
    unsigned str_num;
    unsigned batch_num = 0;

    cout << "sorting starts" << endl;

    auto cmp = [](int left, int right) { return ETS_compare(&unsorted_ary[left][progress_ary[left]], &unsorted_ary[right][progress_ary[right]] ) > 0; };
    priority_queue<unsigned, vector<unsigned>, decltype(cmp)> batch_heap(cmp);

    while (still_reading) {

        str_num = 0;

        input_file >> sorted_ary[str_num].ets_seed;

        for (unsigned i = 0; i < sp_size; ++i)
            input_file >> sorted_ary[str_num].distance[i];
        
        still_reading = !input_file.eof();

        if (still_reading)
            ++str_num;

        while (still_reading && str_num < STR_BATCH) {
         
            input_file >> sorted_ary[str_num].ets_seed;

            for (unsigned i = 0; i < sp_size; ++i)
                input_file >> sorted_ary[str_num].distance[i];

            still_reading = !input_file.eof();

            if (still_reading)
                ++str_num;
        }

        unsigned total_batches = str_num / TRD_BATCH;
   
        #pragma omp parallel for schedule(dynamic, 1)
        for (unsigned i = 0; i < total_batches; ++i) {
            unsigned ary_size = TRD_BATCH;

            if (i == total_batches - 1)
                ary_size = str_num - i * TRD_BATCH;

            memcpy(unsorted_ary[i], sorted_ary + (i * TRD_BATCH),  ary_size * sizeof (ETS_entry) );

            qsort(unsorted_ary[i], TRD_BATCH, sizeof(ETS_entry), ETS_compare);
        }


        for (unsigned i = 0; i < total_batches; ++i) {
            unsigned ary_size = TRD_BATCH;

            if (i == total_batches - 1)
                ary_size = str_num - i * TRD_BATCH;

            progress_ary[i] = 0;
            length_ary[i] = ary_size;

            batch_heap.push(i);
        }

        unsigned sorted_idx = 0;
        
        while (!batch_heap.empty() ) {
            unsigned batch_id = batch_heap.top();
            batch_heap.pop();

            unsigned tr_id = progress_ary[batch_id]++;

            if (progress_ary[batch_id] < length_ary[batch_id])
                batch_heap.push(batch_id);

            memcpy(sorted_ary + sorted_idx, unsorted_ary[batch_id] + tr_id, sizeof (ETS_entry)  );

            sorted_idx++;
        }

        /*
        for (unsigned i = 0; i < total_batches; ++i) {
            unsigned ary_size = TRD_BATCH;

            if (i == total_batches - 1)
                ary_size = str_num - i * TRD_BATCH;

            memcpy(sorted_ary + (i * TRD_BATCH), unsorted_ary[i], ary_size * sizeof (ETS_entry)  );
        }
        */


        ofstream output_file;
        output_file.open(tmp_dir + "__" + to_string(batch_num) );
        
        for (unsigned i = 0; i < str_num; ++i) {
            output_file << sorted_ary[i].ets_seed << " ";
            for (unsigned j = 0; j < sp_size; ++j)
                output_file << sorted_ary[i].distance[j] << " ";
            output_file << endl;
        }

        output_file.close();

        cout << "finished batch: " << batch_num << endl;
        ++batch_num;
        //break;
    }

    input_file.close();
    sp_file.close();

    file_ary = new ifstream[batch_num];
    file_ets_ary = new ETS_entry[batch_num];

    auto cmp_file = [](int left, int right) { return ETS_compare(&file_ets_ary[left], &file_ets_ary[right] ) > 0; };
    priority_queue<unsigned, vector<unsigned>, decltype(cmp_file)> file_heap(cmp_file);

    for (unsigned i = 0; i < batch_num; ++i) {
        file_ary[i].open(tmp_dir + "__" + to_string(i) );

        file_ary[i] >> file_ets_ary[i].ets_seed;

        file_ets_ary[i].distance = new unsigned [sp_size];

        for (unsigned j = 0; j < sp_size; ++j) {
            file_ary[i] >> file_ets_ary[i].distance[j];
        }

        file_heap.push(i);
    }

    ofstream final_output_file;
    final_output_file.open((string)argv[1] + "_fo");

    while (!file_heap.empty() ) {
        unsigned file_id = file_heap.top();
        file_heap.pop();

        final_output_file << file_ets_ary[file_id].ets_seed << " ";

        for (unsigned i = 0; i < sp_size; ++i)
            final_output_file << file_ets_ary[file_id].distance[i] << " ";

        final_output_file << endl;

        file_ary[file_id] >> file_ets_ary[file_id].ets_seed;

        for (unsigned i = 0; i < sp_size; ++i) {
            file_ary[file_id] >> file_ets_ary[file_id].distance[i];
        }

        if (!file_ary[file_id].eof() )
            file_heap.push(file_id);
    }

    for (unsigned i = 0; i < batch_num; ++i) {
        file_ary[i].close();

        remove((tmp_dir + "__" + to_string(i)).c_str() );
    }

    final_output_file.close();

    return 0;
}

