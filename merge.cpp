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
#include <sys/stat.h>
#include <dirent.h>

using namespace std;

struct ETS_entry {
    string ets_seed;
    unsigned *distance;
};

unsigned sp_size = 0;
unsigned spdb_size = 0;
ETS_entry* ets_spdb;

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

    if (argc != 5) {
        cerr << "Not enough arguments. ./disaster_relief input_file sp_db_file sp_file sp_db_size output_file" << endl;
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
    output_file.open(argv[5]);

    string sp;

    sp_file >> sp;

    while (!sp_file.eof() ) {
        ++sp_size;
        sp_file >> sp;
    }

    sp_file.close();

    ets_spdb = new ETS_entry[spdb_size];

    for (unsigned i = 0; i < spdb_size; ++i) {
        spdb_file >> ets_spdb[i].ets_seed;

        for (unsigned j = 0; j < sp_size; ++j) {
            spdb_file >> ets_spdb[i].distance[j];
        }
    }

    cout << "the full database is loaded" << endl;

    ETS_entry input_entry;



    while (!file_heap.empty() ) {
        unsigned file_id = file_heap.top();
        file_heap.pop();

        output_file << file_ets_ary[file_id].ets_seed << " ";

        for (unsigned i = 0; i < sp_size; ++i)
            output_file << file_ets_ary[file_id].distance[i] << " ";

        output_file << endl;

        *file_ary[file_id] >> file_ets_ary[file_id].ets_seed;

        for (unsigned i = 0; i < sp_size; ++i) {
            *file_ary[file_id] >> file_ets_ary[file_id].distance[i];
        }

        if (!file_ary[file_id]->eof() )
            file_heap.push(file_id);
    }

    for (unsigned i = 0; i < file_ary.size(); ++i)
        file_ary[i]->close();

    output_file.close();

    return 0;
}

