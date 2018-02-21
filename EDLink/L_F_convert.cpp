#include<utility>
#include<fstream>
#include<iostream>
#include<cassert>
#include<cstdint>
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include<vector>
#include <queue>
#include<deque>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <omp.h>
#include "common.h"

using namespace std;
using namespace boost::filesystem;

boost::filesystem::ifstream L_result_file;
boost::filesystem::ifstream F_result_file;

void convertNode(uint8_t depth) {
    //result_file << T_offset << ": ";
    uint64_t node_id;
    uint64_t F_size;

    L_result_file.read( (char*) &node_id, sizeof(uint64_t) );
    L_result_file.read( (char*) &F_size, sizeof(uint64_t) );

    cout << node_id << " " << F_size << ": ";
    
    for (uint64_t F_iter = 0; F_iter < F_size; ++F_iter) {
        /*
        */
        uint64_t neighbor_id;
        uint32_t neighbor_id_short;
        uint8_t neighbor_distance;
        uint8_t neighbor_flag;

        if (depth >= __32_64_DEPTH__)
            F_result_file.read( (char*) &neighbor_id, sizeof(uint64_t) );
        else {
            F_result_file.read( (char*) &neighbor_id_short, sizeof(uint32_t) );
            neighbor_id = neighbor_id_short;
        }

        F_result_file.read( (char*) &neighbor_distance, sizeof(uint8_t) );
        F_result_file.read( (char*) &neighbor_flag, sizeof(uint8_t) );
        
        cout << (unsigned long long) neighbor_id << "-" << (unsigned) neighbor_distance << "|" << (unsigned) neighbor_flag << " ";
    }
    cout << endl;
}

void processLevel(unsigned depth, uint64_t L_size) {
    // Iterate through all batches
    for (uint64_t L_iter = 0; L_iter < L_size; ++L_iter) {
        convertNode(depth);
    }
}

int main (int argc, char *argv[]) {
    if (argc != 4) {
        cout << "Usage: ./bin tree_file result_file max_depth" << endl;
        return 0;
    }

    string tree_dir = argv[1];
    string result_dir = argv[2];
    unsigned max_depth = atoi(argv[3]);

    boost::filesystem::ifstream info_file(tree_dir / "info");

    string level_name;
    unsigned long long level_size;

    for (unsigned depth = 0; depth <= max_depth; ++depth) {
        cout << "depth: " << depth << endl;

        info_file >> level_name;
        info_file >> level_size;

        // Proceed to load tree
        if (!info_file.eof() ) {
            L_result_file.open(result_dir / (level_name + ".L") );
            F_result_file.open(result_dir / (level_name + ".F") );

            processLevel(depth, level_size);

            L_result_file.close();
            F_result_file.close();
        }
    }
  
    info_file.close();

    return 0;
}

