#include <string>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/filesystem/fstream.hpp>
#include <sdsl/suffix_trees.hpp>
#include "common.h"

using namespace std;
using namespace boost::filesystem;
using namespace sdsl;

cst_sada<> suffix_tree;
#define TYPE   decltype(suffix_tree.root())

template<typename T>
vector<T> *iter_vec;

int main (int argc, char *argv[]) {
    if (argc != 7) {
        cout << "Usage: ./bin tree_path neighbor_dir max_depth neighbor_distance profile_config result_dir" << endl;
        return 0;
    }

    string tree_path = argv[1];
    string neighbor_path = argv[2];
    unsigned max_depth = atoi(argv[3]);
    unsigned max_distance = atoi(argv[4]);
    string profile_config = argv[5];
    string result_path = argv[6];

    if (exists(result_path) ) {
        for(auto& entry : boost::make_iterator_range(directory_iterator(result_path), {})) {
            cout << "Removing" << entry << "\n";
            remove(entry);
        }
    }
    else
        create_directory(result_path);

    if (!load_from_file(suffix_tree, tree_path) )
    {
        cerr << "ERROR: Could not load the suffix_tree from file!\n";
        return -1;
    }

    //Initialize
    iter_vec<TYPE> = new vector<TYPE> [2];
    unsigned vec_idx = 0;
    unsigned child_vec_idx = 1;
    iter_vec<TYPE>[vec_idx].push_back(suffix_tree.root() );
    unsigned depth_idx = 0;
    vector<uint8_t> *ets_vec = new vector<uint8_t>;
    uint64_t node_count = 0;

    vector<unsigned> depth_vec;
    std::ifstream profile_config_file(profile_config);
    unsigned profile_depth;
    profile_config_file >> profile_depth;
    while (!profile_config_file.eof() ) {
        depth_vec.push_back(profile_depth);
        profile_config_file >> profile_depth;
    }

    ets_vec->resize(suffix_tree.size(suffix_tree.root() ) );
    cout << "ets_vec size " << ets_vec->size() << endl;

    boost::filesystem::ifstream L_file;
    boost::filesystem::ifstream F_file;
    boost::filesystem::ofstream ETS_file;

    for (unsigned depth = 0; depth <= max_depth; ++depth) {
        if (depth == depth_vec[depth_idx]) {
            for (uint64_t vec_idx = 0; vec_idx < ets_vec->size(); ++vec_idx)
                ets_vec->at(vec_idx) = 0;

            L_file.open(neighbor_path / (to_string(depth) + ".L"), std::fstream::binary);
            F_file.open(neighbor_path / (to_string(depth) + ".F"), std::fstream::binary);
            ETS_file.open(result_path / (to_string(depth) + ".ets") );
        }

        cout << "Processing depth " << depth << endl;

        iter_vec<TYPE>[child_vec_idx].clear();

        for (uint64_t iter_idx = 0; iter_idx < iter_vec<TYPE>[vec_idx].size(); ++iter_idx) {

            TYPE cur_iter = iter_vec<TYPE>[vec_idx][iter_idx];

            uint64_t node_id;
            uint64_t neighbor_size;
            uint64_t neighbor_id;
            uint32_t neighbor_id_short;
            uint8_t  neighbor_distance;
            uint8_t  neighbor_flag;

            uint8_t  min_distance = max_distance;

            if (depth == depth_vec[depth_idx]) {
                L_file.read( (char*) &node_id, sizeof(uint64_t) );
                L_file.read( (char*) &neighbor_size, sizeof(uint64_t) );

                assert(node_id == node_count);
           
                for (uint64_t neighbor_idx = 0; neighbor_idx < neighbor_size; ++neighbor_idx) {
                    if (depth >= __32_64_DEPTH__)
                        F_file.read( (char*) &neighbor_id, sizeof(uint64_t) );
                    else
                        F_file.read( (char*) &neighbor_id_short, sizeof(uint32_t) );

                    F_file.read( (char*) &neighbor_distance, sizeof(uint8_t) );
                    F_file.read( (char*) &neighbor_flag, sizeof(uint8_t) );

                    if (neighbor_flag == 0 && neighbor_distance < min_distance)
                        min_distance = neighbor_distance;
                }

                for (uint64_t leaf_idx = suffix_tree.lb(cur_iter); leaf_idx <= suffix_tree.rb(cur_iter); ++leaf_idx) {
                    auto leaf_node = suffix_tree.select_leaf(leaf_idx + 1);
                    uint64_t loc = suffix_tree.sn(leaf_node);
                    ets_vec->at(loc) = min_distance;
                }
            }

            // Have multiple children
            if (suffix_tree.depth(cur_iter) == depth && !suffix_tree.is_leaf(cur_iter) ) {
                
                // Add each child as long as first letter isn't '|'
                for (unsigned child_idx = 1; child_idx <= suffix_tree.degree(cur_iter); ++child_idx) {
                    TYPE child_iter = suffix_tree.select_child(cur_iter, child_idx);

                    if (suffix_tree.edge(child_iter, depth + 1) != '|' && suffix_tree.edge(child_iter, depth + 1) != '\0') {
                        iter_vec<TYPE>[child_vec_idx].push_back(child_iter);
                    }
                }

                //vertex_file << (unsigned long long) child_offset << " "<< child_size << endl;
            }
            else if (suffix_tree.depth(cur_iter) != depth && suffix_tree.edge(cur_iter, depth + 1) != '|') {
                iter_vec<TYPE>[child_vec_idx].push_back(cur_iter);
                //vertex_file << (unsigned long long) child_offset << " 1" << endl;
            }
            else {
                // Simply terminate all further extensions upon eaching the end of a leaf or '|'
                //vertex_file << "0 0" << endl;
            }

            ++node_count;
        }

        vec_idx = 1 - vec_idx;
        child_vec_idx = 1 - vec_idx;

        if (depth == depth_vec[depth_idx]) {
            for (uint64_t loc_idx = 0; loc_idx < ets_vec->size(); ++loc_idx)
                ETS_file << (unsigned) ets_vec->at(loc_idx) << " ";
        }

        if (depth == depth_vec[depth_idx]) {
            L_file.close();
            F_file.close();
            ETS_file.close();

            ++depth_idx;

            if (depth_idx == depth_vec.size() )
                break;
        }
    }
}

