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
#include <boost/range/iterator_range.hpp>
#include <boost/filesystem/fstream.hpp>
#include <omp.h>
#include "common.h"

#define LEFT_SUB        0b00000001
#define RIGHT_SUB       0b00000010
#define LEFT_SUP        0b00000100
#define RIGHT_SUP       0b00001000
#define BOTH_SUB        0b00010000
#define BOTH_SUP        0b00100000
#define LEFT_OVERLAP    0b01000000
#define RIGHT_OVERLAP   0b10000000
#define ERROR_FLAG      0b00000000

#define __32_64_DEPTH__  15
#define __TASK_SIZE__ 32
#define __BATCH_SIZE__ 8192

using namespace std;
using namespace boost::filesystem;

struct TreeNode {
    char            bp;
    uint64_t        child_offset;
    uint8_t         child_size;
};

struct EDLinkNode {
    uint64_t        F_offset;
    uint64_t        F_size;
};

struct F_content {
    uint64_t    n_id;
    uint8_t     n_distance;
    uint8_t     flag;
};

struct Level_struct {
    uint64_t                            offset;
    vector<F_content>                   F_vec;
    vector<EDLinkNode>                  L_vec;
};

struct Thread_struct {
    vector<uint64_t>                    vid_vec;
    vector<F_content>                   F_vec;
    vector<EDLinkNode>                  L_vec;
};

uint8_t ED_threshold;
uint8_t peak_depth;
unsigned continue_depth;
string result_dir;
boost::filesystem::ofstream L_result_file;
boost::filesystem::ofstream F_result_file;

vector<TreeNode>            *T_vec;
Level_struct                *level_obj;
unsigned                    max_thread_num;
Thread_struct               *thread_obj;
uint64_t                    *L_progress;

// Temp files
std::fstream                L_files[2];
std::fstream                F_files[2];
string                      L_names[2] = {"L_0.ets", "L_1.ets"};
string                      F_names[2] = {"F_0.ets", "F_1.ets"};
uint64_t                    L_count[2];

uint8_t                     in_file_ptr;

deque<pair<uint64_t, uint8_t> > T_queue;

int parseLine(char* line) {
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getValue() { //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

unsigned getQueueIdx(uint64_t node_id) {
    for (unsigned i = 1; i < T_queue.size(); ++i) {
        if (T_queue[i].first > node_id)
            return i - 1;
    }
    
    return T_queue.size() - 1;
}

TreeNode* getTreeNode(uint64_t node_id, int8_t &queue_idx) {
    assert(node_id >= T_queue[0].first);

    if (queue_idx < 0) {
        queue_idx = getQueueIdx(node_id);
    }

    uint8_t vec_idx = T_queue[queue_idx].second;
    uint64_t offset = node_id - T_queue[queue_idx].first;

    //assert(offset <= T_vec[vec_idx].size() );
    
    while (offset >= T_vec[vec_idx].size() ) {
        ++queue_idx;
        vec_idx = T_queue[queue_idx].second;
        offset = node_id - T_queue[queue_idx].first;

        assert(queue_idx < T_queue.size() );
    }

    return T_vec[vec_idx].data() + offset;
}

void initializeEDLinks(Level_struct *cur_level) {
    uint64_t F_size = 0;
   
    for (uint8_t depth = 0; depth < ED_threshold; ++depth) {

        for (uint64_t node_iter = 0; node_iter < T_vec[depth].size(); ++node_iter) {
            cur_level->F_vec.push_back(F_content() );
            cur_level->F_vec.back().n_id = F_size + node_iter;
            cur_level->F_vec.back().n_distance = depth;
            cur_level->F_vec.back().flag = LEFT_SUP | RIGHT_SUP;
            //cout << child_iter << " ";
        }
        // Currect root node flag
        cur_level->F_vec[0].flag = LEFT_SUP | RIGHT_SUP | LEFT_SUB | RIGHT_SUB;

        F_size += T_vec[depth].size();
    }

    cur_level->offset = 0;
    cur_level->L_vec.push_back(EDLinkNode() );
    cur_level->L_vec[0].F_offset = 0;
    cur_level->L_vec[0].F_size = F_size;
}

uint8_t obtainDistance(uint64_t &ptr, uint64_t target_id, uint64_t limit, vector<F_content> &F_vec) {
    while (ptr < limit) {
        if (F_vec[ptr].n_id == target_id)
            return F_vec[ptr].n_distance;
        else if (F_vec[ptr].n_id > target_id)
            return ED_threshold;

        ++ptr;
    }
    return ED_threshold;
}

// Must be used in junction with obtainDistance!
uint8_t obtainFlag(uint64_t ptr, uint64_t target_id, uint64_t limit, vector<F_content> &F_vec) {
    if (ptr >= limit || F_vec[ptr].n_id != target_id)
        return ERROR_FLAG;
    else
        return F_vec[ptr].flag;
}

void processNode(EDLinkNode *v, uint64_t v_id, unsigned depth, Level_struct *cur_level, Thread_struct *thread_data) {
    TreeNode *v_Tnode, *u_Tnode, *vC_Tnode, *uC_Tnode;
    EDLinkNode *vC;
    uint64_t vC_id, u_id, uC_id;
    uint8_t dv_u, dvC_u, dv_uC, min_d;
    uint8_t flag_u_v, flag_u_vC, flag_uC_v;
    uint64_t v_u_ptr, vC_u_ptr, v_uC_ptr;
    uint64_t v_ptr_limit;
    int8_t v_queue_idx, u_queue_idx, vC_queue_idx, uC_queue_idx;

    v_queue_idx = (depth < ED_threshold) ? depth: ED_threshold - 1;
    vC_queue_idx = v_queue_idx + 1;

    //cout << "v_id: " << v_id << endl;
   
  //!!! get Tree node 
    //v_Tnode = T_vec->data() + v_id;
    v_Tnode = getTreeNode(v_id, v_queue_idx);

    for (uint64_t vc_i = 0; vc_i < v_Tnode->child_size; ++vc_i) {

        // Initialize vC pointer
        vC_id = v_Tnode->child_offset + vc_i;
        //vC = L_vec->data() + vC_id;
        thread_data->vid_vec.push_back(vC_id);
        thread_data->L_vec.push_back(EDLinkNode() );
        vC = thread_data->L_vec.data() + (thread_data->L_vec.size() - 1);
        vC->F_offset = thread_data->F_vec.size();
       
      //!!! 
        //vC_Tnode = T_vec->data() + vC_id;
        vC_Tnode = getTreeNode(vC_id, vC_queue_idx);
        
        // Reset the ED pointers
        v_u_ptr = v->F_offset;
        v_uC_ptr = v->F_offset;
        vC_u_ptr = vC->F_offset;

        // Reset u_queue_idx
        u_queue_idx = -1;

        v_ptr_limit = v->F_offset + v->F_size;

        // Add root node, if necessary
        if (depth + 1 < ED_threshold) {
            //cout << "pushing 0 with depth: " << (unsigned) T_depth << endl;
            thread_data->F_vec.push_back(F_content() );
            thread_data->F_vec.back().n_id = 0;
            thread_data->F_vec.back().n_distance = depth + 1;
            thread_data->F_vec.back().flag = LEFT_SUB | RIGHT_SUB;

            ++(vC->F_size);
        }

        //cout << "vC_id: " << vC_id << endl;

        for (uint64_t u_i = 0; u_i < v->F_size; ++u_i) {
            // Get dvCu, dvu and move the pointers
            u_id = cur_level->F_vec[v->F_offset + u_i].n_id;
            dv_u  = obtainDistance(v_u_ptr, u_id, v_ptr_limit, cur_level->F_vec);
            flag_u_v = obtainFlag(v_u_ptr, u_id, v_ptr_limit, cur_level->F_vec);

            //cout << u_id << "-"; 
    
            assert(dv_u < ED_threshold);

            dvC_u = obtainDistance(vC_u_ptr, u_id, thread_data->F_vec.size(),  thread_data->F_vec);
            flag_u_vC = obtainFlag(vC_u_ptr, u_id, thread_data->F_vec.size(),  thread_data->F_vec);

          //!!!
            //u_Tnode = T_vec->data() + u_id;
            u_Tnode = getTreeNode(u_id, u_queue_idx);
            
            //cout << "child_size: " << (unsigned) u_Tnode->child_size << "!";

            for (uint64_t uC_i = 0; uC_i < u_Tnode->child_size; ++uC_i) {

                uC_id = u_Tnode->child_offset + uC_i;
                dv_uC = obtainDistance(v_uC_ptr, uC_id, v_ptr_limit, cur_level->F_vec);
                flag_uC_v = obtainFlag(v_uC_ptr, uC_id, v_ptr_limit, cur_level->F_vec);

              //!!!
                //uC_Tnode = T_vec->data() + uC_id;
                uC_queue_idx = u_queue_idx + 1;
                uC_Tnode = getTreeNode(uC_id, uC_queue_idx);

                //cout << (unsigned) uC_id << "-";

                min_d = dv_u;

                if (vC_Tnode->bp != uC_Tnode->bp)
                    ++min_d;

                if (dvC_u + 1 < min_d)
                    min_d = dvC_u + 1;
                
                if (dv_uC + 1 < min_d)
                    min_d = dv_uC + 1;

                //cout << (unsigned) min_d << " ";

                // Check and calculate the flag
                //cout << "v_id " << (unsigned) v_id << " u_id " << (unsigned) u_id << " flag_u_v " << (unsigned) flag_u_v << endl;
                uint8_t flag = 0;
                if ( (flag_u_v & (LEFT_SUB | LEFT_SUP) ) && vC_Tnode->bp == uC_Tnode->bp) {
                    flag = flag | (flag_u_v & (LEFT_SUB | LEFT_SUP) );
                }
                if ( (flag_uC_v & RIGHT_SUB && min_d < ED_threshold) || uC_id == vC_id) {
                    flag = flag | RIGHT_SUB;
                }
                if (flag_u_vC & RIGHT_SUP || uC_id == vC_id) {
                    flag = flag | RIGHT_SUP;
                }
                if (v_id != 0 && uC_id != v_id && flag_uC_v & (BOTH_SUB | LEFT_SUB)) {
                    flag = flag | BOTH_SUB;
                }
                if (u_id != 0 && u_id != vC_id && flag_u_vC & (BOTH_SUP | LEFT_SUP)) {
                    flag = flag | BOTH_SUP;
                }
                if (v_id != 0 && uC_id != v_id && (flag_uC_v & (LEFT_OVERLAP | LEFT_SUP) ) ) {
                    flag = flag | LEFT_OVERLAP;
                }
                if (u_id != 0 && u_id != vC_id && (flag_u_vC & (LEFT_SUB | RIGHT_OVERLAP) ) )  {
                    flag = flag | RIGHT_OVERLAP;
                }

                if (min_d < ED_threshold) {
                    thread_data->F_vec.push_back(F_content() );
                    thread_data->F_vec.back().n_id = uC_id;
                    thread_data->F_vec.back().n_distance = min_d;
                    thread_data->F_vec.back().flag = flag;

                    ++(vC->F_size);
                }
            }
            //cout << endl;
        }

    }
}

void printNode(EDLinkNode *v, uint64_t v_id, Level_struct *level, uint8_t depth) {
    //result_file << T_offset << ": ";
    L_result_file.write( (char*) &v_id, sizeof(uint64_t) );
    L_result_file.write( (char*) &v->F_size, sizeof(uint64_t) );
    
    for (uint64_t F_iter = 0; F_iter < v->F_size; ++F_iter) {
        //result_file << (unsigned long long) level->F_vec[v->F_offset + F_iter].n_id << " " << (unsigned) level->F_vec[v->F_offset + F_iter].n_distance << " " << (unsigned) level->F_vec[v->F_offset + F_iter].flag << " ";
        /*
        */
        uint64_t neighbor_id = level->F_vec[v->F_offset + F_iter].n_id;
        uint32_t neighbor_id_short = neighbor_id;
        uint8_t neighbor_distance = level->F_vec[v->F_offset + F_iter].n_distance;
        uint8_t neighbor_flag = level->F_vec[v->F_offset + F_iter].flag;

        if (depth >= __32_64_DEPTH__)
            F_result_file.write( (char*) &neighbor_id, sizeof(uint64_t) );
        else
            F_result_file.write( (char*) &neighbor_id_short, sizeof(uint32_t) );

        F_result_file.write( (char*) &neighbor_distance, sizeof(uint8_t) );
        F_result_file.write( (char*) &neighbor_flag, sizeof(uint8_t) );
    }
    //result_file << endl;
}

void printNodeToTemp(EDLinkNode *v,  vector<F_content> &F_vec, uint8_t file_ptr) {
    uint64_t F_size = v->F_size;
    uint64_t F_offset = v->F_offset;

    //L_files[file_ptr] << F_size << " ";
    L_files[file_ptr].write( (char*) &F_size, sizeof(uint64_t) );

    for (uint64_t F_idx = 0; F_idx < F_size; ++F_idx) {
        uint64_t node_id = F_vec[F_offset + F_idx].n_id;
        uint8_t node_distance = F_vec[F_offset + F_idx].n_distance;
        uint8_t node_flag = F_vec[F_offset + F_idx].flag;

        F_files[file_ptr].write( (char*) &node_id, sizeof(uint64_t) );
        F_files[file_ptr].write( (char*) &node_distance, sizeof(uint8_t) );
        F_files[file_ptr].write( (char*) &node_flag, sizeof(uint8_t) );
        
        //F_files[file_ptr] << (unsigned long long) F_vec[F_offset + F_idx].n_id << " " << (unsigned) F_vec[F_offset + F_idx].n_distance << " " <<  (unsigned) F_vec[F_offset + F_idx].flag << " ";
    }
}


void populateTreeVec(uint8_t vec_idx, boost::filesystem::ifstream &in_file, uint64_t count) {

    T_vec[vec_idx].clear();

    //cout << in_file.fail() << endl;

    for (uint64_t node_iter = 0; node_iter < count; ++node_iter) {
        T_vec[vec_idx].push_back(TreeNode() );

        in_file.read(&T_vec[vec_idx][node_iter].bp, sizeof(char) );
        in_file.read( (char*) &T_vec[vec_idx][node_iter].child_offset, sizeof(uint64_t) );
        in_file.read( (char*) &T_vec[vec_idx][node_iter].child_size, sizeof(uint8_t) );
        /*
        long long unsigned temp;

        in_file >> T_vec[vec_idx][node_iter].bp;
        in_file >> temp;
        T_vec[vec_idx][node_iter].child_offset = temp;
        in_file >> temp;
        T_vec[vec_idx][node_iter].child_size = temp;
        */
    }
}

void processLevel(unsigned depth) {
    uint8_t out_file_ptr = 1 - in_file_ptr;
    uint64_t batch_count = (L_count[in_file_ptr] + __BATCH_SIZE__ - 1) / __BATCH_SIZE__;
    L_count[out_file_ptr] = 0;
    uint64_t child_offset = level_obj->offset + L_count[in_file_ptr];

    // Open all files
    L_files[in_file_ptr].open(L_names[in_file_ptr], std::fstream::in | std::fstream::binary);
    L_files[out_file_ptr].open(L_names[out_file_ptr], std::fstream::out | std::fstream::binary);
    F_files[in_file_ptr].open(F_names[in_file_ptr], std::fstream::in | std::fstream::binary);
    F_files[out_file_ptr].open(F_names[out_file_ptr], std::fstream::out | std::fstream::binary);
    L_result_file.open(result_dir / (to_string(depth + 1) + ".L") );//, boost::filesystem::ofstream::binary);
    F_result_file.open(result_dir / (to_string(depth + 1) + ".F") );//, boost::filesystem::ofstream::binary);

    // Iterate through all batches
    for (uint64_t batch_idx = 0; batch_idx < batch_count; ++batch_idx) {

        // Clear per thread storage first
        for (unsigned t_id = 0; t_id < max_thread_num; ++t_id) {
            thread_obj[t_id].vid_vec.clear();
            thread_obj[t_id].L_vec.clear();
            thread_obj[t_id].F_vec.clear();

            if (depth >= peak_depth && depth < peak_depth + ED_threshold) {
                vector<uint64_t>().swap(thread_obj[t_id].vid_vec);
                vector<F_content>().swap(thread_obj[t_id].F_vec);
                vector<EDLinkNode>().swap(thread_obj[t_id].L_vec);
            }

            L_progress[t_id] = 0;
        }

        // Feed L_vec and F_vec
        uint64_t batch_size = (batch_idx == batch_count - 1) ? L_count[in_file_ptr] - __BATCH_SIZE__ * batch_idx : __BATCH_SIZE__;
        uint64_t F_offset = 0;
        level_obj->L_vec.clear();
        level_obj->F_vec.clear();

        for (uint64_t L_idx = 0; L_idx < batch_size; ++L_idx) {
            level_obj->L_vec.push_back(EDLinkNode() );
            level_obj->L_vec.back().F_offset = F_offset;
            uint64_t F_size;
            //L_files[in_file_ptr] >> F_size;
            L_files[in_file_ptr].read( (char*) &F_size, sizeof(uint64_t) );
            level_obj->L_vec.back().F_size = F_size;
            F_offset += F_size;

            for (uint64_t F_idx = 0; F_idx < F_size; ++F_idx) {
                level_obj->F_vec.push_back(F_content() );
                /*
                F_files[in_file_ptr] >> level_obj->F_vec.back().n_id;
                unsigned temp;
                F_files[in_file_ptr] >> temp;
                level_obj->F_vec.back().n_distance = temp;
                F_files[in_file_ptr] >> temp;
                level_obj->F_vec.back().flag = temp;
                */
                F_files[in_file_ptr].read( (char*) &level_obj->F_vec.back().n_id, sizeof(uint64_t) );
                F_files[in_file_ptr].read( (char*) &level_obj->F_vec.back().n_distance, sizeof(uint8_t) );
                F_files[in_file_ptr].read( (char*) &level_obj->F_vec.back().flag, sizeof(uint8_t) );
            }
        }

        // Process all nodes of the level
        #pragma omp parallel for schedule(static, __TASK_SIZE__)
        for (uint64_t node_iter = 0; node_iter < level_obj->L_vec.size(); ++node_iter) {
            unsigned t_id = omp_get_thread_num();

            //cout << "t_id: " << t_id << endl;

            processNode(level_obj->L_vec.data() + node_iter, level_obj->offset + node_iter, depth, level_obj, thread_obj + t_id);
        }
        level_obj->offset += level_obj->L_vec.size();
        
        //cout << "done with parallel" << endl;

        /* Merge results from a all threads */
        // Clear level_obj and intialize priority queue
        level_obj->L_vec.clear();
        level_obj->F_vec.clear();

        auto cmp = [](unsigned left, unsigned right) { return (thread_obj[left].vid_vec[L_progress[left] ] > thread_obj[right].vid_vec[L_progress[right] ]);};
        priority_queue<unsigned, vector<unsigned>, decltype(cmp)> tid_heap (cmp);

        // First fill the heap
        for (unsigned t_id = 0; t_id < max_thread_num; ++t_id) {
            if (thread_obj[t_id].vid_vec.size() > 0)
                tid_heap.push(t_id);
        }

        // Pop heap while fill the level_obj vector
        uint64_t L_count_local = 0;
        F_offset = 0;

        while (!tid_heap.empty() ) {
            unsigned t_id = tid_heap.top();
            tid_heap.pop();

            // Push to the level vector
            while (L_progress[t_id] < thread_obj[t_id].vid_vec.size() && L_count_local + child_offset == thread_obj[t_id].vid_vec[L_progress[t_id] ] ) {
                //cout << "L_cout: "  << L_count_local << endl;

                uint64_t t_F_offset = thread_obj[t_id].L_vec[L_progress[t_id]].F_offset;

                level_obj->L_vec.push_back(thread_obj[t_id].L_vec[L_progress[t_id]] );
                level_obj->L_vec[L_count_local].F_offset = F_offset;

                for (uint64_t F_iter = 0; F_iter < level_obj->L_vec[L_count_local].F_size; ++F_iter)
                    level_obj->F_vec.push_back(thread_obj[t_id].F_vec[t_F_offset + F_iter]);
                
                F_offset += level_obj->L_vec[L_count_local].F_size;
                ++L_progress[t_id];
                ++L_count_local;
            }

            //cout << "t_id: " << t_id << endl;
            //cout << "progress: " <<  L_progress[t_id] << endl;
            //cout << "v_id: " << thread_obj[t_id].vid_vec[L_progress[t_id] ] << endl;

            if (L_progress[t_id] < thread_obj[t_id].vid_vec.size() )
                tid_heap.push(t_id);
        }

        //cout << "done with refilling" << endl;

        // Print out children level
        for (uint64_t node_iter = 0; node_iter < level_obj->L_vec.size(); ++node_iter) {
            printNode(level_obj->L_vec.data() + node_iter, child_offset + node_iter, level_obj, depth);
            printNodeToTemp(level_obj->L_vec.data() + node_iter, level_obj->F_vec, out_file_ptr);
        }
        L_count[out_file_ptr] += level_obj->L_vec.size();
        child_offset += level_obj->L_vec.size();
    }

    // Close all files 
    L_files[in_file_ptr].close();
    L_files[out_file_ptr].close();
    F_files[in_file_ptr].close();
    F_files[out_file_ptr].close();
    L_result_file.close();
    F_result_file.close();
    in_file_ptr = 1 - in_file_ptr;
}

void initialize(boost::filesystem::ifstream &info_file, string path) {
    T_vec = new vector<TreeNode> [2 * ED_threshold];
    level_obj = new Level_struct;
    uint64_t total_size = 0;
    max_thread_num = omp_get_max_threads();

    L_progress = new uint64_t [max_thread_num];
    thread_obj = new Thread_struct [max_thread_num];

    level_obj->offset = 0;

    for (unsigned level_iter = 0; level_iter < 2 * ED_threshold; ++level_iter) {
        string level_name;
        unsigned long long level_size;
        info_file >> level_name;
        info_file >> level_size;

        if (stoi(level_name) < (int) continue_depth)
            level_obj->offset += level_size;

        boost::filesystem::ifstream in_file(path / level_name, std::fstream::binary);

        if (!info_file.eof() ) {
            assert(level_iter == stoi(level_name) );

            populateTreeVec(level_iter, in_file, level_size);

            in_file.close();

            T_queue.push_back(make_pair(total_size, level_iter) );
            total_size += level_size;
        }
        // This is the last level! Revise children count
        else {
            uint8_t last_vec_idx = level_iter - 1;
            
            for (uint64_t node_iter = 0; node_iter < T_vec[last_vec_idx].size(); ++node_iter)
                T_vec[last_vec_idx][node_iter].child_size = 0;
                
            in_file.close();

            break;
        }
    }

    if (continue_depth == 0) {
        initializeEDLinks(level_obj);

        L_result_file.open(result_dir / "0.L");//, boost::filesystem::ofstream::binary);
        F_result_file.open(result_dir / "0.F");//, boost::filesystem::ofstream::binary);
        printNode(level_obj->L_vec.data(), 0, level_obj, 0);
        L_result_file.close();
        F_result_file.close();

        in_file_ptr = 0;
        L_files[in_file_ptr].open(L_names[in_file_ptr], std::fstream::out);
        F_files[in_file_ptr].open(F_names[in_file_ptr], std::fstream::out);

        printNodeToTemp(level_obj->L_vec.data(), level_obj->F_vec, in_file_ptr);

        L_files[in_file_ptr].close();
        F_files[in_file_ptr].close();

        L_count[in_file_ptr] = 1;
    }
}

void populateStream(unsigned depth, boost::filesystem::ifstream &L_file, boost::filesystem::ifstream &F_file, uint64_t L_size) {
    L_files[in_file_ptr].open(L_names[in_file_ptr], std::fstream::out | std::fstream::binary);
    F_files[in_file_ptr].open(F_names[in_file_ptr], std::fstream::out | std::fstream::binary);
    L_count[in_file_ptr] = L_size;

    for (uint64_t L_idx = 0; L_idx < L_size; ++L_idx) {
        uint64_t node_id;
        uint64_t neighbor_size;
        uint64_t neighbor_id;
        uint32_t neighbor_id_short;
        uint8_t  neighbor_distance;
        uint8_t  neighbor_flag;

        L_file.read( (char*) &node_id, sizeof(uint64_t) );
        L_file.read( (char*) &neighbor_size, sizeof(uint64_t) );
        L_files[in_file_ptr].write( (char*) &neighbor_size, sizeof(uint64_t) );

        for (uint64_t neighbor_idx = 0; neighbor_idx < neighbor_size; ++neighbor_idx) {
            if (depth >= __32_64_DEPTH__)
                F_file.read( (char*) &neighbor_id, sizeof(uint64_t) );
            else {
                F_file.read( (char*) &neighbor_id_short, sizeof(uint32_t) );
                neighbor_id = neighbor_id_short;
            }

            F_file.read( (char*) &neighbor_distance, sizeof(uint8_t) );
            F_file.read( (char*) &neighbor_flag, sizeof(uint8_t) );

            F_files[in_file_ptr].write( (char*) &neighbor_id, sizeof(uint64_t) );
            F_files[in_file_ptr].write( (char*) &neighbor_distance, sizeof(uint8_t) );
            F_files[in_file_ptr].write( (char*) &neighbor_flag, sizeof(uint8_t) );
        }
    }

    L_files[in_file_ptr].close();
    F_files[in_file_ptr].close();
}

int main (int argc, char *argv[]) {
    if (argc != 6 && argc != 7) {
        cout << "Usage: ./bin input_folder ED_threshold max_depth result_file peak_depth" << endl;
        cout << "Usage: ./bin input_folder ED_threshold max_depth result_file peak_depth continue_depth" << endl;
        return 0;
    }

    cout << "Processing the tree in " << argv[1] << " with " << omp_get_max_threads() << " threads." << endl;

    string path = argv[1];
    ED_threshold = atoi(argv[2]);
    unsigned max_depth = atoi(argv[3]);
    result_dir = argv[4];
    peak_depth = atoi(argv[5]);
    continue_depth = 0;
   
    if (argc == 7) {
        continue_depth = atoi(argv[5]);
        if (continue_depth < ED_threshold) {
            cerr << "Error: continue_depth < ED_threshold!" << endl;
            exit(0);
        }
    }

    if (ED_threshold == 0) {
        cout << "Error: ED_threshold must be greater than 0!" << endl;
        return 0;
    }

    if (exists(result_dir) ) {
        if (continue_depth == 0) {
            for(auto& entry : boost::make_iterator_range(directory_iterator(result_dir), {})) {
                cout << "Removing" << entry << "\n";
                remove(entry);
            }
        }
    }
    else
        create_directory(result_dir);

    boost::filesystem::ifstream info_file(path / "info");

    initialize(info_file, argv[1]);

    bool cleanup = false;

    unsigned continue_L_size;

    for (unsigned depth = 0; depth < max_depth; ++depth) {
        cout << "depth: " << depth << endl;
        cout << "Currently used memory: " << getValue() << "Kb" << endl;

        // Eject from tree queue and load in new tree (optional) for children level 
        if (depth >= ED_threshold) {
            // Eject
            pair<uint64_t, uint8_t> free_node = T_queue.front();
            T_queue.pop_front();

            // Load (optional)
            boost::filesystem::ifstream in_file;
            unsigned long long level_size = 0;
            string level_name;

            info_file >> level_name;
            info_file >> level_size;
            //cout << " level_name: " << level_name << endl;

            if (stoi(level_name) == (int) continue_depth)
                continue_L_size = level_size;

            // Proceed to load tree
            if (!info_file.eof() ) {
                in_file.open(path / level_name, std::fstream::binary);

                free_node.first = T_queue.back().first + T_vec[T_queue.back().second].size();
                T_queue.push_back(free_node);

                uint8_t vec_idx = free_node.second;

                // Populate new tree vector
                populateTreeVec(vec_idx, in_file, level_size);

                in_file.close();
            }
            // The first time we run out of levels! Clear children size!
            else if (T_queue.size() == 2 * ED_threshold - 1) {
                cout << "Tree exhausted at depth: " << depth /*<< " T_queue size:" << T_queue.size()*/ << endl;

                uint8_t last_vec_idx = T_queue.back().second;

                //cout << "last_vec_idx: " << (unsigned) last_vec_idx << endl;
                
                for (uint64_t node_iter = 0; node_iter < T_vec[last_vec_idx].size(); ++node_iter)
                    T_vec[last_vec_idx][node_iter].child_size = 0;

                cleanup = true;
            }

            if (stoi(level_name) < (int) continue_depth)
                level_obj->offset += level_size;
        }

        if (continue_depth != 0 && depth == continue_depth) {
            boost::filesystem::ifstream L_file(result_dir / (to_string(depth) + ".L") );
            boost::filesystem::ifstream F_file(result_dir / (to_string(depth) + ".F") );
            populateStream(depth, L_file, F_file, continue_L_size);
            L_file.close();
            F_file.close();
        }

        if (depth >= continue_depth)
            processLevel(depth);

        // Quit when no more level is available
        if (T_queue.size() == 2)
            break;
    }
  
    info_file.close();
    //result_file.close();

    delete [] T_vec;
    delete level_obj;
    delete [] L_progress;
    delete [] thread_obj;


    return 0;
}

