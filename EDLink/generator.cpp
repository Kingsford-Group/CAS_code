#include <iostream>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <vector>
#include <fstream>
#include <string>
#include <sdsl/suffix_trees.hpp>
#include "RefDB.h"

using namespace seqan;
using namespace std;
using namespace boost::filesystem;
using namespace sdsl;

struct CmdOptions {
    bool build;
    bool load;
    bool process_ref;
    unsigned start_depth;
    unsigned max_depth;
    seqan::CharString fastaFileName;
    string tree_dir;

    CmdOptions() :
        build(false), load(false), process_ref(false) {}

};

cst_sada<> suffix_tree;
#define TYPE   decltype(suffix_tree.root())

template<typename T>
vector<T> *iter_vec;

CmdOptions options;

void profileReference(unsigned max_depth) {
    if (options.start_depth == 0) {
        remove_all(options.tree_dir);
        create_directory(options.tree_dir);
    }

    uint8_t vec_idx = 0;
    uint8_t child_vec_idx = 1;
    uint64_t child_offset = 1;
    boost::filesystem::ofstream info_file(options.tree_dir / "info" );
    boost::filesystem::ofstream vertex_file;
    boost::filesystem::ofstream vertex_freq_file;


    //Initialize
    iter_vec<TYPE>[vec_idx].push_back(suffix_tree.root() );

    for (unsigned depth = 0; depth <= max_depth; ++depth) {
        cout << "Processing depth " << depth << endl;

        iter_vec<TYPE>[child_vec_idx].clear();

        if (depth >= options.start_depth) {
            path p_vertex(options.tree_dir / to_string(depth) );
            path p_freq(options.tree_dir / (to_string(depth) + ".freq") );
            vertex_file.open(p_vertex, std::fstream::binary);
            vertex_freq_file.open(p_freq, std::fstream::binary);
            //vertex_file.open(p);
        }

        for (uint64_t iter_idx = 0; iter_idx < iter_vec<TYPE>[vec_idx].size(); ++iter_idx) {

            TYPE cur_iter = iter_vec<TYPE>[vec_idx][iter_idx];

            if (depth >= options.start_depth) {
                if (cur_iter == suffix_tree.root() )
                    vertex_file.write("R", strlen("R") );
                    //vertex_file << "R" << endl;
                else {
                    char bp = suffix_tree.edge(cur_iter, depth);
                    vertex_file.write(&bp, sizeof(char) );
                    //vertex_file << suffix_tree.edge(cur_iter, depth) << endl;
                }

                uint64_t frequency = suffix_tree.size(cur_iter);
                vertex_freq_file.write( (char*) &frequency, sizeof(uint64_t) );
            }

            uint8_t child_size = 0;

            // Have multiple children
            if (suffix_tree.depth(cur_iter) == depth && !suffix_tree.is_leaf(cur_iter) ) {
                
                // Add each child as long as first letter isn't '|'
                for (unsigned child_idx = 1; child_idx <= suffix_tree.degree(cur_iter); ++child_idx) {
                    TYPE child_iter = suffix_tree.select_child(cur_iter, child_idx);

                    if (suffix_tree.edge(child_iter, depth + 1) != '|' && suffix_tree.edge(child_iter, depth + 1) != '\0') {
                        iter_vec<TYPE>[child_vec_idx].push_back(child_iter);
                        ++child_size;
                    }
                }

                //vertex_file << (unsigned long long) child_offset << " "<< child_size << endl;
            }
            else if (suffix_tree.depth(cur_iter) != depth && suffix_tree.edge(cur_iter, depth + 1) != '|') {
                iter_vec<TYPE>[child_vec_idx].push_back(cur_iter);
                child_size = 1;
                //vertex_file << (unsigned long long) child_offset << " 1" << endl;
            }
            else {
                // Simply terminate all further extensions upon eaching the end of a leaf or '|'
                //vertex_file << "0 0" << endl;
            }

            vertex_file.write( (char*) &child_offset, sizeof(uint64_t) );
            vertex_file.write( (char*) &child_size, sizeof(uint8_t) );
            child_offset += (uint64_t) child_size;
        }

        info_file << depth << " " << iter_vec<TYPE>[vec_idx].size() << endl;
        cout << "Capacity of the child vector: " << iter_vec<TYPE>[child_vec_idx].capacity() << endl;

        if (depth >= options.start_depth) {
            vertex_file.close();
            vertex_freq_file.close();
        }

        vec_idx = 1 - vec_idx;
        child_vec_idx = 1 - vec_idx;
    }
    
    info_file.close();
}


seqan::ArgumentParser::ParseResult
parseCommandLine(CmdOptions & options, int argc, char const ** argv) {
    seqan::ArgumentParser parser("modify_string");

    addOption(parser, seqan::ArgParseOption("b", "fasta-file", "Build index from the fasta file.", seqan::ArgParseArgument::INPUT_FILE) );
    setValidValues(parser, "fasta-file", "fa fasta");

    addOption(parser, seqan::ArgParseOption("l", "index-file", "Load index files.", seqan::ArgParseArgument::INPUT_FILE) );
    setValidValues(parser, "index-file", "");

    addOption(parser, seqan::ArgParseOption("p", "process-ref", "Process the reference genome. Set depth of the suffix tri.", seqan::ArgParseArgument::INTEGER, "INT") );
    setDefaultValue(parser, "process-ref", "2");

    addOption(parser, seqan::ArgParseOption("d", "tree-dir", "Tree directory.", seqan::ArgParseArgument::STRING) );
    setDefaultValue(parser, "tree-dir", "tree_dir");

    addOption(parser, seqan::ArgParseOption("s", "start-depth", "Start depth to process the reference genome.", seqan::ArgParseArgument::INTEGER, "INT") );
    setDefaultValue(parser, "start-depth", "0");

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
    return res;

    options.build = isSet(parser, "fasta-file");

    if (options.build)
        getOptionValue(options.fastaFileName, parser, "fasta-file");

    options.load = isSet(parser, "index-file");

    if (options.load)
        getOptionValue(options.fastaFileName, parser, "index-file");

    options.process_ref = isSet(parser, "process-ref");

    if (options.process_ref)
        getOptionValue(options.max_depth, parser, "process-ref");

    getOptionValue(options.start_depth, parser, "start-depth");
    
    getOptionValue(options.tree_dir, parser, "tree-dir");

    return seqan::ArgumentParser::PARSE_OK;
}

int buildIndices(CharString fastaFileName) {

    CharString baseFileName = fastaFileName;
    while (back(baseFileName) != '.')
        eraseBack(baseFileName);
    eraseBack(baseFileName);
    string treeFileName = string(toCString(baseFileName)) + ".stree";

    RefDB* refdb= new RefDB;

    refdb->loadRefFile(string(toCString(fastaFileName)));

    string* global_ref = refdb->globalReference();
    cout << "global_ref length: " << global_ref->length() << endl;
    cout << "last character: " << (*global_ref)[global_ref->length() - 1] << endl;

    construct_im(suffix_tree, *global_ref, 1);

    if (!store_to_file(suffix_tree, treeFileName) )
    {
        cerr << "ERROR: Could not write the suffix_tree to file!\n";
        return -1;
    }
    cout << "Suffix tree file " << treeFileName << " was successfully created.\n";

    delete global_ref;
    delete refdb;

    return 0;
}


int loadIndices(CharString fastaFileName) {
    CharString baseFileName = fastaFileName;
    while (back(baseFileName) != '.')
        eraseBack(baseFileName);
    eraseBack(baseFileName);
    string treeFileName = string(toCString(baseFileName)) + ".stree";

    if (!load_from_file(suffix_tree, treeFileName) )
    {
        cerr << "ERROR: Could not load the suffix_tree from file!\n";
        return -1;
    }
    cout << "Suffix tree file " << treeFileName << " was successfully loaded.\n";

    return 0;
}


int main(int argc, char const ** argv) {
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res == seqan::ArgumentParser::PARSE_OK) {
        if (options.build) {
            cout << "building index from " << options.fastaFileName << endl;
            buildIndices(options.fastaFileName);
        }
        else if (options.load) {
            cout << "loading index for " << options.fastaFileName << endl;
            loadIndices(options.fastaFileName);
        }

        if (options.process_ref) {
            cout << options.tree_dir << endl;
            iter_vec<TYPE> = new vector<TYPE> [2];
            profileReference(options.max_depth);
            delete [] iter_vec<TYPE>;
        }
    }
    else
        return -1;
    
    return 0;
}

