#include <sdsl/suffix_trees.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "RefDB.h"

template <typename T>
struct test_struct{
    T node;
    unsigned depth;
};

template<typename T>
vector<T> test;

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    cst_sada<> cst;
    csa_sada<> csa;
    csa_wt<> fm_index;
    RefDB refDB;

    if (argc == 3) {
        string ref_name = argv[1];
        refDB.loadRefFile(ref_name);
        string* global_ref = refDB.globalReference();
        construct_im(fm_index, *global_ref, 1);
        store_to_file(fm_index, argv[2]);

        //cout << locate(fm_index, "CCTGCCCCACAGCCTTGCCTGGATTTCTATCTCCCTGGCTTGGTGCCAG")[0] << endl;

        cout << "size: " << fm_index.size() << endl;
    }
    
    string ref = "CTTATTGATACTTACNNAAATAG|";
    construct_im(cst, ref.c_str(), 1);
    construct_im(csa, ref.c_str(), 1);
    construct_im(fm_index, ref.c_str(), 1);

    unsigned degree = cst.degree(cst.root() );
    cout << degree << endl;

    for (unsigned i = 1; i <= degree; ++i)
        cout << extract(cst, cst.select_child(cst.root(), i) ) << endl;

    //load_from_file(cst, "test_file.stree");
    cout << "inner nodes : " << cst.nodes() - cst.csa.size() << endl;
    //auto u = cst.select_child(cst.child(cst.root(), 'A'), 2);
    auto u = cst.child(cst.child(cst.root(), 'T'), 'T');
    auto v = cst.child(cst.child(cst.child(cst.root(), 'A'), 'C'), 'N');
    auto du = cst.depth(u);
    auto dv = cst.depth(v);

    cout << "u : " << du << "-[" << cst.lb(u) << "," << cst.rb(u) << "]" << endl;
    cout << "extract(cst, v) : " << extract(cst, u) << endl;
    cout << "degree: " << cst.degree(u) << endl;
    cout << "is_leaf:" << cst.is_leaf(u) << endl;
    cout << "leftmost_leaf: " << cst.leftmost_leaf(u) << endl;
    cout << "rightmost_leaf: " << cst.rightmost_leaf(u) << endl;
    
    for (uint64_t i = cst.lb(u); i <= cst.rb(u); ++i) {
        auto node = cst.select_leaf(i + 1);
        cout << "leaf[" << i << "]: " << extract(cst, node);
        cout << " at loc: " << cst.sn(node) << endl;
    }

    cout << "lb: " << cst.lb(u) << endl;
    cout << "rb: " << cst.rb(u) << endl;

    cout << "v : " << dv << "-[" << cst.lb(v) << "," << cst.rb(v) << "]" << endl;
    cout << "extract(cst, v) : " << extract(cst, v) << endl;
    cout << "last letter(cst, v) : " << cst.edge(v, cst.depth(v)) << endl;
    cout << "degree: " << cst.degree(v) << endl;
    cout << "is_leaf:" << cst.is_leaf(v) << endl;
    cout << "sa entry of v: " << cst.sn(v) << endl;

    string query = "TT";
    auto occs = locate(cst, query.begin(), query.end() );

    for (uint64_t i = 0; i < count(fm_index, query.begin(), query.end() ); ++i)
        cout << occs[i] << endl;

    /*
#define ITEM   test_struct<decltype(cst.root() )> 

    //vector<test_struct<decltype(v)> > test;
    test<ITEM>.push_back(test_struct<decltype(v)> () );

    test<ITEM>[0].node = u;
    test<ITEM>[0].depth = cst.depth(u);

    cout << "test: " << test<ITEM>[0].depth << " " << extract(cst, test<ITEM>[0].node) << endl;
    */
    store_to_file(cst, "test.stree");

    //delete global_ref;
}
