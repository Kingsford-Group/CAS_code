#include <sdsl/suffix_trees.hpp>
#include <sdsl/suffix_arrays.hpp>
#include "RefDB.h"
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace sdsl;

struct Cell {
	Cell() {
		start = -1;
		end = -1;
		isleaf = false;
		frequency = -1;
		lstart = -1;
		lend = -1;
		lfreq = -1;
		rstart = -1;
		rend = -1;
		rfreq = -1;
	};

	unsigned int start;
	unsigned int end;
	bool isleaf;
	unsigned int frequency;

	unsigned int lstart;
	unsigned int lend;
	unsigned int lfreq;
	unsigned int rstart;
	unsigned int rend;
	unsigned int rfreq;
};

class TreeWraper {
public:
    void loadTree(string file_name) {
        load_from_file(tree, file_name);
    };

    unsigned getHashLength() {return 10;};

    void generateTree(string file_name) {
        refdb.loadRefFile(file_name);
        string* global_ref = refdb.globalReference();
        construct_im(tree, *global_ref, 1);
        store_to_file(tree, file_name + ".stree");
        
        delete [] global_ref;
    };

    void setHashLength(unsigned) {};

    void query(string seed) {
        seed_frequency = count(tree, seed);
    };

    uint64_t frequency() {return seed_frequency;};
    
    bool isLeaf() {return seed_frequency <= 1;};

private:
    csa_wt<> tree;
    RefDB refdb;
    uint64_t seed_frequency;
};

class OptimalSolverLN {
public:
	OptimalSolverLN();
	~OptimalSolverLN();
	void loadTree(string treeFileName);
	void generateTree(string treeFileName);
	void init(int readLength, int seedNum);
	void reset();
	unsigned int solveDNA(string& DNA);
	void fillMatrix(string& DNA);
	unsigned int calculateLastDiv();
	unsigned int calcualteFreq();
	void printSeeds(ostream& stream);
	void printFreqs(ostream& stream);
	void printLength(ostream& stream);
	void printStats(ostream& stream);
	void sortOfFreq();
	void sortOfLength();
	
	void backtrack();
	
	//For debugging
	void setMinLength(int minLength);
	void feedL0();

private:
	//Internal functions
	//Load the first level (level[0])
	void loadL0(string& DNA);
	//Returns the location of the first optimal div in read
	int solveFirstOptimal(int opt_div, int pos, int l);
	//Sort the seeds
	template<class T>
	void sortSeeds(T relation);
	static bool compFreq(Cell left, Cell right);
	static bool compLength(Cell left, Cell right);

	//This is for debugging
	bool L0Loaded;

	//Seed database, currently being a suffix tree	
	TreeWraper tree;

	//Internal data structures
	int readLength;
	int seedNum;
	int minLength;

	Cell* matrix;
	Cell* base;
	vector<Cell*> level;
	vector<Cell*> level0;
	Cell* defaultMatrix;
	Cell* defaultBase;
	int matrixSize;
	int baseSize;

	vector<Cell> seeds;

	unsigned int finalDiv;

	//Statistics
	vector<unsigned long long> divTravel;
	unsigned int processedReads;
};
