#include <sdsl/suffix_trees.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace std;
using namespace sdsl;
using namespace boost::filesystem;

struct Seed {
    unsigned frequency;
    unsigned beg;
    unsigned length;
    uint8_t ets;
};

class EtsSolver {
public:
	EtsSolver();
	~EtsSolver();
	void loadTree(string treeFileName, string ets_config, string ets_path);
	void setSkip(unsigned skip_length);
	void setSeedFreqThreshold(unsigned freq_threshold);
	void solveDNA(string& DNA, uint8_t required_ET);
    int getNumOfSeeds();
    int getTotalETS();
    int getTotalFreq();
	void printSeeds(ostream& stream);

private:
	//Internal functions
	//Load the first level (level[0])
	static bool compFreq(Seed left, Seed right);

	//Seed database, currently being a suffix tree	
	csa_wt<> tree;
    vector<unsigned> ets_profile;
    vector<uint8_t> *ets_distance;

    int skip_length;
	vector<Seed> seeds;
    int seed_count;
    int total_freq;
    unsigned seed_freq_threshold;
    unsigned total_ets;
};

