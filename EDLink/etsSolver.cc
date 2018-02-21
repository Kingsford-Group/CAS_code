#include "etsSolver.h"
#include <cstring>
#include <cassert>
#include <iostream>
#include <climits>
//#define DEBUG

EtsSolver::EtsSolver() {
    ets_distance = NULL;
	skip_length = 2;
    seed_freq_threshold = 10;
}

EtsSolver::~EtsSolver() {
    if (ets_distance != NULL)
        delete [] ets_distance;
}

void EtsSolver::loadTree(string treeFileName, string ets_config, string ets_path) {
	load_from_file(tree, treeFileName);

    std::ifstream profile_config_file(ets_config);
    unsigned profile_depth;
    profile_config_file >> profile_depth;
    while(!profile_config_file.eof() ) {
        ets_profile.push_back(profile_depth);
        profile_config_file >> profile_depth;
    }

    ets_distance = new vector<uint8_t> [ets_profile.size()];

    for (unsigned profile_idx = 0; profile_idx < ets_profile.size(); ++profile_idx) {
        cout << "loading profile: " << to_string(ets_profile[profile_idx]) + ".ets" << endl;

        boost::filesystem::ifstream profile_file(ets_path / (to_string(ets_profile[profile_idx]) + ".ets") );
        ets_distance[profile_idx].resize(tree.size() );
        unsigned distance;
        for (uint64_t loc = 0; loc < tree.size(); ++loc) {
            profile_file >> distance;
            ets_distance[profile_idx][loc] = distance;
        }
        profile_file.close();
    }
}

void EtsSolver::setSkip(unsigned skip_length) {
    this->skip_length = skip_length;
}

void EtsSolver::setSeedFreqThreshold(unsigned freq_threshold) {
    this->seed_freq_threshold = freq_threshold;
}

void EtsSolver::solveDNA(string& DNA, uint8_t required_ET) {
    seeds.clear();
    unsigned beg = 0;
    unsigned end = 1;
    unsigned length = 1000000;
    unsigned frequency = 0;
    uint8_t ets = 1;
    unsigned ets_idx = 0;
    while (end <= DNA.length() ) {
        while (count(tree, DNA.substr(beg, end - beg) ) != 0 && end <= DNA.length() ) {
            if (end - beg == ets_profile[ets_idx]) {
                length = end - beg;
                frequency = count(tree, DNA.substr(beg, length) );
                ets = ets_distance[ets_idx][locate(tree, DNA.substr(beg, length) )[0] ];
                ++ets_idx;
            }
            ++end;
        }
       
        seeds.push_back(Seed() );
        // In case this is the last seed
        if ( (end - beg) > length) {
            seeds.back().beg = beg;
            seeds.back().length = length;
            seeds.back().frequency = frequency;
            seeds.back().ets = ets;
        }
        else {
            seeds.back().beg = beg;
            seeds.back().length = end - 1 - beg;
            seeds.back().frequency = count(tree, DNA.substr(beg, end - 1 - beg) );
            seeds.back().ets = 1;
        }

        if (seeds.back().frequency > seed_freq_threshold)
            seeds.pop_back();

        // Reset for next seed
        end = end + skip_length;
        beg = end - 1;
        ets_idx = 0;
        ets = 1;
    }

	sort(seeds.begin(), seeds.end(), this->compFreq);

    total_ets = 0;
    total_freq = 0;
    for (seed_count = 0; seed_count < (int) seeds.size(); ++seed_count) {
        total_freq += seeds[seed_count].frequency;
        total_ets += seeds[seed_count].ets;
        
        if (total_ets >= required_ET) {
            ++seed_count;
            break;
        }
    }

    if (total_ets < required_ET) {
        total_freq = -1;
        seed_count = 0;
    }
}

int EtsSolver::getNumOfSeeds() {
    return seed_count;
}

int EtsSolver::getTotalFreq() {
    return total_freq;
}

int EtsSolver::getTotalETS() {
    return total_ets;
}

bool EtsSolver::compFreq(Seed left, Seed right) {
    if (left.frequency < right.frequency)
        return true;
    else if (left.frequency == right.frequency)
        return left.ets > right.ets;
    else
	    return false;
}

void EtsSolver::printSeeds(ostream& stream) {
    for (int seed_num = 0; seed_num < seed_count; ++seed_num) {
        stream << "Seed[" << seed_num << "]: start: " << seeds[seed_num].beg;
        stream << " length: " << seeds[seed_num].length;
        stream << " frequency: " << seeds[seed_num].frequency;
        stream << " ets: " << (unsigned) seeds[seed_num].ets << endl;
    }
}

