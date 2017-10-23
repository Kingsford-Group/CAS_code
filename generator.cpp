#include <iostream>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/align.h>
#include <seqan/seq_io.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <deque>
#include <string>

using namespace seqan;

#define __UNIT_LEN__  12
#define __LARGE_NUM__  10000000000
#define __TOGGLE_THRESHOLD__  500
#define __HALT_THRESHOLD__  5
//#define __SKIP_REPEAT__

typedef FastFMIndexConfig<void, uint64_t> TConfig;
typedef Index<StringSet<String<Dna5> >, IndexQGram<UngappedShape<__UNIT_LEN__> > > KIndex;
typedef Index<StringSet<String<Dna5> >, BidirectionalIndex<FMIndex<void, TConfig> > > BWTIndex;
typedef Align<String<Dna5>, ArrayGaps> TAlign;
typedef typename Fibre<BWTIndex, FibreSA>::Type         TSA;
typedef typename Iterator<TSA const, Standard>::Type    TIterator;
typedef ModifiedString<String<Dna5>, ModComplementDna>   TMyComplement;
typedef ModifiedString<TMyComplement, ModReverse> TMyReverseComplement;

struct ETSConfig {
    int unit_len = __UNIT_LEN__;
    int max_units = -1;
    bool load_ets;

    std::vector<int> start_pos;

    ETSConfig() :
        load_ets(false) {}
};


struct CmdOptions {
    bool build;
    bool load;
    bool process_ref;
    bool benchmark;
    bool profile;
    bool interactive;
    int start_length;
    int end_length;
    int target_freq;
    int k_length;
    seqan::CharString fastaFileName;
    seqan::CharString etsFileName;
    seqan::CharString benchmarkFileName;

    CmdOptions() :
        interactive(false), profile(false), build(false), load(false), process_ref(false), benchmark(false) {}

};


struct LocStruct {
    long unsigned seq_no;
    long unsigned seq_offset;
    unsigned support_units;
}; 


struct LocFreqCmp {
    inline bool operator() (const LocStruct & struct1, const LocStruct & struct2) {
        return (struct1.support_units > struct2.support_units);
    }
};


struct LocPosCmp {
    inline bool operator() (const LocStruct & struct1, const LocStruct & struct2) {
        if (struct1.seq_no < struct2.seq_no)
            return true;
        else if (struct1.seq_no > struct2.seq_no)
            return false;
        else
            return (struct1.seq_offset < struct2.seq_offset);
    }
};


struct MultiFreqData {
    long unsigned relative_pos;
    long unsigned size;
};


ETSConfig globalConfig;
BWTIndex fmIndex;
//Finder<BWTIndex> fmFinder;
StringSet<String<Dna5> > genomeSet;
StringSet<String<unsigned> > flagSet;
KIndex kIndex;
//TAlign aligner;
std::vector<LocStruct> potential_locs;
std::vector<unsigned> loc_iters;
std::vector<Infix<const String<Pair<long unsigned, long unsigned, Pack> > >::Type> unit_locs;
std::vector<bool> unit_freq_switch;
std::vector<std::vector<LocStruct>::iterator > freq_unit_locs;
std::vector<unsigned> unit_loc_length;
Score<int, Simple> scoreConfig(0, -1, -1);
AlignConfig<false, true, true, false> alignConfig;
std::vector<LocStruct> multi_freq_data;
//TIterator fm_pos_reference;

std::vector<long unsigned> *streamlinedDir;
std::vector<long unsigned> *streamlinedPosDir;
std::vector<LocStruct> *streamlinedPos;
std::vector<LocStruct> *unitPos;

std::vector<int> base_freq_ary;
std::vector<int> major_freq_ary;

unsigned parseFlag(unsigned & flag, int unit_num) {
   unsigned ets_num = 1;
   while (ets_num < unit_num && !( (1 << ets_num) & globalConfig.start_pos[unit_num - 1]) )
       ++ets_num;

   return ets_num;
};

unsigned parseETS(Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > & bwt_iter, unsigned seed_length) {
    auto occ = getOccurrences(bwt_iter, Fwd() );
    unsigned long seq_no = getSeqNo(occ[0]);
    unsigned long seq_offset = getSeqNo(occ[0]);

    return parseFlag(flagSet[seq_no][seq_offset], seed_length / __UNIT_LEN__);
}

void processReadNaive(String<Dna5> seq, unsigned error_threshold, unsigned & freq_num) {
    unsigned seed_num = error_threshold + 1;
    unsigned seed_length = length(seq) / seed_num;
    freq_num = 0;

    for (int i = 0; i < 2; ++i) {

        for (int seed_iter = 0; seed_iter < seed_num; ++seed_iter) {
            Infix<String<Dna5> >::Type seed = infix(seq, seed_iter * seed_length, (seed_iter + 1) * seed_length);
            Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);
            goDown(bwt_iter, seed, Rev() );

            freq_num += countOccurrences(bwt_iter);
        }

        seq = Dna5StringReverseComplement(seq);
    }
}

void processRead(String<Dna5> seq, unsigned error_threshold, unsigned & seed_num, unsigned & error_tolerant_num, unsigned & freq_num) {
    unsigned seed_num_ary[2] = {0,  0};
    unsigned total_seed_num_ary[2] = {0,  0};
    unsigned ets_num_ary[2] = {0,  0};
    unsigned freq_ary[2] = {0,  0};

    for (int i = 0; i < 2; ++i) {
        Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);

        int bp_iter = length(seq) - 1;
        unsigned unit_length = 0;

        while (bp_iter >= 0) {
            if (goDown(bwt_iter, seq[bp_iter], Fwd() )  && unit_length < globalConfig.max_units * __UNIT_LEN__) {
                ++unit_length;
            }
            else if (unit_length != 0) {

                if (ets_num_ary[i] < error_threshold) {
                    unsigned seed_freq = countOccurrences(bwt_iter);
                    unsigned seed_ets = 1;

                    if (unit_length >= 2 * __UNIT_LEN__) {
                        seed_ets = parseETS(bwt_iter, unit_length);
                    }

                    ets_num_ary[i] += seed_ets;
                    ++seed_num_ary[i];
                    freq_ary[i] += seed_freq;

                //std::cout << representative(bwt_iter) << ": " << seed_freq << " " << seed_ets << std::endl;
                }

                total_seed_num_ary[i]++;

                // reset iterator
                goRoot(bwt_iter);

                if (goDown(bwt_iter, seq[bp_iter], Fwd() ))
                    unit_length = 1;
                else {
                    goRoot(bwt_iter);
                    unit_length = 0;
                }
            }
                
        

            --bp_iter;
        }

        if (unit_length != 0) {

            if (ets_num_ary[i] < error_threshold) {
                unsigned seed_freq = countOccurrences(bwt_iter);
                unsigned seed_ets = 1;

                if (unit_length >= 2 * __UNIT_LEN__) {
                    seed_ets = parseETS(bwt_iter, unit_length);
                }

                ets_num_ary[i] += seed_ets;
                ++seed_num_ary[i];
                
                //std::cout << representative(bwt_iter) << ": " << seed_freq << " " << seed_ets << std::endl;
            }

            total_seed_num_ary[i]++;
        }

        //std::cout << seq << std::endl;
        //std::cout << total_seed_num_ary[i] << " " <<seed_num_ary[i] << " " << ets_num_ary[i] << " " << freq_ary[i] << std::endl;
        seq = Dna5StringReverseComplement(seq);
    }

    unsigned return_idx = 0;
    if (total_seed_num_ary[0] < total_seed_num_ary[1])
        return_idx = 0;
    else
        return_idx = 1;

    seed_num = seed_num_ary[return_idx];
    error_tolerant_num = ets_num_ary[return_idx];
    freq_num = freq_ary[return_idx];
}

void storeStreamline (std::string filename) {
    long unsigned last_dir = streamlinedDir->at(length(indexDir(kIndex) ) - 1);
    long unsigned last_pos = streamlinedPosDir->at(last_dir);

    long unsigned dir_size = last_dir + 1;
    long unsigned pos_size = last_pos + 1;

    std::ofstream output_file;
    output_file.open(filename.c_str(), std::ios::out | std::ios::binary);

    output_file.write((char*) &(dir_size), sizeof(long unsigned) );
    output_file.write((char*) &(pos_size), sizeof(long unsigned) );

    for (long unsigned hash = 0; hash < length(indexDir(kIndex) ); ++hash)
        output_file.write((char*) &(streamlinedDir->at(hash) ), sizeof(long unsigned) );
    
    for (long unsigned dir_idx = 0; dir_idx <= last_dir; ++dir_idx)
        output_file.write((char*) &(streamlinedPosDir->at(dir_idx) ), sizeof(long unsigned) );

    for (long unsigned pos_idx = 0; pos_idx <= last_pos; ++pos_idx) {
        output_file.write((char*) &(streamlinedPos->at(pos_idx).seq_no), sizeof(long unsigned) );
        output_file.write((char*) &(streamlinedPos->at(pos_idx).seq_offset), sizeof(long unsigned) );
        output_file.write((char*) &(streamlinedPos->at(pos_idx).support_units), sizeof(unsigned) );
    }

    output_file.close();
}

void loadStreamline (std::string filename) {
    std::ifstream input_file;
    input_file.open(filename.c_str(), std::ios::in | std::ios::binary);
    
    long unsigned dir_size;
    long unsigned pos_size;

    input_file.read((char*) &dir_size, sizeof(unsigned long) );
    input_file.read((char*) &pos_size, sizeof(unsigned long) );

    streamlinedDir->resize(length(indexDir(kIndex) ) );
    streamlinedPosDir->resize(dir_size);
    streamlinedPos->resize(pos_size);

    for (long unsigned hash = 0; hash < length(indexDir(kIndex) ); ++hash)
        input_file.read((char*) &(streamlinedDir->at(hash) ), sizeof(long unsigned) );
    
    for (long unsigned dir_idx = 0; dir_idx < dir_size; ++dir_idx)
        input_file.read((char*) &(streamlinedPosDir->at(dir_idx) ), sizeof(long unsigned) );

    for (long unsigned pos_idx = 0; pos_idx < pos_size; ++pos_idx) {
        input_file.read((char*) &(streamlinedPos->at(pos_idx).seq_no), sizeof(long unsigned) );
        input_file.read((char*) &(streamlinedPos->at(pos_idx).seq_offset), sizeof(long unsigned) );
        input_file.read((char*) &(streamlinedPos->at(pos_idx).support_units), sizeof(unsigned) );
    }

    input_file.close();
}

unsigned streamlineUnit(Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter, unsigned left_extend, unsigned right_extend, unsigned streamlined_loc_num, unsigned unit_num, unsigned ori_left_extend) {
    auto locs = getOccurrences(bwt_iter, Fwd());
    unsigned count_freq = length(locs);
    unsigned next_loc_num;
    unsigned left_offset = ori_left_extend - left_extend;

    //std::cout << representative(bwt_iter) << std::endl;

    if (right_extend == 0 && left_extend == 0) {
        next_loc_num = streamlined_loc_num + 1;

        if (unitPos->size() < next_loc_num)
            unitPos->resize(1.5 * unitPos->size() );

        unitPos->at(streamlined_loc_num).seq_no = getSeqNo(locs[0]);
        unitPos->at(streamlined_loc_num).seq_offset = getSeqOffset(locs[0]) + left_offset;
        unitPos->at(streamlined_loc_num).support_units = unit_num;

        return next_loc_num;
    }

    if (count_freq <= __HALT_THRESHOLD__) {
        next_loc_num = streamlined_loc_num + count_freq;

        if (unitPos->size() < next_loc_num)
            unitPos->resize(1.5 * unitPos->size() );

        for (unsigned i = 0; i < count_freq; ++i) {
            unitPos->at(streamlined_loc_num + i).seq_no = getSeqNo(locs[i]);
            unitPos->at(streamlined_loc_num + i).seq_offset = getSeqOffset(locs[i]) + left_offset;
            unitPos->at(streamlined_loc_num + i).support_units = 0;
        }

        return next_loc_num;
    }

    //std::cout << representative(bwt_iter) << std::endl;

    if (right_extend != 0) {
        if(!goDown(bwt_iter, Rev() ) ) {
            streamlined_loc_num = streamlineUnit(bwt_iter, left_extend, 0, streamlined_loc_num, unit_num, ori_left_extend);
        }
        else {
            streamlined_loc_num = streamlineUnit(bwt_iter, left_extend, right_extend - 1, streamlined_loc_num, unit_num, ori_left_extend);
            while (goRight(bwt_iter, Rev() ) )
                streamlined_loc_num = streamlineUnit(bwt_iter, left_extend, right_extend - 1, streamlined_loc_num, unit_num, ori_left_extend);
        }
    }
    else {
        if (!goDown(bwt_iter, Fwd() ) ) {
            streamlined_loc_num = streamlineUnit(bwt_iter, 0, 0, streamlined_loc_num, unit_num, ori_left_extend - left_extend);
        }
        else {
            streamlined_loc_num = streamlineUnit(bwt_iter, left_extend - 1, 0, streamlined_loc_num, unit_num, ori_left_extend);
            while (goRight(bwt_iter, Fwd() ) )
                streamlined_loc_num = streamlineUnit(bwt_iter, left_extend - 1, 0, streamlined_loc_num, unit_num, ori_left_extend);
        }
    }
    
    return streamlined_loc_num;
}

void streamlineKIndex(StringSet<String<Dna5> > & ref, unsigned unit_num) {
    unsigned temp = length(indexDir(kIndex) );
    streamlinedDir->resize(length(indexDir(kIndex) ) );
    streamlinedPosDir->resize((long unsigned) (100 + 0.001 * length(ref) ) );
    streamlinedPos->resize((long unsigned) (1000 + 0.1 * length(ref) ) );
    unitPos->resize(1000);

    long unsigned pos_dir_idx = 0;
    long unsigned pos_idx = 0;

    for (long unsigned hash = 0; hash < length(indexDir(kIndex) ) - 1; ++hash) {
        streamlinedDir->at(hash) = pos_dir_idx;

        if (indexDir(kIndex)[hash + 1] - indexDir(kIndex)[hash] >= __TOGGLE_THRESHOLD__) {
            String<Dna5> unit_string;
            unhash(unit_string, hash, __UNIT_LEN__);

            Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);
            goDown(bwt_iter, unit_string, Rev() );

            std::cout << unit_string << ": " << indexDir(kIndex)[hash + 1] - indexDir(kIndex)[hash] << std::endl;

            for (int unit_iter = 0; unit_iter < unit_num; ++unit_iter) {
                streamlinedPosDir->at(pos_dir_idx) = pos_idx;

                unsigned streamlined_loc_num = 0;
                streamlined_loc_num = streamlineUnit(bwt_iter, unit_num + unit_iter * __UNIT_LEN__, unit_num + (unit_num - 1 - unit_iter) * __UNIT_LEN__, streamlined_loc_num, unit_num, unit_num + unit_iter * __UNIT_LEN__);

                //std::cout << streamlined_loc_num << " ";

                sort(unitPos->begin(), unitPos->begin() + streamlined_loc_num, LocPosCmp() );

                unsigned long next_pos_dir_idx = pos_dir_idx + 1;
                unsigned long next_pos_idx = pos_idx + streamlined_loc_num;

                if (streamlinedPosDir->size() <= next_pos_idx + 1)
                    streamlinedPosDir->resize(1.5 * (next_pos_idx + 1) );

                if (streamlinedPos->size() <= next_pos_idx + 1)
                    streamlinedPos->resize(1.5 * (next_pos_idx + 1) );

                std::copy(unitPos->begin(), unitPos->begin() + streamlined_loc_num, streamlinedPos->begin() + pos_idx);

                pos_dir_idx = next_pos_dir_idx;
                pos_idx = next_pos_idx;
            }

            //std::cout << std::endl;
        }
    }

    streamlinedDir->at(length(indexDir(kIndex) ) - 1) = pos_dir_idx;
    streamlinedPosDir->at(pos_dir_idx) = pos_idx;
}


unsigned checkSeqs(Infix<String<Dna5> >::Type & seq, Infix<String<Dna5> >::Type & ref, unsigned error_limit) {
    //assignSource(row(aligner, 0), seq);
    //assignSource(row(aligner, 1), ref);

    //std::cout << seq << std::endl << ref << std::endl;

    //int score = -globalAlignment(aligner, scoreConfig, alignConfig);
    int score = -globalAlignmentScore(seq, ref, scoreConfig, alignConfig, -(int)(2 * error_limit + 1), error_limit + 1);

    //std::cout << score << std::endl;

    if (score < error_limit)
        return score;
    else
        return error_limit;
}


unsigned findFirstFailErrors(Infix<String<Dna5> >::Type & seq, unsigned unit_num, std::vector<LocStruct> & origin_pos_ary, unsigned origin_ary_len) {
    unsigned unit_length = globalConfig.unit_len;
    unsigned total_potential_locs = 0;
    unsigned error_limit = unit_num - 1;
    
    //std::cout << "origin_ary_len: " << origin_ary_len << std::endl;

    assert(length(seq) == unit_num * unit_length);

    for (int unit_iter = 0; unit_iter < unit_num; unit_iter++) {
        Infix<String<Dna5> >::Type unit_seq = infix(seq, unit_iter * unit_length, (unit_iter + 1) * unit_length);

        hash(indexShape(kIndex), begin(unit_seq) );
        unit_locs[unit_iter] = getOccurrences(kIndex, indexShape(kIndex) );
        unit_loc_length[unit_iter] = length(unit_locs[unit_iter]);
        unit_freq_switch[unit_iter] = unit_loc_length[unit_iter] >= __TOGGLE_THRESHOLD__;
        //unit_freq_switch[unit_iter] = false;

        // toggle ultra frequent units
        if (unit_freq_switch[unit_iter]) {
            unsigned long pos_idx = streamlinedPosDir->at(streamlinedDir->at(value(indexShape(kIndex) ) ) + unit_iter);
            freq_unit_locs[unit_iter] = streamlinedPos->begin() + pos_idx;
            unit_loc_length[unit_iter] = streamlinedPosDir->at(streamlinedDir->at(value(indexShape(kIndex) ) ) + unit_iter + 1) - pos_idx;
        }

        //for (int i = 0; i < length(unit_locs[unit_iter]); i++)
            //std::cout << unit_locs[unit_iter][i] << std::endl;

        loc_iters[unit_iter] = 0;

        total_potential_locs += length(unit_locs[unit_iter]);
    }

    if (potential_locs.size() < total_potential_locs)
        potential_locs.resize(total_potential_locs);

    unsigned potential_loc_len = 0;
    bool loc_finished = false;
    long unsigned min_seq_no;
    long unsigned min_seq_offset;
    long unsigned origin_seq_no;
    long unsigned origin_seq_offset;
    unsigned origin_iter = 0;

    while (!loc_finished) {
        loc_finished = true;
        min_seq_no = __LARGE_NUM__;
        min_seq_offset = __LARGE_NUM__;
        
        unsigned support_units = 0;

        for (int unit_iter = 0; unit_iter < unit_num; ++unit_iter) {

            if (loc_iters[unit_iter] < unit_loc_length[unit_iter]) {
                unsigned dist_to_beg = unit_iter * unit_length;
                long unsigned unit_seq_no;
                long unsigned unit_seq_offset;

                unsigned cur_support_units = 0;

                if (unit_freq_switch[unit_iter]) {
                    unit_seq_no = freq_unit_locs[unit_iter][loc_iters[unit_iter]].seq_no;
                    unit_seq_offset = freq_unit_locs[unit_iter][loc_iters[unit_iter]].seq_offset;

                    if (freq_unit_locs[unit_iter][loc_iters[unit_iter]].support_units != 0) {
                        cur_support_units = unit_num;

                       //std::cout << unit_iter << " " << cur_support_units << std::endl;
                    }
                }
                else {
                    unit_seq_no = getSeqNo(unit_locs[unit_iter][loc_iters[unit_iter] ]);
                    unit_seq_offset = getSeqOffset(unit_locs[unit_iter][loc_iters[unit_iter] ]);
                }

                long unsigned corrected_seq_offset = (unit_seq_offset > dist_to_beg) ? unit_seq_offset - dist_to_beg : 0;

                if (unit_seq_no < min_seq_no) {
                    min_seq_no = unit_seq_no; 
                    min_seq_offset = corrected_seq_offset; 
                    support_units = cur_support_units;
                }
                else if (unit_seq_no == min_seq_no) {

                    if (corrected_seq_offset < min_seq_offset) {
                        min_seq_offset = corrected_seq_offset;
                        support_units = cur_support_units;
                    }
                }
            }
        }
       
        bool ignore = false;

        while (origin_iter < origin_ary_len) {
            origin_seq_no = origin_pos_ary[origin_iter].seq_no;
            origin_seq_offset = origin_pos_ary[origin_iter].seq_offset;

           //std::cout << "origin_seq: " << origin_seq_no << "-" << origin_seq_offset << std::endl;

            if (origin_seq_no < min_seq_no || (origin_seq_no == min_seq_no && origin_seq_offset + error_limit < min_seq_offset) ) {
                ++origin_iter;
            }
            else if (origin_seq_no == min_seq_no && origin_seq_offset <= min_seq_offset + error_limit) {
                ++origin_iter;
                ignore = true;
                break;
            }
            else {
                //++origin_iter;
                break;
            }
        }
        
       //std::cout << "min_seq: " << min_seq_no << "x" << min_seq_offset << std::endl;
       //std::cout << "support_units: " << support_units << std::endl;
       //std::cout << "ignore: " << ignore << std::endl;


        for (int unit_iter = 0; unit_iter < unit_num; ++unit_iter) {

            if (loc_iters[unit_iter] < unit_loc_length[unit_iter]) {
                long unsigned unit_seq_no;
                long unsigned unit_seq_offset;

                if (unit_freq_switch[unit_iter]) {
                    unit_seq_no = freq_unit_locs[unit_iter][loc_iters[unit_iter]].seq_no;
                    unit_seq_offset = freq_unit_locs[unit_iter][loc_iters[unit_iter]].seq_offset;
                }
                else {
                    unit_seq_no = getSeqNo(unit_locs[unit_iter][loc_iters[unit_iter] ]);
                    unit_seq_offset = getSeqOffset(unit_locs[unit_iter][loc_iters[unit_iter] ]);
                }

                if (unit_seq_no == min_seq_no) {
                    unsigned dist_to_beg = unit_iter * unit_length;
                    long unsigned corrected_seq_offset = (unit_seq_offset > dist_to_beg) ? unit_seq_offset - dist_to_beg : 0;

                    if (corrected_seq_offset <= min_seq_offset + error_limit) {

                        if (support_units < unit_num)
                            ++support_units;

                        ++loc_iters[unit_iter];

                        // move the iterator beyond the potential loc
                        while (loc_iters[unit_iter] < unit_loc_length[unit_iter] ) {

                            if (unit_freq_switch[unit_iter]) {
                                unit_seq_no = freq_unit_locs[unit_iter][loc_iters[unit_iter]].seq_no;
                                unit_seq_offset = freq_unit_locs[unit_iter][loc_iters[unit_iter]].seq_offset;
                            }
                            else {
                                unit_seq_no = getSeqNo(unit_locs[unit_iter][loc_iters[unit_iter] ]);
                                unit_seq_offset = getSeqOffset(unit_locs[unit_iter][loc_iters[unit_iter] ]);
                            }

                            if (unit_seq_no != min_seq_no)
                                break;
                            else if (unit_seq_offset > min_seq_offset + error_limit + dist_to_beg)
                                break;

                            ++loc_iters[unit_iter];
                        }

                        if (loc_iters[unit_iter] < unit_loc_length[unit_iter] )
                            loc_finished = false;
                    }
                    else
                        loc_finished =false;
                }
                else 
                    loc_finished = false;
            }
        }

        potential_locs[potential_loc_len].seq_no = min_seq_no;
        potential_locs[potential_loc_len].seq_offset = min_seq_offset;
        potential_locs[potential_loc_len].support_units = support_units;

        if (!ignore)
            ++potential_loc_len;
    }

    std::sort(potential_locs.begin(), potential_locs.begin() + potential_loc_len, LocFreqCmp());

    unsigned tolerated_errors = unit_num;

    unsigned loc_iter;

    for (loc_iter = 0; loc_iter < potential_loc_len; ++loc_iter) {

       //std::cout << "loc iter: " << loc_iter << " support units: " << potential_locs[loc_iter].support_units << std::endl;
       //std::cout << potential_locs[loc_iter].seq_no << "x" << potential_locs[loc_iter].seq_offset << std::endl;

        if (unit_num - potential_locs[loc_iter].support_units >=  tolerated_errors)
            break;

       //std::cout << "no break" << std::endl;

//        if (potential_locs[loc_iter].seq_no != origin_seq_no || 
//                ( (potential_locs[loc_iter].seq_no == origin_seq_no) && 
 //                 (potential_locs[loc_iter].seq_offset + error_limit < origin_seq_offset || potential_locs[loc_iter].seq_offset > origin_seq_offset + error_limit) ) ) {
        long unsigned seq_offset = (potential_locs[loc_iter].seq_offset >= error_limit) ? potential_locs[loc_iter].seq_offset - error_limit : 0;
        
        //std::cout << potential_locs[loc_iter].seq_no << ":" << seq_offset << std::endl;

        Infix<String<Dna5> >::Type ref_seq = infix(genomeSet[potential_locs[loc_iter].seq_no], seq_offset, seq_offset + unit_num * unit_length + 2 * error_limit);

        tolerated_errors = checkSeqs(seq, ref_seq, tolerated_errors);
       
       //std::cout << "tolerated_errors: " << tolerated_errors << std::endl;
//        }
    }

    //std::cout << "potential_loc_len: " << potential_loc_len << std::endl;

    /*
    while (loc_iter < potential_loc_len) {
        //std::cout << "loc iter: " << loc_iter << std::endl;
        //std::cout << potential_locs[loc_iter].seq_no << "-" << potential_locs[loc_iter].seq_offset << std::endl;
        ++loc_iter;
    }
    */

    return tolerated_errors;
}


bool checkFlagProcessed(unsigned flag, int unit_num) {
    return flag & 1 << globalConfig.start_pos[unit_num - 1];
}


void updateFlag(unsigned & flag, int unit_num, bool multi_freq, int first_fail_errors) {
    assert(max_toerable_errors <= unit_num);

    if (multi_freq)
        flag = flag | 1 << globalConfig.start_pos[unit_num - 1];
    if (first_fail_errors != unit_num)
        flag = flag | (1 << globalConfig.start_pos[unit_num - 1] + first_fail_errors);
};


void processLoc(unsigned chromo_num, unsigned chromo_pos) {
    // Initialize counter
    int unit_iter = globalConfig.max_units;

    int first_fail_errors;
    
    // Loop until BWT fails
    while (unit_iter >= 2) {
        // check length
        
        //std::cout << "length: " << length(genomeSet[chromo_num]) << std::endl;
        
        if (chromo_pos + unit_iter * globalConfig.unit_len > length(genomeSet[chromo_num]) ) {
            first_fail_errors = 0;
            //updateFlag(flagSet[chromo_num][chromo_pos], unit_iter, false, first_fail_errors);
        }
        else {
            if (!checkFlagProcessed(flagSet[chromo_num][chromo_pos], unit_iter) ) {

                Infix<String<Dna5> >::Type seq = infix(genomeSet[chromo_num], chromo_pos, chromo_pos + unit_iter * globalConfig.unit_len);
                
    /*
                Iterator<BWTIndex, TopDown<> >::Type bwt_iter(fmIndex);
                goDown(bwt_iter, seq);

                if (length(getOccurrences(bwt_iter)) > 1) {
                goBegin(fmFinder);
                clear(fmFinder);
                find(fmFinder, seq);
    */

                Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);
                goDown(bwt_iter, seq, Rev() );

                bool multi_freq = false;
                unsigned freq_ary_len = 0;
                unsigned count_occ = countOccurrences(bwt_iter);
        
                if (count_occ > 1) {
                    multi_freq = true;

#ifdef __SKIP_REPEAT__
                    updateFlag(flagSet[chromo_num][chromo_pos], unit_iter,  multi_freq, unit_iter);
                    --unit_iter;
                    break;
#else

                    // for no-repeat case
                    // for no-repeat case

                    // check for resize
                    if (multi_freq_data.size() < count_occ)
                        multi_freq_data.resize(count_occ);

                    auto occ_ary = getOccurrences(bwt_iter, Fwd() );

                    for (unsigned occ_iter = 0; occ_iter < count_occ; ++occ_iter) {
                        multi_freq_data[freq_ary_len].seq_no = getSeqNo(occ_ary[occ_iter]);
                        multi_freq_data[freq_ary_len].seq_offset = getSeqOffset(occ_ary[occ_iter]);
                        ++freq_ary_len;
                    }

                    // Then sort
                    std::sort(multi_freq_data.begin(), multi_freq_data.begin() + freq_ary_len, LocPosCmp());
#endif
                }
                else {
                    freq_ary_len = 1;
                    multi_freq_data[0].seq_no = chromo_num;
                    multi_freq_data[0].seq_offset = chromo_pos;
                }

               //std::cout << "----------------" << std::endl;

                // process seq
                first_fail_errors = findFirstFailErrors(seq, unit_iter, multi_freq_data, freq_ary_len);

                // Finally update all flags
                for (unsigned freq_iter = 0; freq_iter < freq_ary_len; ++freq_iter) {
                    long unsigned seq_no = multi_freq_data[freq_iter].seq_no;
                    long unsigned seq_offset = multi_freq_data[freq_iter].seq_offset;

                    //std::cout << "updateing flag for " << seq_no << " " << seq_offset <<  " with " << unit_iter << " units. First error: " << first_fail_errors << std::endl;
                    updateFlag(flagSet[seq_no][seq_offset], unit_iter,  multi_freq, first_fail_errors);
                }
            }
        }

        --unit_iter;
    };

#ifdef __SKIP_REPEAT__
    // Finish the rest
    while (unit_iter >= 2) {
        first_fail_errors = 0;
        updateFlag(flagSet[chromo_num][chromo_pos], unit_iter, true, first_fail_errors);
        --unit_iter;
    }
#endif

    // Process single unit
    if (chromo_pos + globalConfig.unit_len > length(genomeSet[chromo_num]) ) {
        updateFlag(flagSet[chromo_num][chromo_pos], 1, true, 1);
    }
    else {
        Infix<String<Dna5> >::Type unit_seq = infix(genomeSet[chromo_num], chromo_pos, chromo_pos + globalConfig.unit_len);
        hash(indexShape(kIndex), begin(unit_seq) );

        if (length(getOccurrences(kIndex, indexShape(kIndex))) > 1)
            updateFlag(flagSet[chromo_num][chromo_pos], 1, true, 1);
        else
            updateFlag(flagSet[chromo_num][chromo_pos], 1, false, 1);
    }
}


void showExtensionFrequencies(String<Dna5> query, unsigned ext_length) {
    //std::cout << "query: " << query << std::endl;

    unsigned query_length = length(query);

    Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);

    goDown(bwt_iter, query, Rev() );

    unsigned rep_len = repLength(_iter(bwt_iter, Rev() ) );

    unsigned total_patterns = 0;
    unsigned total_freq = 0;

    do {
        if (rep_len == ext_length + query_length) {
            // Go process the repeat extension check
            auto iter_range = range(_iter(bwt_iter, Rev() ) );
            unsigned base_freq = iter_range.i2 - iter_range.i1;

            std::cout << representative(bwt_iter) << std::endl;
            std::cout << base_freq << " "; // << std::endl;

            auto occ = getOccurrences(bwt_iter)[0];
            std::cout << getSeqNo(occ) << " " << getSeqOffset(occ) << std::endl;

            // Navigate right or up...
            while (!goRight(bwt_iter, Rev() ) )
                if (!goUp(bwt_iter) )
                    break;

            ++total_patterns;
            total_freq += base_freq;
        }
        //else if (rep_len < ext_length + query_length) {
        else {
            if (!goDown(bwt_iter, Rev() ) && !goRight(bwt_iter, Rev() ) )
                while (goUp(bwt_iter) && !goRight(bwt_iter, Rev() ) )
                    ;
        }

        rep_len = repLength(_iter(bwt_iter, Rev() ) );

        //std::cout << "rep_len: " << rep_len << std::endl;
    } 
    while (rep_len > query_length);

    std::cout << "total_patterns: " << total_patterns << std::endl;
    std::cout << "total_freq: " << total_freq << std::endl;
}


int profileRepeatSequence(Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter, unsigned end_length, unsigned base_freq) {
    decltype (range(_iter(bwt_iter, Rev() ) ) ) iter_range;

    bool keep_going;
    int iter_freq;

    while (repLength(_iter(bwt_iter, Rev() ) ) != end_length) {
        if (!goDown(bwt_iter, Rev() ) ) {
            return -1;
        }

        keep_going = false; 

        do {
            iter_range = range(_iter(bwt_iter, Rev()));
            iter_freq = iter_range.i2 - iter_range.i1;

            if (iter_freq > base_freq * 0.5) {
                keep_going = true;
                break;
            }

        } while (goRight(bwt_iter, Rev() ) );

        if (!keep_going)
            return -1;

    }

    return (iter_range.i2 - iter_range.i1);
}


void profileRepeatGenome(unsigned start_length, unsigned end_length, unsigned target_freq) {
    Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);

    do {
        unsigned rep_len = repLength(_iter(bwt_iter, Rev() ) );

        if (rep_len == start_length) {
            std::cout << representative(bwt_iter) << std::endl;

			/*
            // Go process the repeat extension check
            auto iter_range = range(_iter(bwt_iter, Rev() ) );
            unsigned base_freq = iter_range.i2 - iter_range.i1;

            if (base_freq >= target_freq) {
                int major_freq = profileRepeatSequence(bwt_iter, end_length, base_freq);

                base_freq_ary.push_back(base_freq);
                major_freq_ary.push_back(major_freq);

                std::cout << representative(bwt_iter) << std::endl;
                std::cout << base_freq << std::endl;
                std::cout << major_freq << std::endl;
            }
			*/

            // Navigate right or up...
            while (!goRight(bwt_iter, Rev() ) )
                if (!goUp(bwt_iter) )
                    break;
                
        }
        else if (rep_len < start_length) {
            if (!goDown(bwt_iter, Rev() ) && !goRight(bwt_iter, Rev() ) )
                while (goUp(bwt_iter) && !goRight(bwt_iter, Rev() ) )
                    ;
        }
    } 
    while (!isRoot(bwt_iter));
}


seqan::ArgumentParser::ParseResult
parseCommandLine(CmdOptions & options, int argc, char const ** argv) {
    seqan::ArgumentParser parser("modify_string");

    addOption(parser, seqan::ArgParseOption("b", "fasta-file", "Build index from the fasta file.", seqan::ArgParseArgument::INPUT_FILE) );
    setValidValues(parser, "fasta-file", "fa fasta");

    addOption(parser, seqan::ArgParseOption("l", "index-file", "Load index files.", seqan::ArgParseArgument::INPUT_FILE) );
    setValidValues(parser, "index-file", "");

    addOption(parser, seqan::ArgParseOption("p", "process-ref", "Process the reference genome.") );

    addOption(parser, seqan::ArgParseOption("u", "unit-num", "Please provide the maximum unit Num. Default: 5", seqan::ArgParseArgument::INTEGER) );
    
    addOption(parser, seqan::ArgParseOption("L", "load-ets", "Please provide the error tolerant seed file", seqan::ArgParseArgument::INPUT_FILE) );
    
    addOption(parser, seqan::ArgParseOption("m", "benchmark", "Please provide the fastq file", seqan::ArgParseArgument::INPUT_FILE) );

    addOption(parser, seqan::ArgParseOption("r", "profile", "Please provide the begin_len, end_len, target_freq", seqan::ArgParseArgument::INTEGER, "INT", true, 3) );
    
    addOption(parser, seqan::ArgParseOption("i", "interactive", "Entering the interactive session") );

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

    globalConfig.load_ets = isSet(parser, "load-ets");
    
    if (globalConfig.load_ets)
        getOptionValue(options.etsFileName, parser, "load-ets");

    if (isSet(parser, "unit-num"))
        getOptionValue(globalConfig.max_units, parser, "unit-num");
    else
        globalConfig.max_units = 5;

    options.benchmark = isSet(parser, "benchmark");

    if (options.benchmark)
        getOptionValue(options.benchmarkFileName, parser, "benchmark");

    options.profile = isSet(parser, "profile");

    if (options.profile) {
        options.start_length = std::stoi(getOptionValues(parser, "profile")[0]);
        options.end_length = std::stoi(getOptionValues(parser, "profile")[1]);
        options.target_freq = std::stoi(getOptionValues(parser, "profile")[2]);
    }

    options.interactive = isSet(parser, "interactive");

    globalConfig.start_pos.resize(globalConfig.max_units);
    loc_iters.resize(globalConfig.max_units);
    unit_locs.resize(globalConfig.max_units);
    unit_freq_switch.resize(globalConfig.max_units);
    freq_unit_locs.resize(globalConfig.max_units);
    unit_loc_length.resize(globalConfig.max_units);

    int total_position = 0;

    for (int i = 0; i < globalConfig.max_units; i++) {
        globalConfig.start_pos[i] = total_position;
        total_position += (i + 1);
    }

    assert (total_position < 32);

    return seqan::ArgumentParser::PARSE_OK;
}


int buildIndices(CharString fastaFileName) {

    CharString baseFileName = fastaFileName;
    while (back(baseFileName) != '.')
        eraseBack(baseFileName);
    eraseBack(baseFileName);

    CharString faiFileName = baseFileName;
    append(faiFileName, ".fai");

    CharString fmFileName = baseFileName;
    append(fmFileName, ".fm");

    CharString kFileName = baseFileName;
    append(kFileName, ".kmer");

    CharString streamlinedFileName = baseFileName;
    append(streamlinedFileName, ".st_kmer");

    FaiIndex faiIndex;
    if (!build(faiIndex, toCString(fastaFileName) ) ) {
        std::cerr << "ERROR: Could not build FAI index for file " << fastaFileName << "." << std::endl;
        return -1;
    }

    if (!save(faiIndex, toCString(faiFileName)))
    {
        std::cerr << "ERROR: Could not write the Fai index to file!\n";
        return -1;
    }
    std::cout << "Fai Index file " << faiFileName << " was successfully created.\n";

    std::cout << "Number of indeces: " << numSeqs(faiIndex) << std::endl;

    int chromoNum = numSeqs(faiIndex);

    for (int seqIter = 0; seqIter < chromoNum; ++seqIter) {
        String<Dna5> chromosome;
        readSequence(chromosome, faiIndex, seqIter);
        appendValue(genomeSet, chromosome);
    }

    BWTIndex temp(genomeSet);
    fmIndex = temp;

    Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);
    goDown(bwt_iter, String<Dna5>("CGGCTTGGTGCTCACGCACACAGGAAAGTCCTT"), Rev() );
    std::cout << representative(bwt_iter) << std::endl;

    indexCreate(fmIndex, FibreSA() );
    if (!save(fmIndex, toCString(fmFileName)))
    {
        std::cerr << "ERROR: Could not write the FM index to file!\n";
        return -1;
    }
    std::cout << "FM Index file " << fmFileName << " was successfully created.\n";


    kIndex = KIndex(genomeSet);
    hash(indexShape(kIndex), begin(String<Dna5>("CGGCTTGGTGCTCACGCACACAGGAAAGTCCTT") ) );
    std::cout << length(getOccurrences(kIndex, indexShape(kIndex) ) ) << std::endl;

    if (!save(kIndex, toCString(kFileName)))
    {
        std::cerr << "ERROR: Could not write the K-mer index to file!\n";
        return -1;
    }
    std::cout << "K-mer Index file " << kFileName << " was successfully created.\n";

    std::cout << length(indexDir(kIndex) ) << std::endl;

    //streamlineKIndex(genomeSet, globalConfig.max_units);
    //storeStreamline(toCString(streamlinedFileName) );

/*
    */

    return 0;
}


int loadIndices(CharString fastaFileName, StringSet<String<Dna5> > & genomeSet, CharString etsFileName) {
    CharString baseFileName = fastaFileName;
    while (back(baseFileName) != '.')
        eraseBack(baseFileName);
    eraseBack(baseFileName);

    CharString faiFileName = baseFileName;
    append(faiFileName, ".fai");

    CharString fmFileName = baseFileName;
    append(fmFileName, ".fm");

    CharString kFileName = baseFileName;
    append(kFileName, ".kmer");

    CharString streamlinedFileName = baseFileName;
    append(streamlinedFileName, ".st_kmer");

    FaiIndex faiIndex;
    if (!open(faiIndex, toCString(fastaFileName), toCString(faiFileName) ) ) {
        std::cerr << "ERROR: Could not open FAI index in: " << faiFileName << std::endl;
        return -1;
    }
    std::cout << "Fai Index file " << faiFileName << " was successfully opened.\n";

    int chromoNum = numSeqs(faiIndex);

    for (int seqIter = 0; seqIter < chromoNum; ++seqIter) {
        String<Dna5> chromosome;
        readSequence(chromosome, faiIndex, seqIter);
        appendValue(genomeSet, chromosome);
    }

    if (!open(fmIndex, toCString(fmFileName)))
    {
        std::cerr << "ERROR: Could not load the FM index file: " << fmFileName << "\n";
        return -1;
    }

    Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);
    goDown(bwt_iter, String<Dna5>("CGGCTTGGTGCTCACGCACACAGGAAAGTCCTT"), Rev() );
    std::cout << representative(bwt_iter) << std::endl;

    std::cout << "FM Index " << fmFileName << " was successfully loaded.\n";

    if (!open(kIndex, toCString(kFileName)))
    {
        std::cerr << "ERROR: Could not load the K-mer index to file!\n";
        return -1;
    }
    std::cout << "K-mer Index file " << kFileName << " was successfully loaded.\n";

    //loadStreamline(toCString(streamlinedFileName) );

    if (globalConfig.load_ets) {
        std::ifstream input_file;
        input_file.open(toCString(etsFileName), std::ios::out | std::ios::binary);

        resize(flagSet, length(genomeSet));

        for (unsigned seq_iter = 0; seq_iter < length(genomeSet); ++seq_iter) {
            resize(flagSet[seq_iter], length(genomeSet[seq_iter] ) );

            for (unsigned offset_iter = 0; offset_iter < length(genomeSet[seq_iter]); ++offset_iter)
                input_file.read((char*) &flagSet[seq_iter][offset_iter], sizeof(unsigned) );
        }

    }

    return 0;
}


void generateFlag() {
    unsigned max_unit_num = globalConfig.max_units;
    unsigned unit_length = globalConfig.unit_len;

    for (unsigned seq_iter = 0; seq_iter < length(genomeSet); ++seq_iter) {

        for (unsigned offset_iter = 0; offset_iter < length(genomeSet[seq_iter]); ++offset_iter) {
            processLoc(seq_iter, offset_iter);

            if (offset_iter % 100000 == 0)
                std::cout << seq_iter << ":" << offset_iter << std::endl;
        }

        /*
        for (unsigned offset_iter = 0; offset_iter <= length(genomeSet[seq_iter]) - max_unit_num * unit_length; ++offset_iter) {
            Infix<String<Dna5> >::Type seq = infix(genomeSet[seq_iter], offset_iter, offset_iter + max_unit_num * unit_length);
            flagSet[seq_iter][offset_iter] = findFirstFailErrors(seq, max_unit_num, seq_iter, offset_iter);
        }
        */
    }
}


int main(int argc, char const ** argv)
{
    //allocate global data on heap
    streamlinedDir = new std::vector<long unsigned>;
    streamlinedPosDir = new std::vector<long unsigned>;
    streamlinedPos = new std::vector<LocStruct>;
    unitPos = new std::vector<LocStruct>;
    
    CmdOptions options;

    //resize(rows(aligner), 2);

    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res == seqan::ArgumentParser::PARSE_OK) {
        if (options.build) {
            std::cout << "building index from " << options.fastaFileName << std::endl;
            buildIndices(options.fastaFileName);
        }
        else if (options.load) {
            std::cout << "loading index for " << options.fastaFileName << std::endl;
            loadIndices(options.fastaFileName, genomeSet, options.etsFileName);
        }

        if (options.process_ref) {
            resize(multi_freq_data, 100);
            resize(flagSet, length(genomeSet));

            for (unsigned seq_iter = 0; seq_iter < length(genomeSet); ++seq_iter) {
                resize(flagSet[seq_iter], length(genomeSet[seq_iter] ) );
                for (unsigned offset_iter = 0; offset_iter < length(genomeSet[seq_iter]); ++offset_iter)
                    flagSet[seq_iter][offset_iter] = 0;
            }

            generateFlag();

            CharString baseFileName = options.fastaFileName;
            while (back(baseFileName) != '.')
                eraseBack(baseFileName);

            eraseBack(baseFileName);
            CharString etsFileName = baseFileName;
            append(etsFileName, ".ets");

            std::ofstream output_file;
            output_file.open(toCString(etsFileName), std::ios::out | std::ios::binary);

            for (unsigned seq_iter = 0; seq_iter < length(genomeSet); ++seq_iter) {
                for (unsigned num_iter = 0; num_iter < length(genomeSet[seq_iter]); ++num_iter) {
                    output_file.write((char*) &(flagSet[seq_iter][num_iter]), sizeof(unsigned));
                    //std::cout << (unsigned) flagSet[seq_iter][num_iter] << " ";
                }
                //std::cout << std::endl;
            }

            output_file.close();
/*
            for (unsigned seq_iter = 0; seq_iter < length(genomeSet); ++seq_iter) {

                for (unsigned num_iter = 0; num_iter < length(genomeSet[seq_iter]); ++num_iter)
                    std::cout << (unsigned) flagSet[seq_iter][num_iter] << " ";

                std::cout << std::endl;
            }

            std::ifstream input_file;
            input_file.open(toCString(etsFileName), std::ios::in | std::ios::binary);

            for (unsigned seq_iter = 0; seq_iter < length(genomeSet); ++seq_iter) {
                for (unsigned num_iter = 0; num_iter < length(genomeSet[seq_iter]); ++num_iter) {
                    unsigned flag;
                    input_file.read((char*) &flag, sizeof(unsigned) );
                    std::cout << flag << " ";
                }
                std::cout << std::endl;
            }
                */
        }

        if (options.benchmark) {
            SeqFileIn seqFileIn;

            if (!open(seqFileIn, toCString(options.benchmarkFileName))) {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
            }

            StringSet<CharString> ids;
            StringSet<Dna5String> seqs;
            StringSet<CharString> quals;

            readRecords(ids, seqs, quals, seqFileIn);

            unsigned total_seed_num = 0;
            unsigned total_ets_num = 0;
            unsigned total_freq_num = 0;
            unsigned total_freq_num_naive = 0;

            unsigned read_num = 0;

            unsigned seed_num;
            unsigned ets_num;
            unsigned freq_num;
            unsigned freq_num_naive;

            for (unsigned seq_iter = 0; seq_iter < length(seqs); ++seq_iter) {
                bool have_N = false;

                for (int bp_iter = 0; bp_iter < length(seqs[seq_iter]); ++bp_iter) {

                    if (seqs[seq_iter][bp_iter] == 'N')
                        have_N = true;
                }

                if (have_N)
                    continue;

                processRead(seqs[seq_iter], globalConfig.max_units - 1, seed_num, ets_num, freq_num);
                processReadNaive(seqs[seq_iter], globalConfig.max_units - 1, freq_num_naive);

                //std::cout <<freq_num_naive << std::endl;
                //break;
                //
                total_seed_num += seed_num;
                total_ets_num += ets_num;
                total_freq_num += freq_num;
                total_freq_num_naive += freq_num_naive;

                ++read_num;
            }

            std::cout << "read_num: " << read_num << std::endl;
            std::cout << "total_seed_num: " << total_seed_num << std::endl;
            std::cout << "total_ets_num: " << total_ets_num << std::endl;
            std::cout << "total_freq_num: " << total_freq_num << std::endl;
            std::cout << "total_freq_num_naive: " << total_freq_num_naive << std::endl;

        }

        if (options.profile) {
            profileRepeatGenome(options.start_length, options.end_length, options.target_freq);
/*
            for (int i = 0; i < base_freq_ary.size(); ++i) {
                std::cout << base_freq_ary[i] << std::endl;
                std::cout << major_freq_ary[i] << std::endl;
            }
*/       
        }

        if (options.interactive) {
            std::string temp;
            String<Dna5> query;
            unsigned ext_length;
            while(true) {
                std::cout << "please provide the seed: ";
                std::cin >> temp;
                query = temp;
                std::cout << "please provide further extension length: ";
                std::cin >> ext_length;
                showExtensionFrequencies(query, ext_length);
            }
        }
    }
    else
        return -1;

    //while(find(fmFinder, "GTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTAC"))
        //std::cout << position(fmFinder) << std::endl;

    /*
    goBegin(fmFinder);
    clear(fmFinder);

    while(find(fmFinder, String<Dna5>("GGTGGGCAGTTCTGGAATGGTGCCAGGGGCAGAGGGGGCAATGCCGGGGCCCAGGTC")))
        std::cout << position(fmFinder) << std::endl;

    std::cout << "-----" << std::endl;

    //goBegin(fmFinder);
    //clear(fmFinder);

    Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);
    goDown(bwt_iter, String<Dna5>("GGTGGGCAGTTCTGGAATGGTGCCAGGGGCAGAGGGGG") );

    std::cout << "good so far" << std::endl;

    for (unsigned i; i < length(getOccurrences(bwt_iter) ); i++)
        std::cout << "i: " << i << getOccurrences(bwt_iter)[0] << std::endl;
    */

/*    
    std::cout << "+++++" << std::endl;

    goDown(bwt_iter, String<Dna5>("GTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTAC") );

    for (unsigned i; i < length(getOccurrences(bwt_iter) ); i++)
        std::cout << getOccurrences(bwt_iter)[0] << std::endl;
*/



    return 0;
}
    /*

    DnaString text = "AAAACACAGTTTGA";
    Shape<Dna, UngappedShape<3> > myShape;

    std::cout << hash(myShape, begin(text)) << '\t';
    for (unsigned i = 1; i < length(text) - length(myShape) + 1; ++i)
        std::cout << hashNext(myShape, begin(text) + i) << '\t';

    std::cout << std::endl;

    StringSet<String<Dna5> > genomeSet;
    appendValue(genomeSet, "TTATTAAGCGTATAGCCCTATAAATATAA");
    appendValue(genomeSet, "TATAATATAACCCAATATAACCCAAGGGG");
    typedef FMIndexConfig<void, uint32_t>     TConfig;
    Index<StringSet<String<Dna5> >, FMIndex<void, TConfig> > fmIndex(genomeSet);
    Finder<Index<StringSet<String<Dna5> >, FMIndex<void, TConfig> > > fmFinder(fmIndex);

    while (find(fmFinder, "TATAA"))
    {
        std::cout << position(fmFinder) << std::endl;
    }
    

    goBegin(fmFinder);
    clear(fmFinder);

    while (find(fmFinder, "TATAA"))
    {
        std::cout << fmFinder.range.i2 - fmFinder.range.i1 << std::endl;
        std::cout << '[' << beginPosition(fmFinder) << ',' << endPosition(fmFinder) << ")\t" << infix(fmFinder) << std::endl;
        std::cout << position(fmFinder) << std::endl;
    }

    String<Dna5> read = "AACAACCACACCCCAATTGGGAGCGTAACTGGACA";
    String<Dna5> ref = "CCATATTMGAGAGGGCCAAACCCTAAACCACACCCCAATTGGGAGCGTAACTGGACATTAGGGGACCCAGGACTACCAGTAGGATAC";

    Infix<String<Dna5> >::Type ref_inf = infix(ref, 17,59);
    std::cout << read << " " << length(read) << std::endl;
    std::cout << ref_inf << " " << length(read) << std::endl;
    std::cout << ref << std::endl;

    typedef String<Dna5> TSequence;
    typedef Align<TSequence, ArrayGaps> TAlign;

    TAlign align;

    resize(rows(align), 2);
    assignSource(row(align, 0), read);
    assignSource(row(align, 1), ref_inf);

    int score = globalAlignment(align, Score<int, Simple>(1, -1, -1), AlignConfig<true, false, true, false>());
    std::cout << "Score: " << score << std::endl;
    std::cout << align << std::endl;


    return 0;
}

*/
