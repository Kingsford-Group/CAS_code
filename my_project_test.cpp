#include <iostream>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/align.h>
#include <seqan/seq_io.h>
#include <vector>
#include <algorithm>
#include <string>

//#define __SEGFAULT__
#define __UNIT_LEN__  12
#define __TOGGLE_THRESHOLD__  3
#define __HALT_THRESHOLD__  1


using namespace seqan;

typedef FastFMIndexConfig<void, uint64_t> TConfig;
typedef Index<StringSet<String<Dna5> >, BidirectionalIndex<FMIndex<void, TConfig> > > BWTIndex;
typedef Index<StringSet<String<Dna5> >, IndexQGram<UngappedShape<__UNIT_LEN__> > > KIndex;


struct LocStruct {
    long unsigned seq_no;
    long unsigned seq_offset;
    unsigned support_units;
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



BWTIndex fmIndex;
Finder<BWTIndex> fmFinder;
KIndex kIndex;

std::vector<long unsigned> streamlinedDir;
std::vector<long unsigned> streamlinedPosDir;
std::vector<LocStruct> streamlinedPos;
std::vector<LocStruct> unitPos;

void storeStreamline (std::string filename) {
    long unsigned last_dir = streamlinedDir[length(indexDir(kIndex) ) - 1];
    long unsigned last_pos = streamlinedPosDir[last_dir];

    long unsigned dir_size = last_dir + 1;
    long unsigned pos_size = last_pos + 1;

    std::ofstream output_file;
    output_file.open(filename.c_str(), std::ios::out | std::ios::binary);

    output_file.write((char*) &(dir_size), sizeof(long unsigned) );
    output_file.write((char*) &(pos_size), sizeof(long unsigned) );

    for (long unsigned hash = 0; hash < length(indexDir(kIndex) ); ++hash)
        output_file.write((char*) &(streamlinedDir[hash]), sizeof(long unsigned) );
    
    for (long unsigned dir_idx = 0; dir_idx <= last_dir; ++dir_idx)
        output_file.write((char*) &(streamlinedPosDir[dir_idx]), sizeof(long unsigned) );

    for (long unsigned pos_idx = 0; pos_idx <= last_pos; ++pos_idx) {
        output_file.write((char*) &(streamlinedPos[pos_idx].seq_no), sizeof(long unsigned) );
        output_file.write((char*) &(streamlinedPos[pos_idx].seq_offset), sizeof(long unsigned) );
        output_file.write((char*) &(streamlinedPos[pos_idx].support_units), sizeof(unsigned) );
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

    streamlinedDir.resize(length(indexDir(kIndex) ) );
    streamlinedPosDir.resize(dir_size);
    streamlinedPos.resize(pos_size);

    for (long unsigned hash = 0; hash < length(indexDir(kIndex) ); ++hash)
        input_file.read((char*) &(streamlinedDir[hash]), sizeof(long unsigned) );
    
    for (long unsigned dir_idx = 0; dir_idx < dir_size; ++dir_idx)
        input_file.read((char*) &(streamlinedPosDir[dir_idx]), sizeof(long unsigned) );

    for (long unsigned pos_idx = 0; pos_idx < pos_size; ++pos_idx) {
        input_file.read((char*) &(streamlinedPos[pos_idx].seq_no), sizeof(long unsigned) );
        input_file.read((char*) &(streamlinedPos[pos_idx].seq_offset), sizeof(long unsigned) );
        input_file.read((char*) &(streamlinedPos[pos_idx].support_units), sizeof(unsigned) );
    }

    input_file.close();
}

unsigned streamlineUnit(Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter, unsigned left_extend, unsigned right_extend, unsigned streamlined_loc_num, unsigned unit_num, unsigned ori_left_extend) {
    auto locs = getOccurrences(bwt_iter, Fwd());
    unsigned count_freq = length(locs);
    unsigned next_loc_num;
    unsigned left_offset = ori_left_extend - left_extend;

//    std::cout << representative(bwt_iter) << std::endl;

    if (right_extend == 0 && left_extend == 0) {
        next_loc_num = streamlined_loc_num + 1;

        if (unitPos.size() < next_loc_num)
            unitPos.resize(1.5 * unitPos.size() );

        unitPos[streamlined_loc_num].seq_no = getSeqNo(locs[0]);
        unitPos[streamlined_loc_num].seq_offset = getSeqOffset(locs[0]) + left_offset;
        unitPos[streamlined_loc_num].support_units = unit_num;

        return next_loc_num;
    }

    if (count_freq <= __HALT_THRESHOLD__) {
        next_loc_num = streamlined_loc_num + count_freq;

        if (unitPos.size() < next_loc_num)
            unitPos.resize(1.5 * unitPos.size() );

        for (unsigned i = 0; i < count_freq; ++i) {
            unitPos[streamlined_loc_num + i].seq_no = getSeqNo(locs[i]);
            unitPos[streamlined_loc_num + i].seq_offset = getSeqOffset(locs[i]) + left_offset;
            unitPos[streamlined_loc_num + i].support_units = 0;
        }

        return next_loc_num;
    }

    if (right_extend != 0) {
        goDown(bwt_iter, Rev() );
        streamlined_loc_num = streamlineUnit(bwt_iter, left_extend, right_extend - 1, streamlined_loc_num, unit_num, ori_left_extend);
        while (goRight(bwt_iter, Rev() ) )
            streamlined_loc_num = streamlineUnit(bwt_iter, left_extend, right_extend - 1, streamlined_loc_num, unit_num, ori_left_extend);
    }
    else {
        goDown(bwt_iter, Fwd() );
        streamlined_loc_num = streamlineUnit(bwt_iter, left_extend - 1, right_extend, streamlined_loc_num, unit_num, ori_left_extend);
        while (goRight(bwt_iter, Fwd() ) )
            streamlined_loc_num = streamlineUnit(bwt_iter, left_extend - 1, right_extend, streamlined_loc_num, unit_num, ori_left_extend);

    }
    
    return streamlined_loc_num;
}

void streamlineKIndex(StringSet<String<Dna5> > & ref, unsigned unit_num) {
    streamlinedDir.resize(length(indexDir(kIndex) ) );
    streamlinedPosDir.resize((long unsigned) (100 + 0.001 * length(ref) ) );
    streamlinedPos.resize((long unsigned) (1000 + 0.1 * length(ref) ) );
    unitPos.resize(1000);

    long unsigned pos_dir_idx = 0;
    long unsigned pos_idx = 0;

    for (long unsigned hash = 0; hash < length(indexDir(kIndex) ) - 1; ++hash) {
        streamlinedDir[hash] = pos_dir_idx;

        if (indexDir(kIndex)[hash + 1] - indexDir(kIndex)[hash] >= __TOGGLE_THRESHOLD__) {
            String<Dna5> unit_string;
            unhash(unit_string, hash, __UNIT_LEN__);

            Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);
            goDown(bwt_iter, unit_string, Rev() );

            //std::cout << unit_string << std::endl;

            for (int unit_iter = 0; unit_iter < unit_num; ++unit_iter) {
                streamlinedPosDir[pos_dir_idx] = pos_idx;

                unsigned streamlined_loc_num = 0;
                streamlined_loc_num = streamlineUnit(bwt_iter, unit_num + unit_iter * __UNIT_LEN__, unit_num + (unit_num - 1 - unit_iter) * __UNIT_LEN__, streamlined_loc_num, unit_num, unit_num + unit_iter * __UNIT_LEN__);

                //std::cout << streamlined_loc_num << " ";

                sort(unitPos.begin(), unitPos.begin() + streamlined_loc_num, LocPosCmp() );

                unsigned long next_pos_dir_idx = pos_dir_idx + 1;
                unsigned long next_pos_idx = pos_idx + streamlined_loc_num;

                if (streamlinedPosDir.size() <= next_pos_idx + 1)
                    streamlinedPosDir.resize(1.5 * streamlinedPosDir.size() );

                if (streamlinedPos.size() <= next_pos_idx + 1)
                    streamlinedPos.resize(1.5 * streamlinedPos.size() );

                std::copy(unitPos.begin(), unitPos.begin() + streamlined_loc_num, streamlinedPos.begin() + pos_idx);

                pos_dir_idx = next_pos_dir_idx;
                pos_idx = next_pos_idx;
            }

            //std::cout << std::endl;
        }
    }

    streamlinedDir[length(indexDir(kIndex) ) - 1] = pos_dir_idx;
    streamlinedPosDir[pos_dir_idx] = pos_idx;
}


int main()
{
    StringSet<String<Dna5> > ref;
    appendValue (ref, "TCGTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTAC");
    appendValue (ref, "ATGTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTAC");
    appendValue (ref, "CAGTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTACAAAAAAAAAA");
    appendValue (ref, "AGGTGCTCCTCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTACTTTTTTTTTT");

/*  // UNCOMMENT ME!!!!
    FaiIndex faiIndex;
    build(faiIndex, "chr1.fasta");

    int chromoNum = numSeqs(faiIndex);

    for (int seqIter = 0; seqIter < chromoNum; ++seqIter) {
        String<Dna5> chromosome;
        readSequence(chromosome, faiIndex, seqIter);
        appendValue(ref, chromosome);
    }
*/  // UNCOMMENT ME!!!!
    
    Infix<Dna5String>::Type test = infix(ref[0], 3, 19);

    if (!open(fmIndex, "test.bi_dir") )
	std::cout << "can't load" << std::endl;

#ifndef __SEGFAULT__
    BWTIndex temp(ref);
    fmIndex = temp;
#else
    //fmIndex = BWTIndex(ref);
#endif

    /*
    setHaystack(fmFinder, fmIndex);


    std::cout << test << std::endl;

    goBegin(fmFinder);
    clear(fmFinder);

    while(find(fmFinder, test)) {
         //std::cout << beginPosition(fmFinder) <<  " " << endPosition(fmFinder) << std::endl;
         std::cout << position(fmFinder) << std::endl;
	 long unsigned seq_no = getSeqNo(position(fmFinder));
	 long unsigned seq_offset = getSeqOffset(position(fmFinder));
	 std::cout << infix(ref[seq_no], seq_offset - 2, 23) << std::endl;
    }

    goBegin(fmFinder);
    clear(fmFinder);
    */

    String<Dna5> testRev = test;
    reverse(testRev);
    std::cout << "rev: " << testRev << std::endl;

    Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter(fmIndex);

    goDown(bwt_iter, test, Rev() );

    Iter<BWTIndex, VSTree<TopDown<ParentLinks<> > > > bwt_iter_reserve(bwt_iter);
    
    std::cout << representative(bwt_iter, Fwd()) << std::endl;
    
    goDown(bwt_iter, Rev() );
    goDown(bwt_iter, Rev() );
    //goDown(bwt_iter, Dna5String("TG"), Rev() );
    //goDown(bwt_iter, Dna5String("GA"), Fwd() );
    
    std::cout << representative(bwt_iter, Fwd()) << std::endl;

    while(goRight(bwt_iter, Rev() ) ) {
        std::cout << "going right" << std::endl;
    }
    std::cout << "!!stop going right" << std::endl;

    std::cout << representative(bwt_iter, Fwd()) << std::endl;

    for (unsigned i = 0; i < length(getOccurrences(bwt_iter, Rev())); ++i)
        std::cout << getOccurrences(bwt_iter, Fwd())[i] << std::endl;
    
    //goRight(bwt_iter);
    //goUp(bwt_iter, Rev() );
    //representative(bwt_iter, Rev());

    for (unsigned i = 0; i < length(getOccurrences(bwt_iter, Rev())); ++i)
        std::cout << getOccurrences(bwt_iter, Fwd())[i] << std::endl;

    std::cout << "reserve" << std::endl;
    for (unsigned i = 0; i < length(getOccurrences(bwt_iter_reserve, Rev())); ++i)
        std::cout << getOccurrences(bwt_iter_reserve, Rev())[i] << std::endl;

    std::cout << representative(bwt_iter, Fwd()) << std::endl;
    std::cout << representative(bwt_iter, Rev()) << std::endl;

    std::cout << "count bidirectional iterator: "  << countOccurrences(bwt_iter) << std::endl;

/*
    do {
    	std::cout << representative(bwt_iter, Rev()) << std::endl;
	if (!goDown(bwt_iter, Rev()) && !goRight(bwt_iter, Rev()))
		while (goUp(bwt_iter, Rev()) && !goRight(bwt_iter, Rev()))
			;
    }
    while (true);

    //if (!save(fmIndex, "test.bi_dir") )
	//std::cout << "can't save" << std::endl;


    std::cout << length(getOccurrences(bwt_iter)) << std::endl;

    for (int i = 0; i < length(getOccurrences(bwt_iter)); ++i) {
        long unsigned seq_no = getSeqNo(getOccurrences(bwt_iter)[i]);
        long unsigned seq_off = getSeqOffset(getOccurrences(bwt_iter)[i]);
        std::cout << getOccurrences(bwt_iter)[i] << std::endl;
        std::cout << infix(ref[seq_no], seq_off - 2, seq_off + length(test) + 4 ) << std::endl;
    }
    */

    kIndex = KIndex(ref);
    Infix<String<Dna5> >::Type seq = infix(ref[0], __UNIT_LEN__, 2 * __UNIT_LEN__);

    std::cout << length(ref) << std::endl;

    hash(indexShape(kIndex), begin(seq) );
    //hash(indexShape(kIndex), seq);
    Infix<const String<Pair<long unsigned, long unsigned, Pack> > >::Type locs = getOccurrences(kIndex, indexShape(kIndex) );

    //std::cout << locs << std::endl;
    for (unsigned i = 0; i < length(locs); ++i)
        std::cout << getSeqNo(locs[i]) << " " << getSeqOffset(locs[i]) << std::endl;

    std::cout << "counts" << std::endl;
    std::cout << indexCounts(kIndex) << std::endl;
    
    std::cout << "shape" << std::endl;

    auto hash_val = value(indexShape(kIndex));
    std::cout << hash_val << std::endl;
    std::cout << getBucket(indexBucketMap(kIndex), hash_val) << std::endl;

    std::cout << end(indexDir(kIndex), Standard() ) - begin(indexDir(kIndex), Standard() ) << std::endl;
    std::cout << length(indexDir(kIndex) ) << std::endl;

    String<Dna5> temp_str;
    unhash(temp_str, value(indexShape(kIndex) ), __UNIT_LEN__);
    std::cout << temp_str << std::endl;

    unhash(temp_str, 0, __UNIT_LEN__);
    std::cout << temp_str << std::endl;

    unhash(temp_str, 15, __UNIT_LEN__);
    std::cout << temp_str << std::endl;

    unhash(temp_str, length(indexDir(kIndex) ) - 1, __UNIT_LEN__);
    std::cout << temp_str << std::endl;


    unsigned unit_num = 3;

//    streamlineKIndex(ref, unit_num);

//    storeStreamline("wtf");
    loadStreamline("wtf");

    for (long unsigned hash = 0; hash < length(indexDir(kIndex) ) - 1; ++hash) {
        if (indexDir(kIndex)[hash + 1] - indexDir(kIndex)[hash] >= __TOGGLE_THRESHOLD__) {
            String<Dna5> unit_string;
            unhash(unit_string, hash, __UNIT_LEN__);

            std::cout << unit_string << std::endl;
            std::cout << indexDir(kIndex)[hash + 1] - indexDir(kIndex)[hash] << std::endl;

            long unsigned dir_idx = streamlinedDir[hash];

            for (int i = 0; i < unit_num; ++i)
                std::cout << streamlinedPosDir[dir_idx + i + 1] - streamlinedPosDir[dir_idx + i] << " ";

            std::cout << std::endl;

            long unsigned pos_idx = streamlinedPosDir[dir_idx];

            for (int i = 0; i < streamlinedPosDir[dir_idx + 1] - streamlinedPosDir[dir_idx]; i++)
                std::cout << "<" << streamlinedPos[pos_idx + i].seq_no << ", " << streamlinedPos[pos_idx + i].seq_offset << ">:" << streamlinedPos[pos_idx + i].support_units << " ";

            std::cout << std::endl;
        }
    }
/*
    goBegin(fmFinder);
    clear(fmFinder);
    find(fmFinder, DnaString('A'));

    TIterator initial_fm_pos = fmFinder.range.i1;

    std::cout << fmFinder.range.i2 - fmFinder.range.i1 << std::endl;
    std::cout << fmFinder.range.i1 - initial_fm_pos << std::endl;
    std::cout << fmFinder.range.i2 - initial_fm_pos << std::endl;
    goBegin(fmFinder);
    clear(fmFinder);
    find(fmFinder, DnaString('C'));
    std::cout << fmFinder.range.i2 - fmFinder.range.i1 << std::endl;
    std::cout << fmFinder.range.i1 - initial_fm_pos << std::endl;
    std::cout << fmFinder.range.i2 - initial_fm_pos << std::endl;
    goBegin(fmFinder);
    clear(fmFinder);
    find(fmFinder, DnaString('G'));
    std::cout << fmFinder.range.i2 - fmFinder.range.i1 << std::endl;
    std::cout << fmFinder.range.i1 - initial_fm_pos << std::endl;
    std::cout << fmFinder.range.i2 - initial_fm_pos << std::endl;
    goBegin(fmFinder);
    clear(fmFinder);
    find(fmFinder, DnaString('T'));
    std::cout << fmFinder.range.i2 - fmFinder.range.i1 << std::endl;
    std::cout << fmFinder.range.i1 - initial_fm_pos << std::endl;
    std::cout << fmFinder.range.i2 - initial_fm_pos << std::endl;
    goBegin(fmFinder);
    clear(fmFinder);
    find(fmFinder, DnaString('N'));
    std::cout << fmFinder.range.i2 - fmFinder.range.i1 << std::endl;
    std::cout << fmFinder.range.i1 - initial_fm_pos << std::endl;
    std::cout << fmFinder.range.i2 - initial_fm_pos << std::endl;


    typedef String<Dna5> TSequence;
    typedef Align<TSequence, ArrayGaps> TAlign;

    TAlign align;

    TSequence str1("GTGGCCGGCCCTCGCTCCAGCAGCTGGACC");
    TSequence str2("ACGTGGCCGGCCCTCGCTCCAAGCAGCTGGACCCC");

    resize(rows(align), 2);
    assignSource(row(align, 0), str1);
    assignSource(row(align, 1), str2);

    int score = globalAlignmentScore(str1, str2, Score<int, Simple>(0, -1, -1), AlignConfig<false, true, true, false>(), -3, 2);//, MyersBitVector() );// Hirschberg());
    std::cout << "Score: " << score << std::endl;
    //std::cout << align << std::endl;

	*/

    return 0;
}

