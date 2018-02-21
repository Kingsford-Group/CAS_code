#include "etsSolver.h"
#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cstring>
#include <iostream>
#include <map>

using namespace std;

int main(int argc, const char* argv[]) {
	map <int, unsigned long long> frequencyCounter;
	map <unsigned, unsigned long long> seedNumCounter;

	EtsSolver solver;

	if (argc != 5) {
        cerr << "Usage: ./bin tree_path config_path ets_folder read_file" << endl;
        exit(1);
    }

    cout << "Loading tree " << argv[1] << " into memory.\n";
    solver.loadTree(argv[1], argv[2], argv[3]);
    cout << "Finished loading. \n";

    int errorTolerance;
    int frequency;
    unsigned seedNum;
    //cout << "Please input number of seeds: ";
    cin >> errorTolerance;

    std::ifstream readFile;
    readFile.open(argv[4]);
    string read;

    char benchFileName[80];
    strcpy(benchFileName, argv[4]);
    char* benchName = strtok(benchFileName, ".");

    cout << "benchName: " << benchName << endl;

    std::ofstream output;
    output.open( (string(benchName) + string(".ets_result") ).c_str() );

    while (errorTolerance > 0) {
        cout << "errorTolerance: " << errorTolerance << endl;
        output << "errorTolerance: " << errorTolerance << endl;

        frequencyCounter.clear();
        seedNumCounter.clear();

        readFile.clear();
        readFile.seekg(0, ios::beg);

        //Name
        getline(readFile, read);
        //Read
        getline(readFile, read);
        //cout << "DNA: " << read << endl;

        //Solve
        solver.solveDNA(read, errorTolerance);
        frequency = solver.getTotalFreq();
        //solver.printSeeds(cout);

        if (frequencyCounter.find(frequency) == frequencyCounter.end() )
            frequencyCounter[frequency] = 1;
        else
            frequencyCounter[frequency]++;

        if (frequency >= 0) {
            seedNum = solver.getNumOfSeeds();
            if (seedNumCounter.find(seedNum) == seedNumCounter.end() )
                seedNumCounter[seedNum] = 1;
            else
                seedNumCounter[seedNum]++;
        }

        //Trash
        getline(readFile, read);
        getline(readFile, read);

        //Name
        getline(readFile, read);

        while (!readFile.eof()) {
            //Read
            getline(readFile, read);

            //Solve
            //if (read.find('N') == string::npos) {
                solver.solveDNA(read, errorTolerance);
                frequency = solver.getTotalFreq();
                //solver.printSeeds(cout);

                if (frequencyCounter.find(frequency) == frequencyCounter.end() )
                    frequencyCounter[frequency] = 1;
                else
                    frequencyCounter[frequency]++;

                if (frequency >= 0) {
                    seedNum = solver.getNumOfSeeds();
                    if (seedNumCounter.find(seedNum) == seedNumCounter.end() )
                        seedNumCounter[seedNum] = 1;
                    else
                        seedNumCounter[seedNum]++;
                }

            //}

            //Trash
            getline(readFile, read);
            getline(readFile, read);

            //Name
            getline(readFile, read);
        }

        uint64_t total_freq = 0;
        uint64_t total_seeds = 0;
        uint64_t total_reads = 0;

        for (map<int, unsigned long long>::iterator iter = frequencyCounter.begin(); iter != frequencyCounter.end(); iter++) {
            if (iter->first >= 0) {
                total_reads += iter->second;
                total_freq += (uint64_t) iter->first * iter->second;
            }
            //output << iter->first << ": " << iter->second << endl;
        }

        for (map<unsigned, unsigned long long>::iterator iter = seedNumCounter.begin(); iter != seedNumCounter.end(); iter++) {
            total_seeds += (uint64_t) iter->first * iter->second;
            //output << iter->first << ": " << iter->second << endl;
        }

        cout << "Failed count: ";
        if (frequencyCounter.find(-1) == frequencyCounter.end() )
            cout << 0;
        else
            cout << frequencyCounter[-1];

        cout << endl << total_reads<< " reads with average seed frequency of read: " << (double) total_freq / total_reads << " average seed num of read: " << (double) total_seeds / total_reads << endl;

        //solver.printSeeds(cout);
        
        cin >> errorTolerance;
    }

    output.close();

	return 0;
}
