#include "optimalSolverLN.h"
#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cstring>
#include <iostream>
#include <map>

using namespace std;

int main(int argc, const char* argv[]) {
	map <unsigned int, unsigned long long> frequencyCounter;

	OptimalSolverLN solver;

	if (argc == 3) {
		cout << "Loading tree " << argv[1] << " into memory.\n";
		solver.loadTree(argv[1]);
		cout << "Finished loading. \n";

		int seedNum;
		unsigned int frequency;
		//cout << "Please input number of seeds: ";
		cin >> seedNum;

		ifstream readFile;
		readFile.open(argv[2]);
		string read;

		char benchFileName[80];
		strcpy(benchFileName, argv[2]);
		char* benchName = strtok(benchFileName, ".");

        cout << "benchName: " << benchName << endl;

		ofstream output;
		output.open( (string(benchName) + string(".optimalLN") ).c_str() );
#ifdef FREQ
		ofstream freqOutput;
		freqOutput.open( (string(benchName) + string(".optimalLN_freq") ).c_str() );
#endif

		while (seedNum > 0) {
			cout << "seedNum: " << seedNum << endl;
			output << "seedNum: " << seedNum << endl;
#ifdef FREQ
			freqOutput << "seedNum: " << seedNum << endl;
#endif

			frequencyCounter.clear();

			readFile.clear();
			readFile.seekg(0, ios::beg);

			//Name
			getline(readFile, read);
			//Read
			getline(readFile, read);
			solver.init(read.length(), seedNum);
			//cout << "DNA: " << read << endl;

			//Solve
			frequency = solver.solveDNA(read);
			if (frequencyCounter.find(frequency) == frequencyCounter.end() )
				frequencyCounter[frequency] = 1;
			else
				frequencyCounter[frequency]++;

			//Trash
			getline(readFile, read);
			getline(readFile, read);

			//Name
			getline(readFile, read);

			while (!readFile.eof()) {
				//Read
				getline(readFile, read);
				//			cout << "DNA: " << read << endl;

				//Solve
				//if (read.find('N') == string::npos) {
					frequency = solver.solveDNA(read);
#ifdef FREQ
					solver.backtrack();
					solver.sortOfFreq();
					solver.printFreqs(freqOutput);
					solver.printLength(freqOutput);
#endif
					if (frequencyCounter.find(frequency) == frequencyCounter.end() )
						frequencyCounter[frequency] = 1;
					else
						frequencyCounter[frequency]++;
				//}

				//Trash
				getline(readFile, read);
				getline(readFile, read);

				//Name
				getline(readFile, read);
			}

            uint64_t total_freq = 0;
            uint64_t total_reads = 0;
            for (map<unsigned, unsigned long long>::iterator iter = frequencyCounter.begin(); iter != frequencyCounter.end(); iter++) {
                total_reads += iter->second;
                total_freq += (uint64_t) iter->first * iter->second;
                //output << iter->first << ": " << iter->second << endl;
            }

            cout << total_reads << " reads with average seed frequency of read: " << (double) total_freq / total_reads << endl;

			cin >> seedNum;

            //solver.backtrack();
			//solver.printSeeds(cout);
            
			//solver.printStats(cout);
			//solver.printStats(output);
		}
		output.close();
#ifdef FREQ
		freqOutput.close();
#endif

	}
	else if (argc == 1) {
		cout << "input seedNum and readLength: " << endl;
		int seedNum;
		int readLength;
		cin >> seedNum;
		string dummy = "";
		while (seedNum > 0) {
			cin >> readLength;
			solver.setMinLength(10);
			solver.init(readLength, seedNum);
			solver.feedL0();
			solver.solveDNA(dummy);
			solver.backtrack();
			solver.sortOfFreq();
			solver.printFreqs(cout);
			solver.printStats(cout);
			cin >> seedNum;
		}
	}
	else
	{
		cerr << "USAGE: ./testOptimalSolver <Tree File> <Read File>" << endl;
		exit(1);
	}

	return 0;
}
