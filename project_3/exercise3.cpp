/*
 * File:    exercise3.cpp
 * Authors: Group 5: Guettler, Liebers, Mattes
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <bitset>
#include <math.h>

using namespace std;
/* Function prototypes (instead of header file) */
//string readGenome(string &path);
//void writeOutput(vector<vector<int> > &pos_score2);
/* Function definitions */
string readGenome(string &path) {

	ifstream file;
	string line, input;

	try {
		file.open(path.c_str(), ifstream::in);
		if (!file)
			cout
					<< "File could not be opened. Please check path and file's name..."
					<< endl;
		while (!file.eof()) {
			getline(file, line);
			if (line.at(0) != '>') // if not id line
				input.append(line);
		}
	} catch (exception e) {
		cout << e.what() << endl;
	}

	file.close();
	return input;
}

/* Running time calculation (C style)*/
void time_int(int print) {

	static struct timeval t1; /* var for previous time stamp */
	static struct timeval t2; /* var of current time stamp */
	struct timezone tzp;

	if (gettimeofday(&t2, &tzp) == -1)
		return;

	if (print == 1) {
		double elapsed_seconds = (double) (t2.tv_sec - t1.tv_sec)
				+ ((double) (t2.tv_usec - t1.tv_usec)) / 1000000;
		printf("Time spent [%.2fs] \n", elapsed_seconds);
	}

	t1 = t2;
}

void printVector(vector<int> &v) {

	int length = v.size();

	for (int i = 0; i < length; ++i)
		cout << v[i] << " ";

	cout << endl;
}

void printTable(vector<vector<int> > &v) {

	int length = v.size();
	int width = v[0].size();

	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j)
			cout << v[i][j] << " ";
		cout << endl;
	}
}

void printIntArray(int v[], int size) {

	for (int i = 0; i < size; ++i)
		cout << "[" << i << "] " << v[i] << " " << endl;
}

void printIntTable(int **v, int rows, int columns) {

	for (int i = 0; i < rows; ++i) {
		cout << "[" << i << "] ";
		for (int j = 0; j < columns; j++)
			cout << v[i][j] << " ";
		cout << endl;
	}
}

/* Initializes an integer array with a given value. Implemented since the array
 * initializations with constructor are not reliable */
void initializeIntArray(int *arr, int &nr_of_elements, int value) {

	for (int i = 0; i < nr_of_elements; i++)
		arr[i] = value;
}

/* Function to write the results according to the given format */
/* CHANGE ME!!! */
void writeOutput(vector<vector<int> > &pos_score2) {

	ofstream outfile("hits.result"); // creates a file to write in
	string sep = ","; // symbol for separation
	string nl = "\n"; // new line
	int m = pos_score2.size(); // number of rows

	try {
		for (int i = 0; i < m; i++)
			outfile << "read_" << pos_score2[i][0] << sep << pos_score2[i][1]
					<< sep << pos_score2[i][2] << sep << pos_score2[i][3] << nl;
		outfile.close();
	} catch (exception e) {
		cout << "Error during writing!" << endl;
		cout << e.what() << endl;
	}

}

/* ########################## MAIN ########################## */
int main(int argc, char**argv) {

	time_int(0); // start timing
	string genome_file, outputFile;
	int k, q, b, threshold, w;
	/* Print the arguments */
	if (argc > 1) {
		cout << "Welcome..." << endl << "---- Introduced arguments ----"
				<< endl;
		for (int i = 1; i < argc; i++) {
			switch (i) {
			case 1:
				cout << "FASTA file:" << "\t" << argv[1] << endl;
				break;
			case 2:
				cout << "Output file:" << "\t" << argv[2] << endl;
				break;
			case 3:
				cout << "Mode:" << "\t" << argv[3] << endl;
				break;
			default:
				cout << "undefined input" << "\t" << argv[i] << endl;
				break;
			}
		}

		if (argc < 4)
			cout << "WARNING: not enough arguments given!" << endl;
		/* set given arguments */
		genome_file = argv[1];
		outputFile = argv[2];
		k = atoi(argv[3]); // number of errors
		q = atoi(argv[4]); // length of the q-grams
		b = atoi(argv[5]); // block size for the genome
	}

	/* Read fasta files*/
	cout << "---- Reading ----" << endl;
	string genome = readGenome(genome_file);

	cout << "Nr. of reads: " << m << endl;
	cout << "Text's size: " << genome.size() << endl;
	cout << "Pattern's size: " << reads[0].size() << endl;

	/* check for 1st and 2nd abort condition */
	w = reads[0].size();
	threshold = w + 1 - (k + 1) * q; //lemma 2 (q-gram)
	cout << "---- Values for the computation ----" << endl;
	cout << "Threshold: " << threshold << endl;
	if (b < w) //block size is smaller than the read length
		writeAndTerminate(1); // write a warning and terminate (jump to the end)
	else if (threshold <= 0)
		writeAndTerminate(2);
	else {
		/* Quasar */
		int nr_of_rows = (int) pow(4, q); //expected nr of permutation of a sequence of length q

		// allocate int array and vector for the scan
		int *nr_of_occ = new int[nr_of_rows](); // initialized with zero
		initializeIntArray(nr_of_occ, nr_of_rows, 0); // ensures that the array was initialized with 0
		vector<vector<int> > occurrences(nr_of_rows); // to save the start positions
		cout << "(Scanning genome w.r.t. q-gram permutations ...)" << endl;
		getNrOfOccAndPositions(genome, nr_of_occ, occurrences, q); // scan the given sequence and save nr of occurrences as well their positions

		/* Counting q-grams, division in blocks, etc */
		int nr_of_blocks, min_block_length, min_overlap, seq_length,
				read_length;
		seq_length = genome.size();
		read_length = reads[0].size(); // = w

		// set values for the minimal block length and for the minimal overlap length
		if (b < (read_length + k))
			min_block_length = read_length + k;
		else
			min_block_length = b;
		min_overlap = read_length + k - 1;

		// divide the genome into blocks
		nr_of_blocks = seq_length / min_block_length; // same as floor(). rounded down due to the minimal block length -> last block will contain remaining elements

		// block counting
		cout << "Minimal block length: " << min_block_length << endl;
		cout << "Nr. of blocks: " << nr_of_blocks << endl;
		int *blocks_counter = new int[nr_of_blocks];
		initializeIntArray(blocks_counter, nr_of_blocks, 0);
		int **blocks_range;
		// initialize table
		blocks_range = new int*[nr_of_blocks]; // # rows (m)
		for (int i = 0; i < nr_of_blocks; i++)
			blocks_range[i] = new int[2]; // # columns (n)

		// calculate the limits of each block
		calculateBlocksRange(blocks_range, seq_length, nr_of_blocks,
				min_block_length, min_overlap);

		/* main process for the computation of m reads */
		//Search for matching q-grams and count matches for each block
		vector<vector<int> > pos_score; // score vector for format: (nr. of read, start position, end position, score)
		cout << "---- Starting computation ----" << endl;
		int block_start, block_end;
		int semiglobal_counter = 0;
		for (int readNr = 0; readNr < m; readNr++) { // run ukkonen for readNr many reads, save all scores in pos_score
			calculateOccPerBlock(reads[readNr], occurrences, q, nr_of_blocks,
					min_block_length, blocks_range, blocks_counter); // scan read and increase block countig if it appears in the genome. See function definition for more details
			for (int i = 0; i < nr_of_blocks; i++) {
				if (blocks_counter[i] >= threshold) { // filtering: only the reads with block counter greater/equal threshold will be verified
					//cout<<"ReadNr: "<<readNr<<", block nr: "<<i<<" value: "<<blocks_counter[i]<<endl;//@TODO: delete/comment for the final version
					block_start = blocks_range[i][0];
					block_end = blocks_range[i][1];
					vector<vector<int> > tmp_pos_score; // for the semi-global aligner procedure
					fastUkkonen(tmp_pos_score, k, genome, block_start,
							block_end, reads[readNr], useUkkonenTrick);
					filterHitsAndBacktrack(pos_score, tmp_pos_score, k, genome,
							reads[readNr], readNr, filterResults);
					semiglobal_counter++;
				}
			}
			initializeIntArray(blocks_counter, nr_of_blocks, 0); // reset block counter for the next read
		}
		//allOccurrencesFound(genome,nr_of_occ,q,nr_of_rows);// check whether the expected nr. of occurrences (of q-grams!) was registered

		/*  free allocated space */
		// table (2 dimensional array)
		for (int i = 0; i < nr_of_blocks; i++)
			delete[] blocks_range[i];
		delete[] blocks_range;
		// arrays
		delete[] blocks_counter;
		delete[] nr_of_occ; // may be not necessary since program terminates and the allocated space will be erased anyway

		/* Ending */
		cout << "---- Results ---- " << endl;
		cout << "Nr. of occurrences: " << pos_score.size() << endl;
		cout << "Nr. of semi-global aligner iterations:  " << semiglobal_counter
				<< endl;
		writeOutput(pos_score); // export of the results
	} // end of quasar and semi-global aligning
	time_int(1); // print out elapsed time
	cout << endl;
	return 0;
}
