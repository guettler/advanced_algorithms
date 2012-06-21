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
#include <seqan/index.h> // include path muss be 'relative'
using namespace std;
/* Function prototypes (instead of header file) */
//string readGenome(string &path);
//void writeOutput(vector<vector<int> > &pos_score2);
/* Function definitions */
string readGenome(string &path, string &seq_name) {

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
			else
				seq_name = line;
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

void buildBWT(string &suffixArray, string &sequence, string &bwt) {

	int sequenceLength = sequence.length();
	cout << "Sequence length: " << sequenceLength;
}

/* Auxiliary functions to show content @TODO: remove if not used*/
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
	string input_file, output_file, sequence, seq_name, bwt;
	char mode;
	//./exercise3 <input file> <output file> <mode>
	/* Print the arguments */
	if (argc > 1) {
		cout << "Welcome..." << endl << "---- Introduced arguments ----"
				<< endl;
		for (int i = 1; i < argc; i++) {
			switch (i) {
			case 1:
				cout << "Input file:" << "\t" << argv[1] << endl;
				break;
			case 2:
				cout << "Output file:" << "\t" << argv[2] << endl;
				break;
			case 3:
				cout << "Mode:" << "\t\t" << argv[3] << endl;
				break;
			default:
				cout << "undefined input" << "\t" << argv[i] << endl;
				break;
			}
		}

		if (argc < 4)
			cout << "WARNING: not enough arguments given!" << endl;
		/* set given arguments */
		input_file = argv[1];
		output_file = argv[2];
		mode = argv[3][0]; // specify mode: c (compress) or x (extract). low case sensitive!
		//mode = 'C';
		if (mode != 'c' && mode != 'x')
			cout << "WARNING: given mode is not defined!" << endl;

	}
	/* mode c */
	if (mode == 'c') {
		/* Read fasta files*/
		cout << "---- Reading ----" << endl;
		sequence = readGenome(input_file, seq_name);
		sequence.append("$");

		//    [x] reads a single sequence from a fasta file
		//    [] calculates the BWT of that sequence
		//    [] implements move-to-front encoding and Huffman coding to compress the BWT
		//    [] writes the Huffman code into an outfile
//        string testo = "hello world!";
//        ::seqan::String<char> text = testo;
		::seqan::String<char> text = sequence;
		::seqan::String<char> pattern = "l";
		::seqan::String<unsigned> suffixArray;

		::seqan::resize(suffixArray, ::seqan::length(text));
		::seqan::createSuffixArray(suffixArray, text, ::seqan::Skew7());
		cout << suffixArray[2] << endl;

		buildBWT(suffixArray, sequence, bwt);

	} else if (mode == 'x') {
		/* mode x */

//    [] reads a file containing the BWT compression of some sequence
//    [] writes the uncompressed sequence into a fasta outfile. 
	}
	/* Ending */
	cout << "---- Results ---- " << endl;
	cout << "Name of the sequence: " << seq_name << endl;
	cout << "Sequence: " << sequence.substr(0, 10) << endl;

	time_int(1); // print out elapsed time
	cout << endl;
	return 0;
}
