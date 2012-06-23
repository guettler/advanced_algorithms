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
#include <seqan/index.h>
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
            cout<< "File could not be opened. Please check path and file's name..."<< endl;
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
/************************* NEW FUNCTIONS FOR EXERCISE 3***********************/
/* Construct BWT from the suffix array (see 10.5 Lemma 4)*/
void buildBWT(::seqan::String<unsigned> &suffixArray, string &sequence, string &bwt) {
	for (int i = 0; i < sequence.length(); ++i) {
            if (suffixArray[i] > 0)
                bwt.append(sequence, suffixArray[i] - 1, 1);
            else
                bwt.append("$");
	}
}
/* move-to-front algorithm (according to script p. 11001) */
void moveToFront(string &alphabet, string &L, int *R) {
    int sigma = alphabet.length();
    int *M = new int[sigma];
    for (int i = 0; i < sigma; i++)
        M[i] = i;

    int x;
    for (int i = 0; i < L.length(); i++) {
        x = alphabet.find(L[i]); // since alphabet is sorted, position equals rank
        R[i] = M[x];
        for (int j = 0; j < sigma; j++)
            if (M[j] < M[x])
                M[j] = M[j] + 1;
        M[x] = 0;
    }
    delete[] M;
}

/**
 * Decode function for the move-to-front encoding.
 */
void lToFMapping(string &alphabet, string &L, int *R) {

}

string getUsedSymbols(int *ascii_table, string &bwt) {
    string used_characters;
    // scan characters in bwt and increase corresponding counter in ascii-table
    for (int i = 0; i < bwt.length(); i++)
            ascii_table[(int) bwt[i]]++;
    // scan ascii-table for used characters to build the alphabet
    for (int j = 0; j < 256; j++)
            if (ascii_table[j] != 0)
                    used_characters.append(1, (char) j); // insert character
    return used_characters;
}

/* Writing in an output file the compressed data. */
void writeCompressedOutputfile(string &alphabet, string &HuffmanCodes, string &seq_name) {
    ofstream outputfile("outputlie_c");// Creating a file to wite to.

    try {
        outputfile << "Alphabet" << "\n" << alphabet << "\n" << "HuffmanCodes\n"
                    << 0 << "@todo: code for 0" << "\n" << 1 << "code for 1" << "\n"
                    << 2 << "code for 2" << "\n" << 3 << "code for 3" << "\n"
                    << seq_name << "\n" << "@todo: Huffman code of sequence";

        outputfile.close();
    } catch (exception e) {
        cout << "Error during writing!" << endl;
        cout << e.what() << endl;
    }

}

/* Counterpart of last function: given the name of a compressed file,
 * read and save contained information (e.g. alphabet, sequence name, etc.)
 */
void readCompressedFile(string &file_name, string &alphabet_in,
		string &huffman_codes, string &huffman_seq, string &seq_name_in,
		string &sequence_in) {

}

/* Create a output file (FASTA format-> width=80) for a given sequence. */
void writeUncompressedOutputfile(string &seq_name, string &sequence, int width) {
    int seq_length = sequence.length();
    int nr_of_rows = seq_length / width;
    int rest = seq_length % width;
    ofstream outputfile("outputfile_x.fasta");

    try {
        outputfile << seq_name << "\n";
        for (int i = 0; i < nr_of_rows; i++)
                outputfile << sequence.substr(i * width, width) << "\n";
        if (rest != 0)
                outputfile << sequence.substr(nr_of_rows * width, rest);
        outputfile.close();

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
        cout << "Welcome..." << endl << "---- Introduced arguments ----"<< endl;
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
            //    [x] reads a single sequence from a fasta file
            //    [x] calculates the BWT of that sequence
            //    [] implements [x]move-to-front encoding and []Huffman coding to compress the BWT
            //    [] writes the Huffman code into an outfile (without format)

            /* Read fasta files*/
            cout << "---- Reading ----" << endl;
            sequence = readGenome(input_file, seq_name);
            sequence.append("$");
            // test
            string example = "mississippi";
            example.append("$");
            ::seqan::String<char> text = example;
            //::seqan::String<char> text = sequence;
            ::seqan::String<unsigned> suffixArray;

            ::seqan::resize(suffixArray, ::seqan::length(text));
            ::seqan::createSuffixArray(suffixArray, text, ::seqan::Skew7());
            /* Derive BWT from suffix array */
    //		buildBWT(suffixArray, sequence, bwt);
            buildBWT(suffixArray, example, bwt);

            cout << "BWT(L): " << bwt << endl;
            // Determine alphabet
            int ascii[256] = { 0 }; // extended ASCII table for 256 symbols
            string alphabet = getUsedSymbols(ascii, bwt);
            cout << "Alphabet: " << alphabet << endl;

            /* Move_to_front */
            int *R = new int[bwt.length()];
            moveToFront(alphabet, bwt, R);

            // test: example from script
            int ascii_example[256] = { 0 };
            string L = "aooooaaiioaeieeii";
            int *R_example = new int[L.length()];
            string alphabet_example = getUsedSymbols(ascii_example, L);
            cout << "Move_to_front \nalphabet example: " << alphabet_example
                            << endl;
            moveToFront(alphabet_example, L, R_example);
            cout << "R: " << endl;
            printIntArray(R_example, 17);

            //writeUncompressedOutputfile(seq_name, sequence, 80);
            writeCompressedOutputfile(alphabet, L, seq_name);
    //                delete []R; @todo: to be uncomment at the end

    } else if (mode == 'x') {
		/* mode x */

//    [] reads a file containing the BWT compression of some sequence
//    [] writes the uncompressed sequence into a fasta outfile. 
		//writeUncompressedOutputfile(seq_name, sequence, 80);
	}
	/* Ending */
	cout << "---- Results ---- " << endl;
	cout << "Name of the sequence: " << seq_name << endl;
	cout << "Sequence: " << sequence.substr(0, 10) << endl;

	time_int(1); // print out elapsed time
	cout << endl;
	return 0;
}
