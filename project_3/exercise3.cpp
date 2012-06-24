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
#include <list>
using namespace std;

struct node{
    int frequency;
    char character;
    node *left;
    node *right;
    char path;
    int id, parent_id; // nr within array
bool operator<(node& nd)
    {
    return frequency < nd.frequency;
    }
};
bool nodeComparator(const node& lhs, const node& rhs){
    return lhs.frequency < rhs.frequency;
}


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
void buildBWT(::seqan::String<unsigned> &suffixArray, string &sequence,string &bwt) {
    for (int i = 0; i < sequence.length(); ++i) {
        if (suffixArray[i] > 0)
            bwt.append(sequence, suffixArray[i] - 1, 1);
        else
            bwt.append("$");
    }
}

/**
 *  move-to-front algorithm (according to script p. 11001)
 **/
void moveToFront(string &alphabet, string &L, int *R) {
    int sigma = alphabet.length();
    int *M = new int[sigma];
    for (int i = 0; i < sigma; i++)
        M[i] = i;

    int x;
    for (int unsigned i = 0; i < L.length(); i++) {
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
void moveToFrontDecoding(string &alphabet, int *R, int sizeOfR,
		string decodedL) {

	for (int r = 0; r < sizeOfR; r++) {
		decodedL.push_back(alphabet.at(R[r]));
		string temp;
		temp.push_back(alphabet.at(R[r]));
		alphabet.insert(0, temp);
		alphabet.erase(R[r] + 1, 1);
	}

}

/**
 * Function to revert the BWT transformation and get the original text.
 */
void lfMapping() {


}

void getUsedSymbolsAndFrequencies(int *ascii_table, string &bwt, string &alphabet, vector<int> &frequencies) {
    string used_characters;
    // scan characters in bwt and increase corresponding counter in ascii-table
    for (int unsigned i = 0; i < bwt.length(); i++)
            ascii_table[(int) bwt[i]]++;
    // scan ascii-table for used characters to build the alphabet
    for (int j = 0; j < 256; j++)
        if (ascii_table[j] != 0){
            used_characters.append(1, (char) j); // insert character
            frequencies.push_back(ascii_table[j]);
            }
    alphabet= used_characters;
}



/* Huffman coding */
void getHuffmanTree(string &alphabet, vector<int> &frequencies, node *nodes){
    list<node> nodes_container; // container for the nodes
    int sigma=alphabet.length();
    int processed_nodes_counter;
    int summed_frequencies=0;

    /* Initialize node array */
    for (int i = 0; i < sigma; i++) {
        nodes[i].character = alphabet[i];
        nodes[i].frequency = frequencies[i];
        nodes[i].left=NULL;     // leafs have no children
        nodes[i].right=NULL;
        nodes[i].id=i;
        nodes_container.push_back(nodes[i]);  // after initialization push in container
        summed_frequencies+=frequencies[i];
    }
    processed_nodes_counter=sigma;
    /* Huffman algorithm: description from Wikipedia */
    /* 1. Create a leaf node for each symbol and add it to the priority queue.
       2. While there is more than one node in the queue:
        - Remove the two nodes of highest priority (lowest probability) from the queue
        - Create a new internal node with these two nodes as children and with probability equal to the sum of the two nodes' probabilities.
        - Add the new node to the queue.
        3. The remaining node is the root node and the tree is complete.
     */
    while (nodes_container.size()>1) {
        nodes_container.sort();// sort ascending according frequencies

        node left=nodes_container.front();
        nodes_container.pop_front();
        node right=nodes_container.front();
        nodes_container.pop_front();
        /* Make new node and update references */
        // Parent
        nodes[processed_nodes_counter].id=processed_nodes_counter;
        nodes[processed_nodes_counter].left=&nodes[left.id];
        nodes[processed_nodes_counter].right=&nodes[right.id];
        nodes[processed_nodes_counter].frequency=left.frequency + right.frequency;
        // children
        nodes[left.id].parent_id=processed_nodes_counter;
        nodes[left.id].path='0'; // 0 for left
        nodes[right.id].parent_id=processed_nodes_counter;
        nodes[right.id].path='1'; // 1 for right
        /* add new node to the container*/
        nodes_container.push_front(nodes[processed_nodes_counter]);
        processed_nodes_counter++;
    }
    // set parent of root
    nodes[processed_nodes_counter-1].parent_id=-1;
    nodes[processed_nodes_counter-1].path='R';// root
    // check whether all nodes were in fact processed
    if(summed_frequencies!=nodes[processed_nodes_counter-1].frequency)
        cout<<"WARNING: summed frequencies at the root equals to "<<nodes[processed_nodes_counter-1].frequency<<", expected "<<summed_frequencies<<endl;

}
void getBitCodesHuffmanTree(string &alphabet, node *nodes_array, string bit_codes[]){
    int tmp_parent_id;
    string tmp;
    char tmp_path;
    for (int unsigned i = 0; i < alphabet.length(); i++) {
        if(alphabet[i]!=nodes_array[i].character)
            cout<<"WARNING: characters of the alphabet and the ones stored in tree are not the same!"<<endl;
        tmp_path=nodes_array[i].path;
        tmp_parent_id=nodes_array[i].parent_id;
        tmp.append(1,tmp_path);
        if((unsigned)tmp_parent_id!=2*alphabet.length()-2)// last element is the root
            while (tmp_parent_id!=-1) { // -1 denotes the root of the tree
                tmp_path=nodes_array[tmp_parent_id].path;
                if(tmp_path!='R') //if node is root, do not attach path
                    tmp.append(1,tmp_path);
                tmp_parent_id=nodes_array[tmp_parent_id].parent_id;
            }
        reverse(tmp.begin(),tmp.end()); // reverse order since is bottom up
        bit_codes[i]=tmp;
        tmp=""; // reset temporary string
    }
}

void deriveHuffmanCodeFromR(string bit_codes[], int *R, int nr_of_elements, string &Huffman_code){
    for (int i = 0; i < nr_of_elements; i++) {
        Huffman_code.append(bit_codes[R[i]]);
    }
}

/* Backwards: given Huffman codes and huffmanCode derive R */
/* given char index -> bitcode, invert it to bitcode -> char index */
void createMapOfBitcodes(string bit_codes[], string &alphabet, map<string,int> &char_index_of_bitcode, int &max_bit_length){
    int unsigned max_length=0;
    for (int unsigned i = 0; i < alphabet.length(); i++){ 
         char_index_of_bitcode[bit_codes[i]]= i;
         if(bit_codes[i].length()>max_length)
             max_length=bit_codes[i].length();
    }
    max_bit_length=max_length;
}

//void getNrOfCharactersEncodes(string bit_codes[], string &alphabet, string &Huffman_code, map<string,int> bitcode_of_char_index){
void deriveRFromHuffmanCode(map<string,int> &char_index_of_bitcode,string &Huffman_code,int max_bit_length, vector<int> &R){
    int unsigned start=0;
    string tmp;
    while (start<Huffman_code.length()) {
        for (int i = 1; i <= max_bit_length;  i++) {// not optimal upper bound...
            tmp=Huffman_code.substr(start,i);
            if (char_index_of_bitcode.count(tmp)>0) {// map contains element for given key
                R.push_back(char_index_of_bitcode[tmp]); // add int to R
                start+=tmp.length();
                i=0;
            }
        }
    }
//    int unsigned min_len, max_len; // minimal and maximal length of bit codes
//    min_len=bit_codes[0].length();
//    max_len=bit_codes[0].length();
//    for (int unsigned i = 1; i < alphabet.length(); i++) {
//        if(bit_codes[i].length()<min_len)
//            min_len=bit_codes[i].length();
//        if(bit_codes[i].length()>max_len)
//            max_len=bit_codes[i].length();
//    }
}
/*  ####################### Functions to write/read files #####################################################*/
/* Writing in an output file the compressed data. */
void writeCompressedOutputfile(string &alphabet, string bit_codes[], string &HuffmanCode,string &seq_name) {
    ofstream outputfile("outputlie_c"); // Creating a file to wite to.
    string separator=" ";
    try {
        outputfile << "Alphabet" << "\n" << alphabet << "\n" << "HuffmanCodes\n";
        for (int unsigned i = 0; i < alphabet.length(); i++)
            outputfile<< i<< separator<< bit_codes[i]<<"\n";

        outputfile<< seq_name << "\n";
        // huffman code
        int unsigned width=80;
        if(HuffmanCode.length()>width){
            int nr_of_rows=HuffmanCode.length()/width;
            int rest=HuffmanCode.length()%width;
            for (int i = 0; i < nr_of_rows; i++)
                outputfile << HuffmanCode.substr(i * width, width) << "\n";
            if (rest != 0)
                outputfile << HuffmanCode.substr(nr_of_rows * width, rest);
        }else
            outputfile<<HuffmanCode;
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
            // to uncomment for the final version
//            ::seqan::String<char> text = sequence;
//            ::seqan::String<unsigned> suffixArray;
//            ::seqan::resize(suffixArray, ::seqan::length(text));
//            ::seqan::createSuffixArray(suffixArray, text, ::seqan::Skew7());
            /* Derive BWT from suffix array */
//            buildBWT(suffixArray, sequence, bwt);
            

            //####################################### test @todo: delete!
            string example = "mississippi";
            example.append("$");
            ::seqan::String<char> text = example;
            ::seqan::String<unsigned> suffixArray;

            ::seqan::resize(suffixArray, ::seqan::length(text));
            ::seqan::createSuffixArray(suffixArray, text, ::seqan::Skew7());
            /* Derive BWT from suffix array */
            buildBWT(suffixArray, example, bwt);
//          ################################### end of sector to delete
            cout << "BWT(L): " << bwt << endl;
            
            // Determine alphabet
            int ascii[256] = { 0 }; // extended ASCII table for 256 symbols
            string alphabet;
            vector<int>frequencies;
            getUsedSymbolsAndFrequencies(ascii, bwt, alphabet,frequencies);

            cout << "Alphabet: " << alphabet << endl;
            printVector(frequencies);

            /* Move_to_front */
            int *R = new int[bwt.length()];
            moveToFront(alphabet, bwt, R);
            
            
            // ########################test: example from script @ todo: DELETE
            int ascii_example[256] = { 0 };
            string L = "aooooaaiioaeieeii";
            int *R_example = new int[L.length()];
            string alphabet_example;
            vector<int> frequencies_example;

            getUsedSymbolsAndFrequencies(ascii_example, L,alphabet_example,frequencies_example);

            cout << "Move_to_front \nalphabet example: " << alphabet_example
                            << endl;
            printVector(frequencies_example);
            moveToFront(alphabet_example, L, R_example);
            cout << "R: " << endl;
            printIntArray(R_example, 17);
            // ###############################################################
            /* Huffman coding */
            //@ TODO: uncomment
            //node *nodes_array = new node [2*alphabet.length()-1];  //binary trees have 2N-1 nodes by N leafs
            //getHuffmanTree(alphabet,frequencies,nodes_array);
//            string *huff_codes = new string[alphabet.length()];
//            getBitCodesHuffmanTree(alphabet,nodes_array, huff_codes);
//            string HuffmanCode; // to store the bit code derive from R
//            deriveHuffmanCodeFromR(huff_codes,R,bwt.length(),HuffmanCode);
            
            // test
            node *nodes_example = new node [2*alphabet_example.length()-1];  //binary trees have 2N-1 nodes by N leafs
            getHuffmanTree(alphabet_example,frequencies_example,nodes_example);
            string *bit_codes = new string[alphabet_example.length()];
            getBitCodesHuffmanTree(alphabet_example,nodes_example, bit_codes);
                 for (int unsigned i = 0; i < alphabet_example.length(); i++) 
                        cout<<"path for i: "<<i<<" "<<bit_codes[i]<<endl;

            string HuffmanCode_example;
            deriveHuffmanCodeFromR(bit_codes,R_example,17,HuffmanCode_example);
            cout<<"HuffmanCode: "<<HuffmanCode_example<<endl;
            
            
            //writeUncompressedOutputfile(seq_name, sequence, 80);
            writeCompressedOutputfile(alphabet_example, bit_codes, HuffmanCode_example,seq_name);
            
            map<string,int> probe;
            vector<int> R_probe;
            int max_bit_length;
            createMapOfBitcodes(bit_codes,alphabet_example,probe,max_bit_length);
            
            cout<<"map size: "<<probe.size()<<endl;
            deriveRFromHuffmanCode(probe,HuffmanCode_example,max_bit_length,R_probe);
            cout<<"R probe size: "<<R_probe.size()<<endl;
            printVector(R_probe);
//            delete []R; @todo: to be uncomment at the end
//            delete []nodes_array;
//            delete []huff_codes;
     

    } else if (mode == 'x') {
		/* mode x */

//    [] reads a file containing the BWT compression of some sequence
//    [x] writes the uncompressed sequence into a fasta outfile. 
		//writeUncompressedOutputfile(seq_name, sequence, 80);
	}
	/* Ending */
//	cout << "---- Results ---- " << endl;
//	cout << "Name of the sequence: " << seq_name << endl;
//	cout << "Sequence: " << sequence.substr(0, 10) << endl;

	time_int(1); // print out elapsed time
	cout << endl;
	return 0;
}
