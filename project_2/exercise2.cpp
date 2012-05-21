/*
 * File:    exercise2.cpp
 * Authors: Group 5
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
#define match 0
#define mismatch  1
#define gap 1
#define Q 12    // maximal expected value for q

/* Function prototypes (instead of header file)*/
string readGenome(string &path);
int getNrOfReads(string &path);
void readReads(string &path, vector<string> &reads);
void time_int(int print);
void printVector(vector<int>  &v);
void printTable(vector<vector<int> > &v);
void printIntArray(int [], int );
void printIntTable(int**, int, int);
void initializeIntArray(int *, int &, int);
void fastUkkonen(vector<vector<int> > &, int &k, string &,int &, int &, string &, bool &);
void filterHitsAndBacktrack(vector<vector<int> > &, vector<vector<int> > &,int &, string &, string &, int &, bool &);
void writeOutput(vector<vector<int> > &pos_score2);
void writeAndTerminate(int);
void convBase(char &, int &);
string getQGramFromIndex(int,  int &);
void allOccurrencesFound(string &, int *, int &, int &);
void getNrOfOccAndPositions(string &, int *, vector<vector<int> > &, int &);
void calculateBlocksRange(int**, int &, int &, int &,int &);
void calculateOccPerBlock(string &, vector<vector<int> > &, int &, int &,int &, int **, int *);

/* Function definitions */
string readGenome(string &path){
    ifstream file;
    string line, input;
    try{
        file.open(path.c_str(), ifstream::in);
        if(!file)
            cout<< "File could not be opened. Please check path and file's name..."<< endl;
        while(!file.eof())
        {
            getline(file, line);
            if(line.at(0)!= '>') // if not id line
                input.append(line);
        }
    }
    catch(exception e){
        cout << e.what() << endl;
    }
    file.close();
    return input;
}

int getNrOfReads(string &path){
    ifstream file;
    string line;
    int nr_of_reads=0;
    try{
        file.open(path.c_str(), ifstream::in);
        if(!file)
            cout<< "File could not be opened. Please check path and file's name..."<< endl;
        while(!file.eof())
        {
            getline(file, line);
            if(line.at(0)== '>')
                nr_of_reads++;
        }
        file.close();
    }catch(exception e)
    {
    //		cout << e.what() << endl;
    }
    return nr_of_reads;
}

void readReads(string &path, vector<string> &reads){
    ifstream file;
    string line;
    int read_nr=-1;
    try{
        file.open(path.c_str(), ifstream::in);
        if(!file)
            cout<< "File could not be opened. Please check path and file's name..."<< endl;
        while(!file.eof()){
            getline(file, line);
            if(line.at(0)!='>')
                reads[read_nr].append(line);
            else
                read_nr++;
        }
        file.close();
    }
    catch(exception e){
        //cout << e.what() << endl;
    }
}
/* Running time calculation (C style)*/
void time_int(int print){
    static struct timeval t1; /* var for previous time stamp */
    static struct timeval t2; /* var of current time stamp */
    struct timezone tzp;

    if(gettimeofday(&t2, &tzp) == -1) return;

    if(print == 1){
    double elapsed_seconds=(double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec))/1000000;
    printf("Time spent [%.2fs] \n", elapsed_seconds);
    }
    t1 = t2;
}

void printVector(vector<int>  &v){
    int length = v.size();
    for(int i=0; i<length; ++i)
        cout << v[i]<< " ";
    cout << endl;
}
void printTable(vector<vector<int> > &v){
    int length = v.size();
    int width = v[0].size();

    for(int i=0; i<length; ++i){
        for(int j=0; j<width; ++j)
                cout << v[i][j]<< " ";
        cout << endl;
    }
}
void printIntArray(int v[], int size){
    for(int i=0; i< size; ++i)
        cout << "["<<i<<"] "<<v[i]<< " "<<endl;
}
void printIntTable(int **v, int rows, int columns){
    for(int i=0; i< rows; ++i){
        cout << "["<<i<<"] ";
        for (int j = 0; j < columns; j++) 
            cout<< v[i][j]<< " ";
        cout<<endl;
    }
}
/* Initializes an integer array with a given value. Implemented since the array
 * initializations with constructor are not reliable */
void initializeIntArray(int *arr, int &nr_of_elements, int value){
    for (int i = 0; i < nr_of_elements; i++)
        arr[i] = value;
}

/* Filter the results of fastUkkonen and do backtracking for remaining hits */
void filterHitsAndBacktrack(vector<vector<int> > &pos_score, vector<vector<int> > &tmp_pos_score,
		int &k, string &sequence, string &read, int &readNr, bool &filterResults){
    
    int m = read.size();
    int n = m + k;			// we only calculate an m+k wide part of the dp-matrix for backtracking (worst case: k gaps in the read)

    // filter and backtracking
    if(tmp_pos_score.size() > 2){									// if there are any results (except the two dummy-results added by fastUkkonen)
        for(int unsigned pos=1; pos<(tmp_pos_score.size()-1); pos++){					// then for every result in tmp_pos_score

            // filter to get only the best results in the neighborhood of a hit.
            if(		(tmp_pos_score[pos][0] != tmp_pos_score[pos-1][0] + 1 &&      // if there is no hit in the direct neighborhood
                            tmp_pos_score[pos][0] != tmp_pos_score[pos+1][0] - 1)
                            ||																// OR
                            (tmp_pos_score[pos][0] == tmp_pos_score[pos-1][0] + 1 &&		// if the the hit has two direct neighbors
                            tmp_pos_score[pos][0] == tmp_pos_score[pos+1][0] - 1 &&			// with equal or higher (worse) score --> local score minimum
                            tmp_pos_score[pos][1] <= tmp_pos_score[pos-1][1] &&
                            tmp_pos_score[pos][1] <= tmp_pos_score[pos+1][1]) || !filterResults){ // OR if the filter is turned off

                string B_sequence = sequence.substr(tmp_pos_score[pos][0] - n, n);	// get a m+k long part of the genome sequence from known end-position of the hit
                reverse(B_sequence.begin(), B_sequence.end());				// and reverse it
                string B_read = read;                                                           // make a copy of the read
                reverse(B_read.begin(), B_read.end());                                          // and reverse it

                // run semi-global alignment (normal smith waterman without Ukkonen) for our new sequences (see fastUkkonen for code details)
                int *cj = new int[m+1]; //  g++ does not liked 'int cj[m+1]'
                int cp, cn;
                int min_score = k+1;
                int min_score_pos=-1;

                for(int j=0; j<=m; j++){
                    cj[j] = j;
                }
                for(int i=1; i<=n; i++){
                    cn = 0; cp = 0;
                    for(int j=1; j<=m; j++){
                        if(B_sequence[i-1] == B_read[j-1]){ // match
                            cn = cp;
                        }else if(cp < cn || cj[j] < cn){
                            if(cp < cn){				// mismatch
                                cn = cp;
                                }
                        if(cj[j] < cn){				// gap in read
                            cn = cj[j];
                        }
                        cn++;
                        }else{							// gap in sequence
                            cn++;
                        }
                        cp = cj[j];
                        cj[j] = cn;
                    }
                    if(cj[m] < min_score){				// save the position of the minimal score in bottom row of the DP-matrix
                        min_score = cj[m];
                        min_score_pos = i;
                    }
                }
                delete [] cj;// free reserved space for int-array
                // start_pos of the hit is the known end_pos minus the min_score_pos we just calculated
                int start_pos = tmp_pos_score[pos][0] - min_score_pos;
                int end_pos = tmp_pos_score[pos][0];
                int score = tmp_pos_score[pos][1];
                vector<int> scoreVector(4,0);			// fill a result vector with: (read Nr., start position, end position, score)
                scoreVector[0] = readNr; scoreVector[1] = start_pos; scoreVector[2] = end_pos; scoreVector[3] = score;
                pos_score.push_back(scoreVector);		// save it in pos_score
            }//end if
        }//end for
    }//end if
    cout<<"Procedure 'filterHitsAndBacktrack' for read Nr. "<< readNr<< " done!" <<endl;
}

/* Ukkonen algorithm implemented with one int array and single integer values, see lecture 2 script for pseudo code */
void fastUkkonen(vector<vector<int> > &tmp_pos_score, int &k, string &sequence,int &block_start, int &block_end, string &read, bool &useUkkonenTrick){

	int m,n;
	int lact = k+1;									// initialize last active cell indicator
	m= read.size();									// initialize 'sequence' and 'read' size
	n= block_end+1;

        int *cj=new int[m+1];							// int-vector cj is more or less a column of the DP-matrix. 'vector<int>' was not chosen, becaused it slowed the procedure down
	int cp, cn;                                                                      // single integer values for storage of other DP-matrix cell values

	tmp_pos_score.push_back(vector<int>(2,-5));					// initialize temporary result vector with a "dummy" entry (used for filtering)

	for(int j=0; j<=m; j++){							// initialize our DP-matrix column
            cj[j] = j;
	}

	// Smith-Waterman WITH Ukkonen-trick
	if(useUkkonenTrick){

	for(int i=block_start+1; i<=n; i++){					// for every character of 'sequence'
            cp = 0; cn = 0; 							// set cp and cn zero (correspond to entries [i-1] [0] and [i] [0] of the DP-matrix here

            for(int j=1; j<=lact; j++){					// for characters of 'read' until last active cell
                if(sequence[i-1] == read[j-1]){ // match
                        cn = cp;
                }else{
                    if(cp < cn){				// mismatch
                        cn = cp;
                    }
                    if(cj[j] < cn){				// gap in read
                        cn = cj[j];
                    }
                    cn++;						// else: gap in sequence
                }
                cp = cj[j];
                cj[j] = cn;
            }
            // after calculation of the column actualize the lact indicator
            while(cj[lact] > k){	// reduce lact until last active cell has score lower than k (no need to calculate mor of the column next time)
                lact--;
            }
            if(lact == m){								// if lact == m we have a hit !!
                int pos_end = i;						// end pos of the match is i
                int score = cj[lact];						// score is saved in the last(bottom) cell of cj
                vector<int> tmp_scoreVector(2,0);				// fill a temporary result vector with this info
                tmp_scoreVector[0] = pos_end; tmp_scoreVector[1] = score;
                tmp_pos_score.push_back(tmp_scoreVector);
            }else{										// if lact < m we have to increase it by one
                lact++;
            }
	}
	// Smith-Waterman WITHOUT Ukkonen-trick
	}else{
            for(int i=block_start+1; i<=n; i++){					// for every character of 'sequence'

                cn = 0; cp = 0;						// set cp and cn zero (correspond to entries [i-1] [0] and [i] [0] of the DP-matrix here
                for(int j=1; j<=m; j++){				// for all characters of 'read'
                    if(sequence[i-1] == read[j-1]){ // match
                        cn = cp;
                    }else{
                        if(cp < cn){				// mismatch
                            cn = cp;
                        }
                        if(cj[j] < cn){				// gap in read
                            cn = cj[j];
                        }
                        cn++;						// else: gap in sequence
                    }
                    cp = cj[j];
                    cj[j] = cn;
                }
                if(cj[m] <= k){								// if cj[m] <= k we have a hit !!
                    int pos_end = i;						// end pos of the match is i
                    int score = cj[m];						// score is saved in the last(bottom) cell of cj
                    vector<int> tmp_scoreVector(2,0);				// fill a temporary result vector with this info
                    tmp_scoreVector[0] = pos_end; tmp_scoreVector[1] = score;
                    tmp_pos_score.push_back(tmp_scoreVector);
                }
            }
	}//end else
	tmp_pos_score.push_back(vector<int>(2,-5));		// add another "dummy" entry at the end of the temporary result vector
        delete [] cj; // free reserved space for int-array
}

/* Function to write the results according to the given format */
void writeOutput(vector<vector<int> > &pos_score2){
    ofstream outfile("hits.result"); // creates a file to write in
    string sep = ","; // symbol for separation
    string nl = "\n"; // new line
    int m=pos_score2.size(); // number of rows
    try{
        for (int i = 0; i < m; i++) 
            outfile<<"read_"<<pos_score2[i][0]<<sep<<pos_score2[i][1]<<sep<<pos_score2[i][2]<<sep<<pos_score2[i][3]<<nl;
        outfile.close();
    }catch(exception e){
            cout << "Error during writing!" << endl;
            cout << e.what() << endl;
    }
}
/* write a warning message in the result file according to the given case */
void writeAndTerminate(int termination_case){
    ofstream outfile("hits.result"); // creates a file to write in
    if(termination_case==1){
        outfile<<"read length exceeds block size";
        cout<<"WARNING: 1st abort condition detected => program will terminate..."<<endl;
    }
    else if(termination_case==2){
        outfile<<"Bad choice of k and q for input read length";
        cout<<"WARNING: 2nd abort condition detected => program will terminate..."<<endl;
    }
    outfile.close();
}
/* Given a character (nucleotide base) return assigned number (0-3). If an invalid
 character is introduced -1 is returned */
void convBase(char &input, int &output){
    switch(input){
        case 'A':
            output=0;
            break;
        case 'C':
            output=1;
            break;
        case 'G':
            output=2;
            break;
        case 'T':
            output=3;
            break;
        default:
            output=-1;
    }
}
/* NOT USED
 * Given a number (permutation id) compute the corresponding string(nucleotide 
 * sequence) of length q. With this function one could construct the q-gram table
 * (works up to q=12)*/
string getQGramFromIndex(int id,  int &q){
    string sigma="ACGT";
    bitset< (2*Q) > permutation(id);
    bitset<2> tmp; // bit representation of one nucleotide
    string seq;
    int j;
    for (int i =(2*q-1); i>0; i-=2) {
        tmp[1]=permutation[i];
        tmp[0]=permutation[i-1];

        j=tmp.to_ulong();
        seq.append(sigma,j,1);
    }
    return seq;
}
/* check whether all expected occurrences were found */
void allOccurrencesFound(string &seq, int *nr_occ, int &q, int &nr_rows){
    int counter=0;
    for (int i = 0; i < nr_rows; i++) {
        counter+= nr_occ[i];
    }
    int expected= seq.size()-q+1;
    if(expected==counter)
        cout<<"All occurrences of the q-grams were found!"<<endl;
    else
        cout<<"Warning: expected: "<<expected<<", founded: "<<counter<<endl;
}
/* Scan an input sequence counting the occurrence of the found permutations and 
 saving their positions in a vector */
void getNrOfOccAndPositions(string &seq, int *output, vector<vector<int> > &occ, int &q){
    bitset<2*Q> bit_array;  // all entries equal 0
    int base_id; // nr in the range (0,3)
    int permutation_id; // integer encoding the permutation
    /* read initial (q-1) positions of the sequence, transforming them in bits */
    for (int i = 0; i < (q-1); i++) {
        convBase(seq[i],base_id);
        bit_array=bit_array<<2 |bitset<2*Q>(base_id);
    }
    /* continue the reading of the sequence from position equal q. After the transformation
     * in bits, convert them into integer and use it as index */
    for (unsigned int i = (q-1); i < seq.size(); i++) {
        convBase(seq[i],base_id);
        bit_array=bit_array<<2 |bitset<2*Q>(base_id);
        if(q<Q){        // since Q=12, a neutralization of 2 positions is needed after shifting, if q is smaller than Q
            bit_array.set((q*2),0);
            bit_array.set((q*2+1),0);
        }
        permutation_id=bit_array.to_ulong();    // convert bits in integer
        output[permutation_id]++;               // increase nr. of occurrences of the corresponding permutation
        occ[permutation_id].push_back(i-q+1);   // save start position in a vector
    }
}

/* Fills the table of dimension nr_of_blocks x 2 with the limits of each block
 * according to the block length and the minimal overlap */ 
void calculateBlocksRange(int **blocks_limits,int &seq_length, int &nr_of_blocks, int &block_length,int &min_overlap){
    int expansion_left, expansion_right;
    expansion_right = min_overlap/2;        //same as floor()
    expansion_left = min_overlap - expansion_right;
    // 1st block limits
    blocks_limits[0][0] = 0;
    blocks_limits[0][1] = block_length-1 +expansion_right;
    // limits of the remaining blocks
    for (int i = 1; i < nr_of_blocks; i++) {
        blocks_limits[i][0] = i*block_length - expansion_left;
        blocks_limits[i][1] = (i+1)*block_length-1 + expansion_right;
        }
    // correction/adjust of the ending position of the last block
    blocks_limits[nr_of_blocks-1][1]= seq_length-1;
}
/* similar to getNrOfOccAndPositions in reading/transforming a nucleotide sequence
 into bits. After the transformation of the read sequence it checks, whether the 
 given permutation has an occurrence in the genome (table occ). If yes, then with help
 * of the start position calculates the block(s) containing it and increases their counting number */
void calculateOccPerBlock(string &read, vector<vector<int> > &occ, int &q, int &nr_of_blocks, 
                          int &block_length, int **blocks_range, int *block_counter){
    bitset<2*Q> bit_array;  // all entries equal 0
    int base_id; // nr in the range (0,3)
    int permutation_id, nr_of_occ, q_start,q_end;
    int expected_block_index,block_begin, block_end;

    /* read initial (q-1) positions of the sequence*/
    for (int i = 0; i < (q-1); i++) {
        convBase(read[i],base_id);
        bit_array=bit_array<<2 |bitset<2*Q>(base_id);
    }
    /* continue the reading of the sequence from position equal q. After the transformation
     * in bits, convert them into integer and use it as index */
    for (unsigned int i = (q-1); i < read.size(); i++) {
        convBase(read[i],base_id);
        bit_array=bit_array<<2 |bitset<2*Q>(base_id);
        if(q<Q){
            bit_array.set((q*2),0);
            bit_array.set((q*2+1),0);
        }

        permutation_id=bit_array.to_ulong();
        nr_of_occ=occ[permutation_id].size(); // number of occurrences of the permutation equals the size of the vector containing start positions
        if(nr_of_occ>0){
            for (int j = 0; j < nr_of_occ; j++) {
                q_start= occ[permutation_id][j];
                q_end= q_start+q-1;
                expected_block_index= q_start/block_length;     // block which possibly contains the q-gram
                if(expected_block_index==nr_of_blocks)         // depending on the length of the block, i.e. the minimal length, it could happen
                    expected_block_index--;                             // that the last block contains almost twice elements as the others. A correction of the calculated block index is needed

                /* increase counter of the corresponding block, for simplicity all 3 cases
                are checked: actual block, overlap left and overlap right */
                block_begin=blocks_range[expected_block_index][0];
                block_end=blocks_range[expected_block_index][1];

                if(q_start>=block_begin && q_end<=block_end)        // current block
                    block_counter[expected_block_index]++;

                if(expected_block_index>0){                         // block left
                    block_begin=blocks_range[expected_block_index-1][0];
                    block_end=blocks_range[expected_block_index-1][1];
                    if(q_start>=block_begin && q_end<=block_end)
                        block_counter[expected_block_index-1]++;
                }
                if(expected_block_index<nr_of_blocks-1){            // block right
                    block_begin=blocks_range[expected_block_index+1][0];
                    block_end=blocks_range[expected_block_index+1][1];
                    if(q_start>=block_begin && q_end<=block_end)
                        block_counter[expected_block_index+1]++;
                }
            }
        }
    }
}
/* ########################## MAIN ########################## */
int main(int argc, char**argv) {
    time_int(0); // start timing
    string genome_file, reads_file;
    int k,q,b,threshold,w; 
    b=0; // initialization required by linux
    bool useUkkonenTrick=true; //indicates whether the ukkonen trick will be used or not. 
    bool filterResults = true; // default value for filtering
    /* Print the arguments */
    if (argc > 1) {
        cout << "Welcome..."<<endl << "---- Introduced arguments ----" << endl;
        for (int i = 1; i < argc; i++) {
            switch (i){
                case 1:
                    cout <<"Genome's file:" << "\t" << argv[1] << endl;
                    break;
                case 2:
                    cout <<"Reads file:" << "\t" << argv[2] << endl;
                    break;
                case 3:
                    cout <<"Nr. of allowed errors:" << "\t" << argv[3] << endl;
                    break;
                case 4:
                    cout <<"Length of the q-grams:" << "\t" << argv[4] << endl;
                    break;
                case 5:
                    cout <<"Block size:" << "\t\t" << argv[5] << endl;
                    break;
                default:
                    cout<<"undefined input"<< "\t" <<argv[i]<<endl;
                    break;
            }
        }
        if(argc<6)
            cout<<"WARNING: not enough arguments given!"<<endl;
        /* set given arguments */
        genome_file = argv[1];
        reads_file = argv[2];
        k = atoi(argv[3]); // number of errors
        q = atoi(argv[4]);// length of the q-grams
        b = atoi(argv[5]); // block size for the genome
    } 

    /* Read fasta files*/
    cout <<"---- Reading ----" << endl;
    string genome = readGenome(genome_file);
    int m = getNrOfReads(reads_file); // number of reads 
    vector<string> reads(m, "" );
    readReads(reads_file, reads);

    cout <<"Nr. of reads: "<< m<<endl;
    cout <<"Text's size: "<< genome.size()<<endl;
    cout <<"Pattern's size: "<< reads[0].size()<<endl;

    /* check for 1st and 2nd abort condition */
    w=reads[0].size();
    threshold = w+1-(k+1)*q; //lemma 2 (q-gram)
    cout <<"---- Values for the computation ----" << endl;
    cout<<"Threshold: "<<threshold<<endl;
    if(b<w) //block size is smaller than the read length
        writeAndTerminate(1); // write a warning and terminate (jump to the end)
    else if(threshold<=0)
        writeAndTerminate(2);
    else{
        /* Quasar */
        int nr_of_rows=(int)pow(4,q); //expected nr of permutation of a sequence of length q

        // allocate int array and vector for the scan
        int *nr_of_occ=new int[nr_of_rows](); // initialized with zero
        initializeIntArray(nr_of_occ,nr_of_rows,0); // ensures that the array was initialized with 0
        vector<vector<int> > occurrences(nr_of_rows); // to save the start positions
        cout <<"(Scanning genome w.r.t. q-gram permutations ...)" << endl;
        getNrOfOccAndPositions(genome,nr_of_occ,occurrences,q); // scan the given sequence and save nr of occurrences as well their positions

        /* Counting q-grams, division in blocks, etc */
        int nr_of_blocks,min_block_length,min_overlap,seq_length, read_length;
        seq_length=genome.size();
        read_length=reads[0].size();    // = w
        
        // set values for the minimal block length and for the minimal overlap length
        if(b<(read_length+k))
            min_block_length=read_length+k;
        else
            min_block_length=b;
        min_overlap=read_length+k-1;

        // divide the genome into blocks
        nr_of_blocks=seq_length/min_block_length; // same as floor(). rounded down due to the minimal block length -> last block will contain remaining elements
                
        // block counting
        cout<<"Minimal block length: "<<min_block_length<<endl;
        cout<<"Nr. of blocks: "<<nr_of_blocks<<endl;
        int *blocks_counter=new int[nr_of_blocks]; 
        initializeIntArray(blocks_counter,nr_of_blocks,0);
        int **blocks_range;
        // initialize table 
        blocks_range= new int*[nr_of_blocks]; // # rows (m)
        for (int i = 0; i < nr_of_blocks; i++)
            blocks_range[i] = new int[2]; // # columns (n)
        
        // calculate the limits of each block
        calculateBlocksRange(blocks_range,seq_length,nr_of_blocks,min_block_length,min_overlap);

        /* main process for the computation of m reads */
        //Search for matching q-grams and count matches for each block
        vector<vector<int> > pos_score;					// score vector for format: (nr. of read, start position, end position, score)
        cout <<"---- Starting computation ----" << endl;
        int block_start,block_end;
        for(int readNr = 0; readNr < m; readNr++){				// run ukkonen for readNr many reads, save all scores in pos_score
            calculateOccPerBlock(reads[readNr],occurrences,q,nr_of_blocks,min_block_length,blocks_range,blocks_counter); // scan read and increase block countig if it appears in the genome. See function definition for more details
            for (int i = 0; i < nr_of_blocks; i++) {
                if(blocks_counter[i]>=threshold){       // filtering: only the reads with block counter greater/equal threshold will be verified
                    //cout<<"ReadNr: "<<readNr<<", block nr: "<<i<<" value: "<<blocks_counter[i]<<endl;//@TODO: delete/comment for the final version
                    block_start=blocks_range[i][0];
                    block_end=blocks_range[i][1];
                    vector<vector<int> > tmp_pos_score;         // for the semi-global aligner procedure
                    fastUkkonen(tmp_pos_score, k, genome,block_start,block_end, reads[readNr], useUkkonenTrick);
                    filterHitsAndBacktrack(pos_score, tmp_pos_score, k, genome, reads[readNr], readNr, filterResults);
                }
            }
            initializeIntArray(blocks_counter,nr_of_blocks,0); // reset block counter for the next read
        }
        //allOccurrencesFound(genome,nr_of_occ,q,nr_of_rows);// check whether the expected nr. of occurrences (of q-grams!) was registered
        
        /*  free allocated space */
        // table (2 dimensional array)
        for (int i = 0; i < nr_of_blocks; i++)
            delete[] blocks_range[i];
        delete[] blocks_range;
        // arrays
        delete [] blocks_counter; 
        delete [] nr_of_occ; // may be not necessary since program terminates and the allocated space will be erased anyway 

        /* Ending */
        cout<<"Nr. of occurrences: "<<pos_score.size()<<endl;
        writeOutput(pos_score); // export of the results
    }// end of quasar and semi-global aligning 
    time_int(1); // print out elapsed time
    return 0;
}
