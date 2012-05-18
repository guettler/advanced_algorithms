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
void fastUkkonen(vector<vector<int> > &tmp_pos_score, int &k, string &sequence, string &read, bool &useUkkonenTrick);
void filterHitsAndBacktrack(vector<vector<int> > &pos_score, vector<vector<int> > &tmp_pos_score,
		int &k, string &sequence, string &read, int &readNr, bool &filterResults);
void writeOutput(vector<vector<int> > &pos_score2);
void firstTerminateWrite();
/* @TODO: complete section with the declaration of new functions */


/* Function definitions */
string readGenome(string &path)
{
    ifstream file;
    string line, input;
    try
    {
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
    catch(exception e)
    {
        cout << e.what() << endl;
    }
    file.close();
    return input;
}

int getNrOfReads(string &path)
{
    ifstream file;
    string line;
    int nr_of_reads=0;
    try
    {
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

void readReads(string &path, vector<string> &reads)
{
    ifstream file;
    string line;
    int read_nr=-1;
    try
    {
        file.open(path.c_str(), ifstream::in);
        if(!file)
            cout<< "File could not be opened. Please check path and file's name..."<< endl;
        while(!file.eof())
        {
            getline(file, line);
            if(line.at(0)!='>')
                reads[read_nr].append(line);
            else
                read_nr++;
        }
        file.close();
    }
    catch(exception e)
    {
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

void printVector(vector<int>  &v)
{
    int length = v.size();
    for(int i=0; i<length; ++i)
        cout << v[i]<< " ";
    cout << endl;

}
void printTable(vector<vector<int> > &v)
{
    int length = v.size();
    int width = v[0].size();

    for(int i=0; i<length; ++i)
    {
        for(int j=0; j<width; ++j)
                cout << v[i][j]<< " ";
        cout << endl;
    }
}

/* Filter the results of fastUkkonen and do backtracking for remaining hits */
void filterHitsAndBacktrack(vector<vector<int> > &pos_score, vector<vector<int> > &tmp_pos_score,
		int &k, string &sequence, string &read, int &readNr, bool &filterResults){

    
    int m = read.size();
    int n = m + k;			// we only calculate an m+k wide part of the dp-matrix for backtracking (worst case: k gaps in the read)

    // filter and backtracking
    if(tmp_pos_score.size() > 2){									// if there are any results (exept the two dummy-results added by fastUkkonen)
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

                // run semiglobal alignment (normal smith waterman without ukkonen) for our new sequences (see fastUkkonen for code details)
                int *cj = new int[m+1]; //  g++ doesnt liked 'int cj[m+1]'
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
void fastUkkonen(vector<vector<int> > &tmp_pos_score, int &k, string &sequence, string &read, bool &useUkkonenTrick){

	int m,n;
	int lact = k+1;									// initialize last active cell indicator
	m= read.size();									// initialize 'sequence' and 'read' size
	n= sequence.size();

        int *cj=new int[m+1];							// int-vector cj is more or less a column of the DP-matrix. 'vector<int>' was not chosen, becaused it slowed the procedure down
	int cp, cn;                                                                      // single integer values for storage of other DP-matrix cell values

	tmp_pos_score.push_back(vector<int>(2,-5));					// initialize temporary result vector with a "dummy" entry (used for filtering)

	for(int j=0; j<=m; j++){							// initialize our DP-matrix column
            cj[j] = j;
	}

	// Smith-Waterman WITH Ukkonen-trick
	if(useUkkonenTrick){

	for(int i=1; i<=n; i++){					// for every character of 'sequence'
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

            for(int i=1; i<=n; i++){					// for every character of 'sequence'

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
void writeOutput(vector<vector<int> > &pos_score2)
{
    ofstream outfile("hits.result"); // creates a file to write in
    string sep = ","; // symbor for separtion
    string nl = "\n"; // new line
    int m=pos_score2.size(); // number of rows
    try
    {
        for (int i = 0; i < m; i++) 
            outfile<<"read_"<<pos_score2[i][0]<<sep<<pos_score2[i][1]<<sep<<pos_score2[i][2]<<sep<<pos_score2[i][3]<<nl;
        outfile.close();
    }catch(exception e)
    {
            cout << "Error during writing!" << endl;
            cout << e.what() << endl;
    }
}
/* write a warning message in the result file according to the given case */
void writeAndTerminate(int termination_case){
    ofstream outfile("hits.result"); // creates a file to write in
    if(termination_case==1)
    {
        outfile<<"read length exceeds block size";
        cout<<"WARNING: 1st abort condition detected => program will terminate..."<<endl;
    }
    else if(termination_case==2)
    {
        outfile<<"Bad choice of k and q for input read length";
        cout<<"WARNING: 2nd abort condition detected => program will terminate..."<<endl;
    }
    outfile.close();
}
void convBase(char &input, int &output)
{
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
            break;
    }
}
/* Given a number (permunation id) compute the corresponding string(nucleotide 
 * sequence) of length q
 * works up to q=12*/
string getQGramFromIndex(int id,  int &q)
{
    string sigma="ACGT";
    bitset< (2*Q) > permutation(id);
//    cout<<"permuataion: "<<permutation<<endl;
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
/* check whether all expected occurences were found */
void allOccurencesFound(string &seq, int nr_occ[],int &q, int &nr_rows)
{
    int counter=0;
    for (int i = 0; i < nr_rows; i++) {
        counter+= nr_occ[i];
    }
    int expected= seq.size()-q+1;
    if(expected==counter)
        cout<<"All occurences were found!"<<endl;
    else
        cout<<"Warning: expected: "<<expected<<", founded: "<<counter<<endl;

}
/* Scan an input sequence counting the occurrence of the found permutations and 
 saving their positions in a vector */
void getNrOfOccAndPositions(string &seq, int *output, vector<vector<int> > &occ, int &q){
    
    bitset<2*Q> bit_array;  // all entries equal 0
    int base_id; // nr in the range (0,3)
    int permutation_id;
    /* read initial q positions of the sequence*/
    for (int i = 0; i < q; i++) {
        convBase(seq[i],base_id);
        bit_array=bit_array<<2 |bitset<2*Q>(base_id);
    }
    permutation_id=bit_array.to_ulong();
    output[permutation_id]++; // increment the occurence of the given q-gram
    occ[permutation_id].push_back(0);
//    cout<<permutation_id<<endl;
    
    /* same as above, for the rest of the sequence */
    for (unsigned int i = q; i < seq.size(); i++) {
        convBase(seq[i],base_id);
        bit_array=bit_array<<2 |bitset<2*Q>(base_id);
//        cout<<"i: "<<i<< "bit_array: "<<bit_array<<endl;
        if(q<Q)
        {
            bit_array.set((q*2),0);
            bit_array.set((q*2+1),0);
//            cout<<"bitarray: "<<bit_array<<endl;
        }
        permutation_id=bit_array.to_ulong();
//        cout<<"Permutation id :"<<permutation_id<<endl;
        output[permutation_id]++;
        occ[permutation_id].push_back(i-q+1);
//        cout<<"Output at: "<<output[permutation_id]<<endl;
    }
}
/* Initializes an integer array with a given value. Implemented since the array
 * initializations with constructor are not reliable */
void initializeIntArray(int *arr, int &nr_of_elements, int value)
{
    for (int i = 0; i < nr_of_elements; i++)
        arr[i] = value;
}
void printIntArray(int v[], int size)
{
    for(int i=0; i< size; ++i)
        cout << "["<<i<<"] "<<v[i]<< " "<<endl;
}
/* ########################## MAIN ########################## */
int main(int argc, char**argv) {
    time_int(0); // start timing
    string genome_file, reads_file;
    int k,q,b,threshold,w; 
    bool useUkkonenTrick=true; //indicates whether the ukkonen trick will be used or not. 
    bool filterResults = true; // default value for filtering
// format: ./exercise2 <genome.fasta> <reads.fasta> <k> <q> <b> 
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
                    cout <<"Nr. of allowed erros:" << "\t" << argv[3] << endl;
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
    /* ------ Input block for the working phase.------------------------------        
     * -------To be erased before checking -----*/
    else if(argc==1) // if no arguments are given by running the progam. 
    {
        string fileNames[]={"random10M.fasta","random10M_reads50_100.fasta","random10M_reads50_1k.fasta",
        "random10M_reads100_100.fasta","random10M_reads100_1k.fasta",
        "random10M_reads400_100.fasta","random10M_reads400_1k.fasta"}; // 7 given test files

        genome_file=fileNames[0];
        reads_file=fileNames[6];
        k=0;
        q= 0;
        b=400;    
    }
    /* ------------------------------------------------------------------------*/
   
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
    threshold = w+1-(k+1)*q; //lemma 2 (q-gramm)
    cout<<"Computed threshold: "<<threshold<<endl;
    if(b<w) //block size is smaller than the read length
        writeAndTerminate(1); // write a warning and terminate (jump to the end)
    else if(threshold<=0)
        writeAndTerminate(2);
   
    else
    {
        /* MAIN */
        int nr_of_rows=(int)pow(4,q); //expected nr of permutation of a sequence of length q

        // allocate int array and vector for the scan
        int *nr_of_occ=new int[nr_of_rows](); // initialized with zero
        initializeIntArray(nr_of_occ,nr_of_rows,0); // ensures that the array was initialized with 0
        vector<vector<int> > occurrences(nr_of_rows); // to save the start positions
        
        getNrOfOccAndPositions(genome,nr_of_occ,occurrences,q); // scan the given sequence and save nr of occurrences as well their positions
        
        
        /* @TODO: counting q-grams, blocking, etc*/
        
        
        
        /* main process for the computation of m reads */
        // @TODO: to be adapted according to: "For each block where the matches 
        // exceed the threshold (q-gram Lemma) use the semi-global aligner from Exercise 1 for verification." 
//        vector<vector<int> > pos_score;					// score vector for format: (nr. of read, start position, end position, score)
//        cout <<"---- Starting computation ----" << endl;
//        for(int readNr = 0; readNr < m; readNr++){				// run ukkonen for readNr many reads, save all scores in pos_score
//            vector<vector<int> > tmp_pos_score;
//            fastUkkonen(tmp_pos_score, k, genome, reads[readNr], useUkkonenTrick);
//            filterHitsAndBacktrack(pos_score, tmp_pos_score, k, genome, reads[readNr], readNr, filterResults);
//        }
//
//        /* Ending */
//        cout<<"Nr. of occurences: "<<pos_score.size()<<endl;
//        writeOutput(pos_score); // export of the results
        
        allOccurencesFound(genome,nr_of_occ,q,nr_of_rows);
        //delete [] nr_of_occ; // may be not necessary since program terminates and the allocated space will be erased anyway 
    }
    time_int(1); // print out elapsed time
    return 0;
}
