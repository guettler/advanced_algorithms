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

using namespace std;
#define match 0
#define mismatch  1
#define gap 1

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
void firstTerminateWrite(){
    ofstream outfile("hits.result"); // creates a file to write in
    outfile<<"read length exceeds block size";
    outfile.close();
}

/* ########################## MAIN ########################## */
int main(int argc, char**argv) {
    time_int(0); // start timing
    string genome_file, reads_file;
    int k,q,b; 
    bool useUkkonenTrick=true; //indicates whether the ukkonen trick will be used or not. 
    bool filterResults = true; // default value for filtering
//./exercise2 <genome.fasta> <reads.fasta> <k> <q> <b> 
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
    }
    if(argc<6)
        cout<<"WARNING: not enough arguments given!"<<endl;
    /* set given arguments */
    genome_file = argv[1];
    reads_file = argv[2];
    k = atoi(argv[3]); // number of errors
    q=atoi(argv[4]);// length of the q-grams
    b = atoi(argv[5]); // block size for the genome
//    useUkkonenTrick = (atoi(argv[4])==1); // Only value 1 sets it true, all other input values make it false.
    
    /* Read fasta files*/
    cout <<"---- Reading ----" << endl;
    string genome = readGenome(genome_file);
    int m = getNrOfReads(reads_file); // number of reads 
    vector<string> reads(m, "" );
    readReads(reads_file, reads);

    cout <<"Nr. of reads: "<< m<<endl;
    cout <<"Text's size: "<< genome.size()<<endl;
    cout <<"Pattern's size: "<< reads[0].size()<<endl;

    /* check for 1st and 2nd? abort condition */
    if(b<reads[0].size()) //block size is smaller than the read length
        firstTerminateWrite(); // write a warning and terminate (jump to the end)
    else
    {

        /* main process for the computation of m reads */
        vector<vector<int> > pos_score;					// score vector for format: (nr. of read, start position, end position, score)
        cout <<"---- Starting computation ----" << endl;
        for(int readNr = 0; readNr < m; readNr++){				// run ukkonen for readNr many reads, save all scores in pos_score
            vector<vector<int> > tmp_pos_score;
            fastUkkonen(tmp_pos_score, k, genome, reads[readNr], useUkkonenTrick);
            filterHitsAndBacktrack(pos_score, tmp_pos_score, k, genome, reads[readNr], readNr, filterResults);
        }

        /* Ending */
        cout<<"Nr. of occurences: "<<pos_score.size()<<endl;
        writeOutput(pos_score); // export of the results
    }
    time_int(1); // print out elapsed time
    return 0;
}
