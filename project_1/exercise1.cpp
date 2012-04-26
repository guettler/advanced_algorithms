/* 
 * File:   exercise1.cpp
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
void getMinValue(int &cell, char &nucleotide1, char &nucleotide2,int &val_d, int &val_h, int &val_v);
void semiGlobalWithout(vector<pair<int,int> > &pos_score, int &k, string &sequence, string &read);


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
///* Method without precalculation of the number of reads, misses to save one read */
//void readReads(string &path, vector<string> &reads)
//{
//    ifstream file;
//    string line, tmp;
//
//    int read_nr=-1;
//    try
//    {
//        file.open(path.c_str(), ifstream::in);
//        if(!file)
//        {
//            cout<< "File could not be opened. Please check path and file's name..."<< endl;
//
//        }
//        while(!file.eof())
//        {
//            getline(file, line);
//
//            if(line.at(0)!='>')
//            {
//                tmp.append(line);
//            }
//            else
//            {   
//    
//                if (read_nr!=-1)
//                {                   
//                    reads.push_back(tmp);
//                    tmp="";
//                }                                
//                read_nr++;
//                cout<<read_nr<<endl;
//            }
//        }
//
//        file.close();
//    }
//    catch(exception e)
//    {
//        cout << e.what() << endl;
//    }
//   
//
//}
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
		{
			cout << v[i][j]<< " ";
		}
		cout << endl;
	}
}
void getMinValue(int &cell,char &nucleotide1, char &nucleotide2,int &val_d, int &val_h, int &val_v){
    int diagonal, vertical, horizontal;
    if (nucleotide1==nucleotide2)
        diagonal= val_d +match;
    else
        diagonal= val_d + mismatch;
    vertical=val_v + gap;
    horizontal= val_h + gap;
    
    cell = min(min(horizontal,vertical),diagonal);
}

/* Calculates the scores of the alignment, i.e., last row of the dp-matrix*/
//void semiGlobalWithout(vector<int> &scores, vector<pair<int,int> > &pos_score, int &k, string &sequence, string &read){
void semiGlobalWithout(vector<pair<int,int> > &pos_score, int &k, string &sequence, string &read){
    int m,n;
    m= read.size();
    n= sequence.size();

    /* Dp-matrix: a (m+1)x 2 matrix initialized with 0 */
    vector<vector<int> > dp(m+1,vector<int> (2,0));


    /* Initialize the first column -> initial gap penalty for sequence. The 
     * initial gap penalties for the reads amout to 0 */
    for (int i = 1; i <= m; i++)
        dp[i][0] = dp[i-1][0]+ gap;

//    scores[0]=dp[m][0];
    /* modified Smith-Waterman */
    for (int j = 1; j <= n ; j++) {
        for (int i = 1; i <= m; i++) {
//            dp[i][1] = getMinValue(sequence[j-1], read[i-1],dp[i-1][0], dp[i][0],dp[i-1][1]);
            getMinValue(dp[i][1],sequence[j-1], read[i-1],dp[i-1][0], dp[i][0],dp[i-1][1]);
            
            /* shift one to the right: replace the first value with the second 
             * and write 0 in the second position, same as efect as  adding a 
             * column and deleting the first one */
            dp[i-1][0]=dp[i-1][1]; 
            dp[i-1][1]= 0;

        }


        dp[m][0]=dp[m][1]; // last row has yet to be updated
        dp[m][1]= 0;
        /* save the score after finilized the column */
//        scores[j]=dp[m][0]; // probably we do not need to save all scores, but only the ones <=k, see next line
        
        /* save column position and score if it is less/equal k */
        if (dp[m][0]<=k){
            pos_score.push_back(make_pair(j,dp[m][0]));
//            cout<<"Occurence were found!"<<endl;
        }
//        cout<<"Iteration nr: "<<j <<endl;
    }
    cout<<"Procedure 'semiGlobalWithout' done!" <<endl;
}
/* ########################## MAIN ########################## */
int main(int argc, char**argv) {
    time_int(0); // start timing
    string genome_file, reads_file;
    int k, ukkonen_on;
    
    // Prints welcome message...
    cout << "Welcome ..." << endl;
    
    /* --- read the arguments. To be uncomment for the final version ----  */
    // Prints arguments...
//    if (argc > 1) {
//        cout << endl << "Arguments:" << endl;
//        for (int i = 1; i < argc; i++) {
//            cout << i << ": " << argv[i] << endl;
//        }
//    }
//    genome_file = argv[1];
//    reads_file = argv[2];
//    k = atoi(argv[3]); // number of errors
//    ukkonen_on = atoi(argv[4]); // indicates whether the ukkonen trick will be used or not
    /* --------------------- end of block ----------------------------------*/
    
    /* ------ Input block for the working phase.------------------------------ 
     * -------To be erased before checking -----*/
    string fileNames[]={"random10M.fasta","random10M_reads50_100.fasta","random10M_reads50_1k.fasta",
    "random10M_reads100_100.fasta","random10M_reads100_1k.fasta",
    "random10M_reads400_100.fasta","random10M_reads400_1k.fasta"}; // 7 given test files
    
    genome_file=fileNames[0];
    reads_file=fileNames[6];
    k=3;
    ukkonen_on=0;
    /* ------------------------------------------------------------------------*/

    
    
    /* Read fasta files*/
    string genome = readGenome(genome_file);
    int m = getNrOfReads(reads_file);
    vector<string> reads(m, "" );
//    vector<string> reads;// other way without m (does not work perfectly...)
    readReads(reads_file, reads);
    
    cout <<"Nr. of Reads: "<< m<<endl;
    cout <<"1st sequence's size: "<< genome.size()<<endl;
    cout <<"2nd sequence's size: "<< reads[0].size()<<endl;
 
    
    /* ################ test ################################################*/
//    cout<< reads.size()<<endl;
//    cout<< reads[0]<<endl;
    /* Vector for the scores: probably not needed */
//    vector<int> scores (genome.size()+1,0);
//    vector<int> scores (reads[1].size()+1,0);
    
    /* Vector to save the pairs: column position and respective score*/
    vector<pair<int,int> > pos_score;
    
    semiGlobalWithout(pos_score, k,reads[1],reads[1]);// match read with itself, it works.


//     semiGlobalWithout(pos_score, k,genome,reads[0]); // needed 31 minutes with netbeans, but 35 seconds under linux
    

    
    cout<<"Nr. of occurences: "<<pos_score.size()<<endl;
//    cout<<pos_score.back().first <<pos_score.back().second <<endl;
    time_int(1); // print out elapsed time
    return 0;
}
