/* 
 * File:   exercise1.cpp
 * Author: tbd
 */



#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/time.h>

using namespace std;

/* Function prototypes (instead of header file)*/
string readGenome(string &path);
int getNrOfReads(string &path);


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
void printTable(vector<vector<int> > &v)
{
	int length = v.size();
	int width = v[0].size();

	for(int i=0; i<length; ++i)
	{
		for(int j=0; j<width; ++j)
		{
			cout << v[i][j];
		}
		cout << endl;
	}
}
int getMaxValue(char &nucleotide1, char &nucleotide2,int &match,int &mismatch, int &gap,int &val_d, int &val_h, int &val_v){
    int diagonal, vertical, horizontal,cell;
    if (nucleotide1==nucleotide2)
        diagonal= val_d +match;
    else
        diagonal= val_d + mismatch;
    vertical=val_v + gap;
    horizontal= val_h + gap;
    
    cell = max(max(horizontal,vertical),diagonal);
    return(cell);
}

/* Calculates the scores of the alignment, i.e., last row of the dp-matrix*/
void semiGlobalWithout(vector<int> &scores, string &sequence, string &read){
    int m,n,match,mismatch,internal_gap,initial_gap;
    m= read.size();
    n= sequence.size();
    match=1;
    mismatch=-1;
    initial_gap=-1;
    internal_gap=-2;
    /* Dp-matrix: a (m+1)x2 matrix initialized with 0 */
    vector<vector<int> > dp(m+1,vector<int> (2,0));


    /* Initialize the first column -> initial gap penalty for sequence. The 
     * initial gap penalties for the reads amout to 0 */
    for (int i = 1; i <= m; i++)
        dp[i][0] = dp[i-1][0]+ initial_gap;

    scores[0]=dp[m][0];
    /* modified Smith-Waterman */
    for (int j = 1; j <= n ; j++) {
        for (int i = 1; i <= m; i++) {
            dp[i][1] = getMaxValue(sequence[j-1], read[i-1],match, mismatch, internal_gap, dp[i-1][0], dp[i][0],dp[i-1][1]);
            
            /* shift one to the right: replace the first value with the second 
             * and write 0 in the second position, same as efect as  adding a 
             * column and deleting the first one */
            dp[i-1][0]=dp[i-1][1]; 
            dp[i-1][1]= 0;

        }

        scores[j-1]=dp[m][1];
        dp[m][0]=dp[m][1]; // last row has yet to be updated
        dp[m][1]= 0;
    }
    scores[n+1]=dp[m][1];

}
/* ########################## MAIN ########################## */
int main(int argc, char**argv) {
    time_int(0); // start timing
    // Prints welcome message...
    cout << "Welcome ..." << endl;

    // Prints arguments...
    if (argc > 1) {
        cout << endl << "Arguments:" << endl;
        for (int i = 1; i < argc; i++) {
            cout << i << ": " << argv[i] << endl;
        }
    }
    string genome_file = argv[1];
    string reads_file = argv[2];
    int k = atoi(argv[3]); // number of errors
    int ukkonen_on = atoi(argv[4]); // indicates wheter the ukkonen trick will be used or not
    
//    cout << genome_file[2]<< endl;
    
    string genome = readGenome(genome_file);
    int m = getNrOfReads(reads_file);
    cout << m<<endl;

    vector<string> reads(m, "" );
//    vector<string> reads;// other way without m (does not work perfectly...)
    readReads(reads_file, reads);
  
    
    time_int(1); // print out elapsed time
    
    /* ################ test ################################################*/
//    cout<< reads.size()<<endl;
//    cout<< reads[0]<<endl;
    /* Vector for the scores */
//    vector<int> scores (reads[1].size()+1,0);
//    semiGlobalWithout(scores,reads[1],reads[1]);
    
    string seq1="ABCD";
    vector<int>scores (seq1.size()+1,0);
    semiGlobalWithout(scores,seq1,seq1);
    cout<<scores[1]<<endl;
    return 0;
}
