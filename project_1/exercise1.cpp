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
#include <map>

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
void getMinValueAndPredecessor(int &cell,char &predecessor, char &nucleotide1, char &nucleotide2,int &val_d, int &val_h, int &val_v);
void semiGlobalWithout(vector<pair<int,int> > &pos_score, map<pair<int,int>, char> &traces,int &k, string &sequence, string &read);


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
void printMapTable(map<pair<int,int>,char> &m)
{
	map<pair<int,int>,char>::iterator it;
	//int length = m.size();
	
  for ( it=m.begin() ; it != m.end(); it++ )// in part taken from http://www.cplusplus.com/reference/stl/map/clear/
    cout << (*it).first.first<<","<<(*it).first.second<< " => " << (*it).second << endl;
}


//void getMinValue(int &cell,char &nucleotide1, char &nucleotide2,int &val_d, int &val_h, int &val_v){
//    int diagonal, vertical, horizontal;
//    if (nucleotide1==nucleotide2)
//        diagonal= val_d +match;
//    else
//        diagonal= val_d + mismatch;
//    vertical=val_v + gap;
//    horizontal= val_h + gap;
//    
//    cell = min(min(horizontal,vertical),diagonal);
//}
void getMinValueAndPredecessor(int &cell,char &predecessor,char &nucleotide1, char &nucleotide2,int &val_d, int &val_h, int &val_v){
    int diagonal, vertical, horizontal, h_vs_v, h_vs_v_vs_d;
    if (nucleotide1==nucleotide2)
        diagonal= val_d +match;
    else
        diagonal= val_d + mismatch;
    vertical=val_v + gap;
    horizontal= val_h + gap;
        /* two ways @TODO: after runtime comparison one of them has to be erased*/ 
//    // 1st
//    if (horizontal <= diagonal)
//    {
//        cell=horizontal;
//        predecessor="h";
//    }else{
//        cell=diagonal;
//        predecessor="d";
//    }
//    if (horizontal <= vertical)
//    {
//        cell=horizontal;
//        predecessor="h";
//    }else{
//        cell=vertical;
//        predecessor="v";
//    }
//    if (diagonal <= vertical)
//    {
//        cell=diagonal;
//        predecessor="d";
//    }else{
//        cell=vertical;
//        predecessor="v";
//    }
    
    // 2nd
    h_vs_v=min(horizontal,vertical);
    h_vs_v_vs_d=min(h_vs_v,diagonal);
    if(h_vs_v_vs_d==diagonal)
    {
        cell=diagonal;
        predecessor='d';// 
    }else if (h_vs_v==horizontal) {
        cell=horizontal;
        predecessor='h'; // horizonal, i.e., left from current cell
    }else
    {
        cell=vertical;
        predecessor='v';// vertical, i.e., above the current cell
    }
/* Remark: due to the order of comparisions the default resulting assigned 
 * values are:
 * - If diagonal=horizontal=vertical -> diagonal 
 * - if horizontal=vertical, both < diagonal -> horizontal */

    
}

/* Calculates the scores of the alignment, i.e., last row of the dp-matrix*/
//void semiGlobalWithout(vector<int> &scores, vector<pair<int,int> > &pos_score, int &k, string &sequence, string &read){
void semiGlobalWithout(vector<pair<int,int> > &pos_score, map<pair<int,int>, char> &traces,int &k, string &sequence, string &read){
    int m,n;
    char predecessor;
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
            getMinValueAndPredecessor(dp[i][1],predecessor,sequence[j-1], read[i-1],dp[i-1][0], dp[i][0],dp[i-1][1]);
            if(dp[i][1]<=k)
                traces.insert(pair<pair<int,int>,char >(make_pair(i,j),predecessor));//save coordinates and maximum's type
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
/* Function to write the results. DRAFT: to be completed according the its argument*/
void writeOutput()
{
    ofstream outfile("hits.result"); // creates a file to write in

    string sep = ","; // symbor for separtion
    string nl = "\n"; // new line

    try
    {
        string output = ""; 
//      outfile<<"<id>, <start>, <end>, <errors>"<<nl; // headline
        
        for (int i = 0; i < 10; i++) // example
        {
        outfile<<"read_"<<i<<sep<<"start_pos"<<sep<<"end_pos"<<sep<<"errors"<<nl;


        }

        outfile << output << endl;
        outfile.close();

    }catch(exception e)
    {
            cout << "Error during writing!" << endl;
            cout << e.what() << endl;
    }
}
/* Backtracking */
int getStartPosition(int final_row, int final_column, map<pair<int,int>, char> &traces)
{
    int i,j;
    i=final_row;
    j= final_column;
    char predecessor;
    while (i >1 ) {
        predecessor=traces.find(make_pair(i,j))->second;
        switch(predecessor)
        {
            case 'd':
                i--;
                j--;
                break;
            case 'h':
                j--;
                break;
            case 'v':
                i--;
                break;
        }


    }
    return(j); // return only starting column since starting row is always 1
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
    k=2;
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
    map<pair<int,int>,char> traces;
    
    
    semiGlobalWithout(pos_score,traces, k,reads[1],reads[1]);// match read with itself, it works.



//     semiGlobalWithout(pos_score,traces, k,genome,reads[0]); // needed 31 minutes with netbeans, but 35 seconds under linux
    
    cout<<"Nr. of occurences: "<<pos_score.size()<<endl;
    /* testing backtracking (example lecture)*/
//    string seq1="ANNEALING";
//    string seq2="ANNUAL";
//    semiGlobalWithout(pos_score,traces, k,seq1,seq2);
//    printMapTable(traces);
//    int pos = getStartPosition(6,5,traces);
    int pos = getStartPosition(400,400,traces);
    cout<<pos<<endl;
    time_int(1); // print out elapsed time
    return 0;
}
