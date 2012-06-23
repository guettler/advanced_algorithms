#include<iostream>
#include<string>
#include<vector>
using namespace std;;
using std::string;


int main(){

void move2front(string); //prototype
void m2f_decoding(int r[],string);


//Testing ===============//
string s ="ipssm$pissii"; //teststring
int R[] ={0,3,0,0,0,1,0,3,0,2,2,3,3,1,0,1,0};
//string s ="aooooaaiioaeieeii";
string Alp ="aeio"; //should be in input file
//move2front(s); //function calls
m2f_decoding(R,Alp);
    

system("PAUSE");
return(0);
} //end main


//====== functions ===========================//
//============================================//

void move2front(string L)
{
int laenge= L.length();
string Y=""; //dummy string

//========= preprocessing ===================//

   for(int i = 0;i <= L.length()-1; i++) // scan text for chars with growing dummy string
   { 
             int found =Y.find(L.at(i));
             if(found == -1){Y=Y+L.at(i);} //-1 for not in string                                    
   }
   int sigma = Y.length();
   //std::sort(Y.begin(), Y.end()); // sort Y alphabetically if needed
   cout<< "Die Zeichenkette beinhaltet folgende Buchstaben: " << Y << endl;
   cout<< "Die Laenge betraegt: " << sigma << endl;       
    
//========== m2f ============================//


int R[laenge-1]; //R-table

    for(int i=0;i<=laenge-1;i++)
    {
    R[i] = Y.find(L.at(i));
         //cout <<"Anfang: " << Y<<endl;    
         if(R[i]!=0){ 
         Y = Y.erase(R[i],1); //delete one char from Y at pos R[i]
         string temp;
         temp.push_back(L[i]); // make dummy string for insertion (cause of conversion problems O_o)
         //cout << temp;
         Y = Y.insert(0,temp); //copy deleted char at begining
         temp =""; // reset temp
         }
    //cout << i<< "-te "<<"Schleife: "  << Y << endl;
    cout << R[i];
            
    }

 
} //end m2f

//========== m2f decoding ====================//

void m2f_decoding(int r[],string alphabet)
{
string L;

int laenge=17; //length of R, determine in main function or m2f_decoding

       for(int i=0;i<=laenge-1;i++)
       {
       L.push_back(alphabet.at(r[i]));
            if(r[i]!=0){
            string temp;
            temp.push_back(alphabet.at(r[i]));
            alphabet = alphabet.erase(r[i],1); //delete one char from alphabet at pos R[i]
            alphabet = alphabet.insert(0,temp); //copy deleted char at begining
            temp =""; // reset temp
            }
       }

cout << L;

} //end m2f decoding





