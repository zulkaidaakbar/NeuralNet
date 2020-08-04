#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>


using namespace std; 



double temp[50000];

void write_text()
{


char name[100];
char name1[100];
int epoch = 20000;


ifstream inputfile;
sprintf(name,"data_july14c.txt"); 
inputfile.open(name,ios::in);
cout<<" Get the data from "<<name<<endl<<endl;


  
for (int i = 0; i < 9720 ; i++) 
 {
  inputfile >> temp[i];
 } 

//
ofstream myfile;
myfile.open ("rnd_test2.txt");
int top1= 8;
int top2= 15;
int top3= 3;
myfile<<"topology:"<<" "<<top1<<" "<<top2<<" "<<top3<<endl;
for (int mm = 1; mm<2; mm++)
 {
  for (int oo=0; oo<epoch; oo++)
    {
     for (int nn=0; nn<36; nn++)
     {
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(temp[mm*648+nn*18+7],0.05*temp[mm*648+nn*18+7]);
   
      myfile<<"in:"<<" "<<temp[mm*648+nn*18+2]<<" "<<temp[mm*648+nn*18+3]<<" "<<temp[mm*648+nn*18+4]<<" "<<temp[mm*648+nn*18+5]<<" "<<temp[mm*648+nn*18+6]<<" "<<temp[mm*648+nn*18+9]
<<" "<<temp[mm*648+nn*18+10]<<" "<<temp[mm*648+nn*18+11]<<endl;
      myfile<<"out:"<<" "<<distribution(generator)<<" "<<distribution(generator)<<" "<<distribution(generator)<<endl;
    }
   
  }
 }


}
