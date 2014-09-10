#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>

#include "tools.h"
#include "AR.h"

#define MAXCOEFF 1000

using namespace std;
char *inputfilename=0;
char *modelfilename=0;
unsigned int order = 0;
int method = 0;
#define usage(msg, args...)			\
  {						\
  fprintf(stderr, msg, ##args);			\
  fprintf(stderr,						       \
	  "Usage:\n\n"							\
	  "./ARTrain -i inputfilename -m modelfilename -d order [-l method)]\n\n" \
	  "-i      input file name\n"\
	  "-m      model file name \n"\
	  "-d      order of AR model [1 .. 100]\n"\
	  "-l      method\n"\
	  "        0 Least Squares (default)\n"\
	  "        1 Max Entropy\n"\
	  "\n"); \
  exit(1);								\
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) inputfilename=strdup(argv[argi]);
      else usage("-i requires input file argument\n");
    } else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) modelfilename=strdup(argv[argi]);
      else usage("-m requires model file argument\n");
    } else if(strcmp(argv[argi],"-d")==0){
      if(++argi<argc) order=atoi(argv[argi]);
      else usage("-d requires order argument\n");
    } else if(strcmp(argv[argi],"-l")==0){
      if(++argi<argc) method=atoi(argv[argi]);
      else usage("-l requires method argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!inputfilename) usage("No input file specified\n");
  if(!modelfilename) usage("No model file specified\n");
  if(order == 0) usage("Order of AR model not specified\n");
  if(method == -1) usage("Method not specified\n");
}

int main(int argc,char **argv)
{
  parseArgs(argc,argv);

  vector<vector<float> > data;
  double *coefficients,*series;
  if (order >= MAXCOEFF) {
    fprintf(stderr,"Maximum degree is %d\n",MAXCOEFF-1);
    exit(-1);
  }
  
   if (method == 0) {
     method = MAXENTROPY;
   } else if (method == 1) {
     method = LEASTSQUARES;
   } else {
     fprintf(stderr,"Didn't get a valid method\n");
     exit(-1);
   }
   
   ifstream in;
   in.open(inputfilename);
   if(in.fail()) {cerr << "input file does not exist "<<inputfilename<<"\n"; exit(0);}
   readSequences(data,in);
   in.close();
   int length = 0;
   for(unsigned int i = 0; i < data.size(); i++) length += data[i].size();
   //prepare an array for the AutoRegression function
   series = (double *) malloc(sizeof(double)*length);
   unsigned int k = 0;
   for(unsigned int i = 0; i < data.size(); i++){
     for(unsigned int j = 0; j < data[i].size(); j++){
       series[k] = data[i][j];
       k++;
     }
   }
   coefficients = (double *) malloc(sizeof(double)*order);
   
   // Calculate and print the coefficients
   if (!AutoRegression(series,length,order,coefficients,method)) {
     cerr << "AR routine failed\n";
     exit(-1);
   }
   //print out the AR model to a file
   ofstream out(modelfilename);
   for (unsigned int i=0;i<order;i++) out<<coefficients[i]<<"\n";
   out.close();
}
