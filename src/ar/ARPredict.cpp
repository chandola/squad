#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>

#include "tools.h"
#include "AR.h"

using namespace std;
char *inputfilename=0;
char *modelfilename=0;
char *outputfilename=0;
// int nflag = 0;
#define usage(msg, args...)			\
  {						\
  fprintf(stderr, msg, ##args);			\
  fprintf(stderr,						       \
	  "Usage:\n\n"							\
	  "./ARPredict -i inputfilename -o outputfilename -m modelfilename\n\n"); \
  exit(1);								\
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) inputfilename=strdup(argv[argi]);
      else usage("-i requires input file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) outputfilename=strdup(argv[argi]);
      else usage("-o requires output file argument\n");
    } else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) modelfilename=strdup(argv[argi]);
      else usage("-m requires model file argument\n");
//     } else if(strcmp(argv[argi],"-n")==0){
//       if(++argi<argc) nflag=atoi(argv[argi]);
//       else usage("-n requires normalization argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!inputfilename) usage("No input file specified\n");
  if(!outputfilename) usage("No output file specified\n");
  if(!modelfilename) usage("No model file specified\n");
}

int main(int argc,char **argv)
{
  parseArgs(argc,argv);
  vector<vector<float> > series;
  vector<float> arcoeff(0);
  ifstream in;
  in.open(inputfilename);
  if(in.fail()) {cerr << "input file does not exist "<<inputfilename<<"\n"; exit(0);}
  readSequences(series,in);
  in.close();
//   //subtract mean if the file is unnormalized
//   if(nflag){
//     for(unsigned int i = 0; i < series.size(); i++){
//       float sum=0.0;
//       for(unsigned int j = 0; j < series[i].size(); j++) sum += series[i][j];
//       sum /= series[i].size();
//       for(unsigned int j = 0; j < series[i].size(); j++) series[i][j] -= sum;
//     }
//   }
  // Read the coefficient file
  in.open(modelfilename);
  if(in.fail()) {cerr << "model file does not exist "<<modelfilename<<"\n"; exit(0);}
  vector<vector<float> > tmp;
  readSequences(tmp,in);
  in.close();
  for(unsigned int i = 0; i < tmp.size(); i++)
    arcoeff.push_back(tmp[i][0]);
  unsigned int order = arcoeff.size();
  tmp.clear();tmp.resize(0);
  vector<vector<float> > scores;
  //predict
  for(unsigned int i = 0; i < series.size(); i++){
    vector<float> tmp1(0);
    for(unsigned int j = order; j < series[i].size(); j++){
      float est = 0;
      for (unsigned int k=0;k<order;k++) est += arcoeff[k] * series[i][j-k-1];
      tmp1.push_back((est-series[i][j])*(est-series[i][j]));
    }
    scores.push_back(tmp1);
  }
  ofstream out;
  out.open(outputfilename);
  printSequences(scores,out);
  out.close();
}
