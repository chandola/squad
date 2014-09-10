#include <iostream>
#include <fstream>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include "svm.h"
#include "sax.h"
#include "tools.h"

using namespace std;
char *infilename;
char *modelfilename;
char *outputfilename;
int len = 0;
int hop = 1;
int segment = 1;
struct svm_model* model;
#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n\n" \
"./SVRPredict -i test_sequence_filename -m model_filename -o output_filename -l length [-x hop -s segment length]\n\n" \
"   -i        -> input file\n" \
"   -m        -> SVM model file\n" \
"   -o        -> output results file\n" \
"   -l        -> length of window\n" \
"   -x        -> number of hops (default 1)\n"\
"   -s        -> segment length (default 1)\n"\
"\n");\
  exit(1); \
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) infilename=strdup(argv[argi]);
      else usage("-i requires input file argument\n");
    } else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) modelfilename=strdup(argv[argi]);
      else usage("-m requires model file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) outputfilename=strdup(argv[argi]);
      else usage("-o requires output file argument\n");
    }else if(strcmp(argv[argi],"-l")==0){
      if(++argi<argc) len=atoi(argv[argi]);
      else usage("-l requires measure argument\n");
    }else if(strcmp(argv[argi],"-x")==0){
      if(++argi<argc) hop=atoi(argv[argi]);
      else usage("-x requires hop argument\n");
    }else if(strcmp(argv[argi],"-s")==0){
      if(++argi<argc) segment=atoi(argv[argi]);
      else usage("-s requires segment argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!infilename) usage("No input sequence filename specified\n");
  if(!modelfilename)   usage("No model filename specified\n");
  if(!outputfilename)   usage("No output filename specified\n");
  if(len==0)     usage("Window length not specified\n");
}

int main (int argc, char *argv[]){ 
  
  parseArgs(argc,argv);

  //read SVM Model File
  if((model=svm_load_model(modelfilename))==0){
    fprintf(stderr,"can't open model file %s\n",modelfilename);
    exit(1);
  }

  // Read testing set
  vector<vector<float> > ts,ts_segs;
  ifstream in(infilename);
  readSequences(ts,in);
  in.close();
  if(segment > 1) TStoPAA(ts,ts_segs,segment);
  else ts_segs = ts; 
  vector<vector<float> > windows;
  vector<unsigned int> lens;getWindows_hop(ts_segs,windows,len,hop,lens);
  if(windows.size() == 0){
    cerr<<"Error: No windows created.\n";
    exit(0);
  }
  ts.clear();ts_segs.clear();

  //predict each window
  struct svm_node *x;
  vector<float> scores;
  for(unsigned int i = 0; i < windows.size(); i++){
    x = (struct svm_node *) malloc((windows[i].size()+1)*sizeof(struct svm_node));
    for(unsigned int j = 0; j < windows[i].size(); j++){
      x[j].value = (double) windows[i][j];
      x[j].index = (int) j+1;
    }
    x[windows[i].size()].index = -1;
    double pred = svm_predict(model,x);
    scores.push_back(pred);
  }
  svm_destroy_model(model);

  vector<vector<float> > scores1;
  int currind = 0;
  for(unsigned int i = 0; i < lens.size();i++){
    vector<float>tmp(scores.begin()+currind,scores.begin()+currind+lens[i]);
    currind += lens[i];
    scores1.push_back(tmp);
  }
  ofstream out;
  out.open(outputfilename);
  printSequences(scores1,out);
  out.close();
}
