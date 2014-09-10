//winc.cpp - Varun Chandola

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <queue>
#include <sys/stat.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "tools.h"

using namespace std;

char *trainseqfile = 0;//file containing input training sequences
char *testseqfile = 0;//file containing input test sequences
char *outputfile = 0;//file to which results will be written
int n = 2;//size of the window
int l = 1;//sliding step
int nn = 10;//nearest neighbors
#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n\n" \
"./WINC -i test_sequence_filename -t train_sequence_filename -o outputfile [-n window -l hop -nn nearestneighbors]\n\n"); \
  exit(1); \
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) testseqfile=strdup(argv[argi]);
      else usage("-i requires test file argument\n");
    } else if(strcmp(argv[argi],"-t")==0){
      if(++argi<argc) trainseqfile=strdup(argv[argi]);
      else usage("-t requires training file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) outputfile=strdup(argv[argi]);
      else usage("-o requires output file argument\n");
    } else if(strcmp(argv[argi],"-n")==0){
      if(++argi<argc) n=atoi(argv[argi]);
      else usage("-n requires window length argument\n");
    }else if(strcmp(argv[argi],"-l")==0){
      if(++argi<argc) l=atoi(argv[argi]);
      else usage("-l requires hop length size argument\n");
    }else if(strcmp(argv[argi],"-nn")==0){
      if(++argi<argc) nn=atoi(argv[argi]);
      else usage("-nn requires nearest neighbors argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!testseqfile)     usage("No test file specified\n");
  if(!trainseqfile)     usage("No train file specified\n");
  if(!outputfile)   usage("No output file specified\n");
  if(n <= 0) usage("Invalid value for n. Should be greater than equal to 1.");
  if(l <= 0) usage("Invalid value for l. Should be greater than equal to 1.");
}

int main(int argc, char *argv[]){
  parseArgs(argc,argv);
  vector<vector<float> > testsequences,trainsequences;
  ifstream in;
  in.open(testseqfile);readSequences(testsequences,in);in.close();
  in.open(trainseqfile);readSequences(trainsequences,in);in.close();
  ofstream out(outputfile);
  for(unsigned int i = 0; i < testsequences.size(); i++){
    //extract windows of length n moving at a step of l
    if(testsequences[i].size() < (unsigned int) n){cerr<<"Line is shorter than windowsize. Exiting ... \n"; exit(0);}
    vector<float> tmpScore;
    for(unsigned int it = 0; it < testsequences[i].size()-n+1;it += l){
      vector<float> vectest(testsequences[i].begin()+it,testsequences[i].begin()+it+n);
      priority_queue<float,vector<float> > q;
      for(unsigned int j = 0; j < trainsequences.size(); j++){
	if(trainsequences[j].size() < (unsigned int) n){cerr<<"Line is shorter than windowsize. Exiting ... \n"; exit(0);}
        for(unsigned int it1 = 0; it1 < trainsequences[j].size() - n+1; it1 += l){
	  vector<float> vectrain(trainsequences[j].begin()+it1,trainsequences[j].begin()+it1+n);
	  float dist = EUC(vectest,vectrain);
	  if(q.size() < (unsigned int) nn){ q.push(dist);}
	  else{
	    if(dist < q.top()){
	      q.pop();
	      q.push(dist);
	    }
	  }
	}
      }
      tmpScore.push_back(q.top());
    }
    printSequence(tmpScore,out);
    out<<endl;
  }
  out.close();
}
