//FSAPredict.cpp - Varun Chandola

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "FSA.h"

using namespace std;

char *seqfile = 0;//file containing input test sequences
int max_alphabet=INT_MAX;//maximum alphabet
char *modelfile = 0;//file containing FSA model
char *outputfile = 0;//file to which results will be written
int n=0;//length of the history
int l=1;//size of current state (default set to 1)

#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n\n" \
"FSAPredict -i sequence_filename -m model_filename -o outputfile -a max_alphabet\n\n"); \
  exit(1); \
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) seqfile=strdup(argv[argi]);
      else usage("-i requires input file argument\n");
    } else if(strcmp(argv[argi],"-a")==0){
      if(++argi<argc) max_alphabet=atoi(argv[argi]);
      else usage("-a requires max alphabet argument\n");
    } else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) modelfile=strdup(argv[argi]);
      else usage("-m requires model file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) outputfile=strdup(argv[argi]);
      else usage("-o requires output file argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!seqfile)     usage("No test file specified\n");
  if(!modelfile)   usage("No model file specified\n");
  if(!outputfile)   usage("No output file specified\n");
}


int main(int argc, char *argv[]){
    parseArgs(argc,argv);    
    //create a map of nodes and edges
    map<int,node> nodesMap;
    int totaledges = 0;

    map<vector<int>,int> str2intnodemap;
    map<int,vector<int> > int2strnodemap;
    map<vector<int>,int> str2intedgemap;
    map<int,vector<int> > int2stredgemap;
    ofstream out;ifstream in;
    in.open(modelfile);
    if(in.fail()) {cerr << "model file does not exist\n"; exit(0);}
    //load model
    loadFSA(nodesMap,in,str2intnodemap,int2strnodemap,str2intedgemap,int2stredgemap,&n,&l,&totaledges);
    in.close();in.clear();
    //predict sequences
    out.open(outputfile);
    vector<vector<vector<int> > > windows;
    readSequences(windows,max_alphabet,seqfile,n+l);
    predictSequences(nodesMap,str2intnodemap,int2strnodemap,str2intedgemap,int2stredgemap,windows,out,n,l,totaledges,-1);
    //sequences predicted
    in.close();in.clear();    out.close();out.clear();
}
