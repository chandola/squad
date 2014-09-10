//FSATrain.cpp - Varun Chandola

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "FSA.h"
using namespace std;

char *seqfile = 0;//file containing train sequences
int max_alphabet=INT_MAX;//maximum alphabet
char *modelfile = 0;//file to which model will be written
char *oldmodelfile = 0;//file containing old model
int n=0;//length of the history
int l=1;//size of current state (default set to 1)

#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n\n" \
"./FSATrain -i training_sequences -o fsa_model_filename -n history_length [-m old_modelfile -a max_alphabet -l state_length]\n\n"); \
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
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) modelfile=strdup(argv[argi]);
      else usage("-o requires model file argument\n");
    } else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) oldmodelfile=strdup(argv[argi]);
      else usage("-m requires old model file argument\n");
    } else if(strcmp(argv[argi],"-n")==0){
      if(++argi<argc) n=atoi(argv[argi]);
      else usage("-n requires sequence length argument\n");
    }else if(strcmp(argv[argi],"-l")==0){
      if(++argi<argc) l=atoi(argv[argi]);
      else usage("-l requires state length size argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!seqfile)     usage("No test file specified\n");
  if(!modelfile)   usage("No output file specified\n");
  
  if(!oldmodelfile){
    if((n <= 0) || (n > 1000)) usage("Invalid value for n. Should be between 1 and 1000");
    if((l <= 0) || (l > 1000)) usage("Invalid value for l. Should be between 1 and 1000");
  }
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
    //load model if it exists
    if(oldmodelfile){
      in.open(oldmodelfile);
      loadFSA(nodesMap,in,str2intnodemap,int2strnodemap,str2intedgemap,int2stredgemap,&n,&l,&totaledges);
      in.close();in.clear();
    }
    //read sequences and update the model
    vector<vector<vector<int> > > windows;
    readSequences(windows,max_alphabet,seqfile,n+l);
    trainSequences(nodesMap,str2intnodemap,int2strnodemap,str2intedgemap,int2stredgemap,windows,n,l,&totaledges);

    out.open(modelfile);
    printFSA(nodesMap,out,str2intnodemap,int2strnodemap,str2intedgemap,int2stredgemap,n,l,totaledges);
    out.close();out.clear();
}
