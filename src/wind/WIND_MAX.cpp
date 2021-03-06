//discord_max.cpp - Varun Chandola

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
char *alphamapfile = 0;//file containing alphabet map
char *outputfile = 0;//file to which results will be written
int n = 2;//size of the window
int l = 1;//sliding step
int nn = 10;//nearest neighbors
int x = 1;//top discords within each test sequence to be reported.
map<pair<int,int>,float> alphaMap;

#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n\n" \
"./DISCORDD_MAX -i test_sequence_filename -t train_sequence_filename -o outputfile [-m alphamap_filename -n window -l hop -nn nearestneighbors -x max]\n\n"); \
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
    } else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) alphamapfile=strdup(argv[argi]);
      else usage("-m requires alphabet map file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) outputfile=strdup(argv[argi]);
      else usage("-o requires output file argument\n");
    } else if(strcmp(argv[argi],"-n")==0){
      if(++argi<argc) n=atoi(argv[argi]);
      else usage("-n requires sequence length argument\n");
    }else if(strcmp(argv[argi],"-l")==0){
      if(++argi<argc) l=atoi(argv[argi]);
      else usage("-l requires hop length size argument\n");
    }else if(strcmp(argv[argi],"-x")==0){
      if(++argi<argc) x=atoi(argv[argi]);
      else usage("-x requires number of top discords argument\n");
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
  if(n <= 0 ) usage("Invalid value for n. Should be greater than 0\n");
  if(l <= 0 ) usage("Invalid value for l. Should be greater than 0\n");
  if(l > n) usage("Window length should be greater than equal to hop length\n");
}


int main(int argc, char *argv[]){
  parseArgs(argc,argv);
  if(alphamapfile){
    ifstream amfilein; amfilein.open(alphamapfile);
    if(amfilein.fail()) {cerr << "alphabet map file does not exist\n"; exit(0);}
    loadMap(alphaMap,amfilein);
    amfilein.close();
  }
  vector<vector<int> > testsequences,trainsequences;
  ifstream in;
  in.open(testseqfile);readSequences(testsequences,in);in.close();
  in.open(trainseqfile);readSequences(trainsequences,in);in.close();
  ofstream out(outputfile);
  for(unsigned int i = 0; i < testsequences.size(); i++){
    //extract windows of length n moving at a step of l
    if(testsequences[i].size() < (unsigned int) n){cerr<<"Line is shorter than windowsize. Exiting ... \n"; exit(0);}
    priority_queue<float,vector<float> > qq;
    for(unsigned int it = 0; it < testsequences[i].size()-n+1;it += l){
      vector<int> vectest(testsequences[i].begin()+it,testsequences[i].begin()+it+n);
      bool flag = false;
      priority_queue<float,vector<float>,greater<float> > q;
      for(unsigned int j = 0; j < trainsequences.size(); j++){
	if(trainsequences[j].size() < (unsigned int) n){cerr<<"Line is shorter than windowsize. Exiting ... \n"; exit(0);}
        for(unsigned int it1 = 0; it1 < trainsequences[j].size() - n+1; it1 += l){
	  vector<int> vectrain(trainsequences[j].begin()+it1,trainsequences[j].begin()+it1+n);
	  float sim = 0;
	  if(alphamapfile) sim = WSMC(vectest, vectrain, alphaMap);
	  else sim = SMC(vectest,vectrain);
	  if(q.size() < (unsigned int) nn){ 
	    q.push(sim);
	  }else{
	    if(sim > q.top()){
	      q.pop();
	      q.push(sim);
	    }
	    //prune here
	    if(qq.size() >= (unsigned int) x){
	      if(q.top() > qq.top()){
		flag = true;
		break;
	      }
	    }
	  }
	}
	if(flag) break;
      }
      if(flag) continue;
      if(qq.size() < (unsigned int) x){
	qq.push(q.top());
      }else{
	if(q.top() < qq.top()){
	  qq.pop();
	  qq.push(q.top());
	}
      }
    }
    while(!qq.empty()){
      out<<qq.top()<<" ";
      qq.pop();
    }
    out<<"\b\n";
  }
  out.close();
}
