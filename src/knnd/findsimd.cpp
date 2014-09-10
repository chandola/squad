#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <queue>
#include <math.h>
#include <time.h>

#include "tools.h"

using namespace std;

char *test_sequence_filename=0;
char *train_sequence_filename=0;
char *output_filename=0;
char *alphamapfile=0;
int measure = 0;
int level = 0;
#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n" \
"./findsimd -i test_sequence_filename -t train_sequence_filename -o output_filename -m measure [-a alphamapfile -l level]\n\n" \
"   -i        -> test sequence file \n" \
"   -t        -> train sequence file \n" \
"   -o        -> similarity output file \n" \
"   -m        -> measure\n" \
"             ->  1 SMC\n" \
"             ->  2 WSMC\n" \
"             ->  3 LCS_DP\n" \
"             ->  4 LCS_HY\n" \
"             ->  5 BITMAP\n" \
"   -a        -> alphabet map (required when m = 2)" \
"   -l        -> level (required when m = 5)" \
"    \n\n"); \
  exit(1); \
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) test_sequence_filename=strdup(argv[argi]);
      else usage("-i requires input file argument\n");
    }else  if(strcmp(argv[argi],"-t")==0){
      if(++argi<argc) train_sequence_filename=strdup(argv[argi]);
      else usage("-t requires input file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) output_filename=strdup(argv[argi]);
      else usage("-o requires output file argument\n");
    }else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) measure=atoi(argv[argi]);
      else usage("-m requires measure argument\n");
    }else if(strcmp(argv[argi],"-a")==0){
      if(++argi<argc) alphamapfile=strdup(argv[argi]);
      else usage("-a requires alphamap argument\n");
    }else if(strcmp(argv[argi],"-l")==0){
      if(++argi<argc) level=atoi(argv[argi]);
      else usage("-l requires level argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!test_sequence_filename) usage("No test sequence file specified\n");
  if(!train_sequence_filename) usage("No train sequence file specified\n");
  if(!output_filename)   usage("No output file specified\n");
  if((measure == 2) &&  !(alphamapfile)) usage("Measure 2 (WSMC) requires alphamapfile\n");
  if((measure == 5) &&  (level==0)) usage("Measure 5 (BITMAP) requires level value\n");
}

int main(int argc, char *argv[]){
  parseArgs(argc,argv);
  //read sequences
  vector<vector<int> > trainSequences;
  vector<vector<int> > testSequences;
  //select numclusters random points
  ifstream in(train_sequence_filename);
  if(in.fail()) {cerr << "train sequence file does not exist\n"; exit(0);}
  readSequences(trainSequences,in);
  in.close();in.clear();
  in.open(test_sequence_filename);
  if(in.fail()) {cerr << "test sequence file does not exist\n"; exit(0);}
  readSequences(testSequences,in);
  in.close();
  //find pairwise similarity between all sequences
  vector<vector<float> > pairwise;
  for(unsigned int i = 0; i < testSequences.size();i++){
    pairwise.push_back(vector<float>(trainSequences.size(),0.0));
  }
  map<pair<int,int>,float> alphaMap;
  if(measure == 2){
    ifstream mapin(alphamapfile);
    if(mapin.fail()) {cerr << "Error: alphabet map file does not exist\n"; exit(0);}
    loadMap(alphaMap,mapin);
    mapin.close();
  }
  //calculate max_symbol and max_occ_any_symbol
  int max_symbol = -1;int max_occ_any_symbol = -1;
  if(measure == 4){
    findMax(testSequences,&max_symbol,&max_occ_any_symbol);
    findMax(trainSequences,&max_symbol,&max_occ_any_symbol);
  }
  for(unsigned int i = 0; i < testSequences.size();i++){
    for(unsigned int j = 0; j < trainSequences.size();j++){
      switch(measure){
      case 1:
	pairwise[i][j] = SMC(testSequences[i],trainSequences[j]);
	break;
      case 2:
	pairwise[i][j] = WSMC(testSequences[i],trainSequences[j],alphaMap);
	break;
      case 3:
	pairwise[i][j] = LCS_DP(testSequences[i],trainSequences[j]);
	break;
      case 4:
	pairwise[i][j] = LCS_HY(testSequences[i],trainSequences[j],max_symbol,max_occ_any_symbol);
	break;
      case 5:
	pairwise[i][j] = -1*BITMAP(testSequences[i],trainSequences[j],level);
	break;
      default:
	cerr<<"Unsupported measure "<<measure<<"\n";
	exit(0);
	break;
      }
    }
  }
 //print out pairwise similarity
  ofstream out(output_filename);
  printSequences(pairwise,out);
  out.close();
}
