#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <queue>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"

using namespace std;

char *sequence_filename=0;
char *clusters_filename=0;
char *output_filename=0;
char *simfile=0;
int measure = 0;
#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n" \
"./CLUSTERDPredict -i sequence_filename -c clusters_filename -o output_filename -m measure [-f alphabetmap]\n" \
"   -i        -> test sequence file \n" \
"   -c        -> cluster  file \n" \
"   -o        -> scores output file \n" \
"   -m        -> measure\n" \
"             ->  1 SMC\n" \
"             ->  2 WSMC\n" \
"             ->  3 LCS_DP\n" \
"             ->  4 LCS_HY\n" \
"   -f        -> alphabet map (required when m = 2)" \
"    \n\n"); \
  exit(1); \
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) sequence_filename=strdup(argv[argi]);
      else usage("-i requires input file argument\n");
    } else if(strcmp(argv[argi],"-c")==0){
      if(++argi<argc) clusters_filename=strdup(argv[argi]);
      else usage("-o requires clusters file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) output_filename=strdup(argv[argi]);
      else usage("-c requires output file argument\n");
    }else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) measure=atoi(argv[argi]);
      else usage("-m requires measure argument\n");
    }else if(strcmp(argv[argi],"-f")==0){
      if(++argi<argc) simfile=strdup(argv[argi]);
      else usage("-f requires alphabet file argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!sequence_filename) usage("No input sequence file specified\n");
  if(!clusters_filename) usage("No clustering file specified\n");
  if(!output_filename)   usage("No output file specified\n");
  if((measure < 1) || (measure > 4)) usage("Incorrect similarity measure specifie\n");
  if((measure == 2) && (!simfile)) usage("Measure 3 requires an alphabet map file\n");
}

int main(int argc, char *argv[]){
  parseArgs(argc,argv);
  //read sequences
  vector<vector<int> > sequences;
  ifstream in;

  in.open(sequence_filename);
  if(in.fail()) {cerr << "sequence file does not exist\n"; exit(0);}
  readSequences(sequences,in);
  in.close();
  //read cluster centroids
  vector<vector<int> > clusters;
  in.open(clusters_filename);
  if(in.fail()) {cerr << "clusters file does not exist\n"; exit(0);}
  readSequences(clusters,in);
  in.close();
  map<pair<int,int>,float> alphaMap;
  if(measure == 2){
    ifstream mapin(simfile);
    if(mapin.fail()) {cerr << "alphabet map file does not exist\n"; exit(0);}
    loadMap(alphaMap,mapin);
    mapin.close();
  }
  //calculate max_symbol and max_occ_any_symbol
  int max_symbol = -1;int max_occ_any_symbol = -1;
  if(measure == 4){
    findMax(sequences,&max_symbol,&max_occ_any_symbol);
  }

  ofstream out(output_filename);
  for(unsigned int i = 0; i < sequences.size(); i++){
      float maxsim = 0.0;
      for(unsigned int j = 0; j < clusters.size(); j++){
	float sim;
	switch(measure){
	case 1:
	  sim = SMC(sequences[i],clusters[j]);
	  break;
	case 2:
	  sim = WSMC(sequences[i],clusters[j],alphaMap);
	  break;
	case 3:
	  sim = LCS_DP(sequences[i],clusters[j]);
	  break;
	case 4:
	  sim = LCS_HY(sequences[i],clusters[j],max_symbol,max_occ_any_symbol);
	  break;
	default:
	  cerr<<"Unsupported measure "<<measure<<"\n";
	  exit(0);
	  break;
	}
	if(sim > maxsim) maxsim = sim;
      }
      out<<maxsim<<"\n";
  }
  out.close();
}
