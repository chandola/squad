#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <queue>
#include <math.h>

#include "tools.h"
#include "fftw3.h"


using namespace std;

char *testsequence_filename=0;
char *trainsequence_filename=0;
char *output_filename=0;

int measure = 0;
/* Data structures required if measure is sim_corrcoef */
fftw_plan plan_forward,plan_backward;
double *fin, *bout;
fftw_complex *fout, *bin;
#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n" \
"./findsimc -i test_sequence_filename -t train_sequence_filename -o output_filename [-m measure]\n\n" \
"   -i        -> test sequence file \n" \
"   -t        -> train sequence file \n" \
"   -o         -> similarity output file \n" \
"   -m         -> measure\n" \
"              ->  1 DIST_EUC  (default)\n" \
"              ->  2 DIST_DTW\n" \
"              ->  3 DIST_DTW_SC\n" \
"              ->  4 DIST_DTW_LB\n" \
"              ->  5 SIM_CROSSCORR\n"); \
  exit(1); \
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-t")==0){
      if(++argi<argc) trainsequence_filename=strdup(argv[argi]);
      else usage("-t requires traininput file argument\n");
    }else  if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) testsequence_filename=strdup(argv[argi]);
      else usage("-i requires testinput file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) output_filename=strdup(argv[argi]);
      else usage("-o requires output file argument\n");
    }else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) measure=atoi(argv[argi]);
      else usage("-m requires measure argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!testsequence_filename) usage("No test sequence file specified\n");
  if(!trainsequence_filename) usage("No train sequence file specified\n");
  if(!output_filename)   usage("No output file specified\n");
  if((measure < 1) || (measure > 5) ) usage("Unsupported Measure\n");
}

int main(int argc, char *argv[]){
  parseArgs(argc,argv);
  //read sequences
  vector<vector<float> > testsequences;
  vector<vector<float> > trainsequences;
  ifstream in(testsequence_filename);
  if(in.fail()) {cerr << "test sequence file does not exist\n"; exit(0);}
  readSequences(testsequences,in);
  in.close();
  in.open(trainsequence_filename);
  if(in.fail()) {cerr << "train sequence file does not exist\n"; exit(0);}
  readSequences(trainsequences,in);
  in.close();

  /* Prepare FFT plans and allocate data structures if measure is sim_crosscorr*/
  int len = (int) trainsequences[0].size();
  if(measure == 5){
    fin  = (double *) malloc(sizeof(double) * len);
    fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
    bin  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
    bout = (double *) malloc(sizeof(double) * len);
    plan_forward  = fftw_plan_dft_r2c_1d(len, fin, fout, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_1d(len, bin, bout, FFTW_ESTIMATE);
  }

  //find pairwise distance between all sequences
  vector<vector<float> > pairwise;
  if(measure != 5){
    pairwiseDist(testsequences,trainsequences,pairwise,measure);
  }else{
    for(unsigned int i = 0; i < testsequences.size(); i++){
      vector<float> tmp(trainsequences.size(),0.0); pairwise.push_back(tmp);
      for(unsigned int j = 0; j < trainsequences.size(); j++){
	pairwise[i][j] = -1*CROSSCORR(testsequences[i],trainsequences[j],plan_forward,plan_backward,fin,fout,bin,bout);
      }
    }
  }
  //print out pairwise distance
  if(measure == 5){
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(fout); 
    fftw_free(bin);
    free(fin);free(bout);
  }
  ofstream out(output_filename);
  printSequences(pairwise,out);
  out.close();
}
