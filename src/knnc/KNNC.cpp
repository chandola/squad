#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <queue>
#include <math.h>
#include <time.h>

#include "tools.h"
#include "fftw3.h"

/* Data structures required if measure is sim_crosscorr */
fftw_plan plan_forward,plan_backward;
double *fin, *bout;
fftw_complex *fout, *bin;

using namespace std;
char *train_filename = 0;
char *test_filename = 0;
char *output_filename = 0;
int nn=0;
int measure=0;
vector<vector<float> > pairwise;
#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n\n" \
"./KNNC -i test_filename -t train_filename -o output_filename -k nearest_neighbors [-m measure]\n\n"\
"-i     test file name\n"\
"-t     train file name\n"\
"-o     output file for anomaly scores\n"\
"-k     number of nearest neighbors\n"\
"-m     distance measure\n"\
"       1 - Euclidean (Default)\n"\
"       2 - DTW\n"\
"       3 - DTW_SC\n"\
"       4 - DTW_SC_LB_KEOGH\n"\
"       5 - sim_crosscorr\n"\
); \
  exit(1); \
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-t")==0){
      if(++argi<argc) train_filename=strdup(argv[argi]);
      else usage("-t requires train file argument\n");
    } else if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) test_filename=strdup(argv[argi]);
      else usage("-i requires test file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) output_filename=strdup(argv[argi]);
      else usage("-o requires output file argument\n");
    } else if(strcmp(argv[argi],"-k")==0){
      if(++argi<argc) nn=atoi(argv[argi]);
      else usage("-k requires number of nearest neighbors argument\n");
    }else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) measure=atoi(argv[argi]);
      else usage("-m requires measure argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!train_filename) usage("No training file specified.\n");
  if(!test_filename) usage("No testing file specified\n");
  if(!output_filename)   usage("No output file specified\n");
  if(nn == 0) usage("Nearest Neighbors not specified\n");
  if((measure < 1) || (measure > 5) ) usage("Unsupported Measure\n");
}
bool compare_dist(float a, float b){return a > b;}

int main(int argc, char *argv[]){
  parseArgs(argc,argv);
  ifstream trainin(train_filename);
  ifstream testin(test_filename);
  if(trainin.fail()) {cerr << "training file does not exist\n"; exit(0);}
  if(testin.fail()) {cerr << "testing file does not exist\n"; exit(0);}
  vector<vector<float> > trainTS;
  vector<vector<float> > testTS;
  readSequences(trainTS,trainin);
  readSequences(testTS,testin);
  trainin.close();testin.close();
  
  /* Prepare FFT plans and allocate data structures if measure is sim_crosscorr*/
  int len = (int) trainTS[0].size();
  if(measure == 5){
    fin  = (double *) malloc(sizeof(double) * len);
    fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
    bin  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
    bout = (double *) malloc(sizeof(double) * len);
    plan_forward  = fftw_plan_dft_r2c_1d(len, fin, fout, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_1d(len, bin, bout, FFTW_ESTIMATE);
  }

  ofstream out;out.open(output_filename);
  for(unsigned int i = 0; i < testTS.size(); i++){
    if((measure != 2) && (measure != 3)){
      if(testTS[i].size() != (unsigned int) len){
	fprintf(stderr,"ERROR : Test and training Sequences should be of equal length for this measure. Exiting ... \n");
	break;
      }
    }
    if(measure != 4){
      priority_queue<float,vector<float> > q;
      for(unsigned int j = 0; j < trainTS.size(); j++){
	float d;
	if(measure != 5){
	  d = DIST(testTS[i],trainTS[j],measure);
	}else {
	  d = -1*CROSSCORR(testTS[i],trainTS[j],plan_forward,plan_backward,fin,fout,bin,bout); 
	}
	if(q.size() < (unsigned int) nn){
	  q.push(d);
	}else{
	  if(d < q.top()){
	    q.pop();
	    q.push(d);
	  }
	}
      }
      //find the k^th nearest neighbor
      out<<q.top()<<"\n";
    }else{
      unsigned int index = 0;
      float d = queryNN(testTS[i],trainTS,nn,&index);
      out<<d<<"\n";
    }
  }
  out.close();
  if(measure == 5){
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(fout); 
    fftw_free(bin);
    free(fin);free(bout);
  }
}
