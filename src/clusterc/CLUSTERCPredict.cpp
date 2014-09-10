#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <queue>
#include <math.h>

#include "tools.h"
#include "fftw3.h"

using namespace std;

char *sequence_filename=0;
char *clusters_filename=0;
char *output_filename=0;
int measure = 0;

/* Data structures required if measure is sim_crossorr */
fftw_plan plan_forward,plan_backward;
double *fin, *bout;
fftw_complex *fout, *bin;

#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n\n" \
"./CLUSTERCPredict -i sequence_filename -c clustersfilename -o output_filename [-m measure]\n\n" \
"-i     input file name\n"\
"-o     output file for anomaly scores\n"\
"-c     clusters file\n"\
"-m     distance measure\n"\
"       1 - Euclidean (Default)\n"\
"       2 - DTW\n"\
"       3 - DTW_SC\n"\
"       4 - DTW_SC_LB_KEOGH\n"\
"       5 - sim_crossorr\n"\
"\n"); \
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
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!sequence_filename) usage("No input sequence file specified\n");
  if(!clusters_filename) usage("No clustering file specified\n");
  if(!output_filename)   usage("No output file specified\n");
  if((measure < 1) || (measure > 5)) usage("Unsupported Measure\n");
}

int main(int argc, char *argv[]){
  parseArgs(argc,argv);
  //read sequences
  vector<vector<float> > sequences;
  ifstream in;

  in.open(sequence_filename);
  if(in.fail()) {cerr << "sequence file does not exist\n"; exit(0);}
  readSequences(sequences,in);
  in.close();
  //read cluster centroids
  vector<vector<float> > clusters;
  in.open(clusters_filename);
  if(in.fail()) {cerr << "clusters file does not exist\n"; exit(0);}
  readSequences(clusters,in);
  in.close();

  /* Prepare FFT plans and allocate data structures if measure is sim_crosscorr*/
  int len = (int) sequences[0].size();
  if(measure == 5){
    fin  = (double *) malloc(sizeof(double) * len);
    fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
    bin  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
    bout = (double *) malloc(sizeof(double) * len);
    plan_forward  = fftw_plan_dft_r2c_1d(len, fin, fout, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_1d(len, bin, bout, FFTW_ESTIMATE);
  }

  clock_t start,end;
  ofstream out(output_filename);
  start = clock();

  for(unsigned int i = 0; i < sequences.size(); i++){
    float mindist=0.0;
    if(measure == 5){
      float maxsim = 0;
      for(unsigned int j = 0; j < clusters.size(); j++){
	float s = CROSSCORR(sequences[i],clusters[j],plan_forward,plan_backward,fin,fout,bin,bout);
	if(s > maxsim) maxsim = s;	
      }
      mindist = -1*maxsim;
    }else{
      if(measure != 4){
	mindist = DIST(sequences[i],clusters[0],measure);
	for(unsigned int j = 1; j < clusters.size(); j++){
	  float d = DIST(sequences[i],clusters[j],measure);
	  if(d < mindist) mindist = d;
	}
      }else{
	unsigned int id = 0;
	mindist = queryNN(sequences[i],clusters,1,&id);	
      }
    }
    out<<mindist<<"\n";
  }
  end = clock();
  cerr<<(float) (end-start)/CLOCKS_PER_SEC<<"\n";
  out.close();
  if(measure == 5){
    fftw_destroy_plan(plan_forward);    fftw_destroy_plan(plan_backward);
    fftw_free(fout); fftw_free(bin);free(fin);free(bout);
  }
}
