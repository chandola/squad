#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <queue>
#include <math.h>
#include <time.h>

#include "tools.h"
#include "fftw3.h"

#define MAXITERS 100
#define MAXQUEUE 5

using namespace std;

char *sequence_filename=0;
char *output_filename=0;
int numclusters=0;
float eps=0.001;
int measure = 0;

/* Data structures required if measure is sim_crosscorr */
fftw_plan plan_forward,plan_backward;
double *fin, *bout;
fftw_complex *fout, *bin;

#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n\n" \
"./CLUSTERCTrain -i sequence_filename -o output_filename -c numclusters [-e eps -m measure]\n\n" \
"-i     input file name\n"\
"-o     output file for clusters\n"\
"-c     number of clusters [1 ... ]\n"\
"-e     eps for convergence of kmeans [Default 0.001]\n"\
"-m     distance measure\n"\
"       1 - Euclidean (Default)\n"\
"       2 - DTW\n"\
"       3 - DTW_SC\n"\
"       4 - DTW_SC_LB_KEOGH\n"\
"       5 - sim_crosscorr\n"\
"\n"); \
  exit(1); \
}

/* This function clusters the sequences into clusters defined by the vector ids */
void findClusters(vector<vector<float> > &sequences,vector<vector<float> > &centroids,vector<vector<unsigned int> > &assignments){
  for(unsigned int j = 0; j < centroids.size(); j++) assignments.push_back(vector<unsigned int>(0));
  for(unsigned int i = 0; i < sequences.size(); i++){
    if(measure == 5){
      //sim_crosscorr
      unsigned int maxsimid = 0;
      float maxsim = 0;
      for(unsigned int j = 0; j < centroids.size(); j++){
	float s = CROSSCORR(sequences[i],centroids[j],plan_forward,plan_backward,fin,fout,bin,bout);
	if(s > maxsim){
	  maxsim = s;
	  maxsimid = j;
	}
      }
      assignments[maxsimid].push_back(i);
    }else{
      unsigned int mindistid = 0;
      float mindist = 0;
      if(measure != 4){
	//find the closest centroid to the ith sequence
	mindist = DIST(sequences[i],centroids[0],measure);
	for(unsigned int j = 1; j < centroids.size(); j++){
	  float d = DIST(sequences[i],centroids[j],measure);
	  if(d < mindist){
	    mindist = d;
	    mindistid = j;
	  }
	}
      }else{
	mindist = queryNN(sequences[i],centroids,1,&mindistid);
      }
      assignments[mindistid].push_back(i);
    }
  }
}

/* This function finds centroid for a given cluster */
void findCentroids(vector<vector<float> > &sequences,vector<vector<float> > &centroids,vector<vector<unsigned int> > &assignments){
  centroids.resize(0);
  for(unsigned int i = 0; i < assignments.size(); i++){
    vector<float> centroid(sequences[0].size(),0.0);
    for(unsigned int k = 0; k <sequences[0].size(); k++){
      float sum = 0;
      for(unsigned int j = 0; j <assignments[i].size(); j++) sum += sequences[assignments[i][j]][k];
      sum /= assignments[i].size();
      centroid[k] = sum;
    }
    centroids.push_back(centroid);
  }
}

/* this function evaluates the clustering as the sum of intracluster distances */
float objFun(vector<vector<float> > &sequences,vector<vector<float> > &centroids,vector<vector<unsigned int> > &assignments){
    float sum = 0;
    //do it for every cluster
    for(unsigned int i = 0; i < assignments.size(); i++){
	//for every sequence in the cluster find its distance to the centroid
	for(unsigned int j = 0; j < assignments[i].size(); j++){
	  if(measure == 5){
	    sum += -1*CROSSCORR(sequences[assignments[i][j]],centroids[i],plan_forward,plan_backward,fin,fout,bin,bout);
	  }else{
	    int measure1 = measure;	  
	    if(measure == 4) measure1 = 3;
	    sum += DIST(sequences[assignments[i][j]],centroids[i],measure1);
	  }
	}
    }
    return sum;
}

bool stoppingCondition(queue<float> &objvals){
    vector<float> vals;float minval=objvals.front();float maxval=objvals.front();
    while(!objvals.empty()){
	float val = objvals.front();
	objvals.pop();vals.push_back(val);
	if(val < minval) minval = val;
	if(val > maxval) maxval = val;
    }
    for(unsigned int i = 0; i < vals.size(); i++) objvals.push(vals[i]);
    if(((maxval - minval)/maxval) < eps) return true;
    else return false;
}

void unifrnd(int N, unsigned int k, vector<unsigned int> &s){
  // generated k random numbers in the interval [0,N-1]
  // if N > k then generates unique random numbers only
  unsigned long r;
  unsigned int i;
  map<int,int> ids_map;
  srand( (unsigned)time( NULL ) );
  while(s.size()!=k){
    r = rand(); r <<= 15; r += rand();
    i = r%N;
    if(N > (int) k){
      if(ids_map.find(i) == ids_map.end()){
	s.push_back(i);
	ids_map[i] = i;
      }
    }else{
      s.push_back(i);
    }
  }
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) sequence_filename=strdup(argv[argi]);
      else usage("-i requires input file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) output_filename=strdup(argv[argi]);
      else usage("-o requires output file argument\n");
    } else if(strcmp(argv[argi],"-c")==0){
      if(++argi<argc) numclusters=atoi(argv[argi]);
      else usage("-c requires number of clusters argument\n");
    }else if(strcmp(argv[argi],"-e")==0){
      if(++argi<argc) eps=atof(argv[argi]);
      else usage("-e requires eps argument\n");
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
  if(!output_filename)   usage("No output file specified\n");
  if(numclusters==0)     usage("Number of clusters not specified\n");
  if((measure < 1) || (measure > 5) ) usage("Unsupported Measure\n");
}

int main(int argc, char *argv[]){
  parseArgs(argc,argv);
  //read sequences
  vector <vector <float> > sequences;
  ifstream in(sequence_filename);
  if(in.fail()) {cerr << "sequence file does not exist\n"; exit(0);}
  readSequences(sequences,in);
  in.close();

  if(sequences.size() < (unsigned int) numclusters){cerr<<"Data has fewer sequences than number of clusters required.\n";exit(0);}
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
  //select random cluster IDs
  vector<unsigned int> ids;
  ids.resize(0);
  unifrnd(sequences.size(),numclusters,ids);
  vector<vector<float> > centroids(0);
  for(unsigned int i = 0; i < (unsigned int) numclusters;i++) centroids.push_back(sequences[ids[i]]);

  //start clustering
  int numiters = 0;
  queue<float> objvals;
  bool stoppingFlag = false;
  clock_t start, end;
  start = clock();
  while((numiters != MAXITERS) && !(stoppingFlag)){
    vector<vector<unsigned int> > assignments;
    assignments.clear();assignments.resize(0);
    findClusters(sequences,centroids,assignments);

    centroids.clear();centroids.resize(0);
    findCentroids(sequences,centroids,assignments);
    float objval = objFun(sequences,centroids,assignments);
    if(objvals.size() == MAXQUEUE) objvals.pop();
    objvals.push(objval);
    if(objvals.size() >= MAXQUEUE) stoppingFlag = stoppingCondition(objvals);
    numiters++;
    if((numiters == MAXITERS) || (stoppingFlag)){
      ofstream out(output_filename);
      printSequences(centroids,out);
      out.close();
      cerr<<numclusters<<"\t"<<objval<<"\t";
      break;
    }
  }
  end = clock();
  cerr<<(float) (end-start)/CLOCKS_PER_SEC<<"\t"<<numiters<<"\n";
  if(measure == 5){
    fftw_destroy_plan(plan_forward);    fftw_destroy_plan(plan_backward);
    fftw_free(fout); fftw_free(bin);free(fin);free(bout);
  }
}
