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

#include "CLUSTERDTrain.h"
#include "tools.h"

#define MAXITERS 100
#define MAXQUEUE 5

using namespace std;

char *sequence_filename=0;
char *output_filename=0;
char *simfile = 0;
int numclusters=0;
float sampling=1.0;
float eps=0.001;
int measure = 0;
int verbose = 0;
vector<vector<float> > pairwise;
#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n" \
"./CLUSTERDTrain -i sequence_filename -o output_filename -c numclusters -m measure [-s sampling-ratio -e eps -f similarity_matrix/alphamap -v verbosity]\n"\
"   -i        -> training sequence file\n"\
"   -o        -> cluster output file\n"\
"   -e [0..1] -> epsilon for convergence (default 0.001)\n"\
"   -c [2...] -> number of clusters\n"\
"   -s [0..1] -> sampling ratio (default 1.0)\n"\
"   -m        -> measure\n" \
"             ->  1 SMC\n" \
"             ->  2 WSMC\n" \
"             ->  3 LCS_DP\n" \
"             ->  4 LCS_HY\n "\
"             ->  5 SIMILARITY MATRIX\n " \
"   -f        -> similarity matrix (required when m = 5)\n" \
"             -> alphabet map (required when m = 2)\n" \
"   -v [0,1]  -> verbosity level (default 0)" \
"    \n\n"); \
  exit(1); \
}

/* This function clusters the sequences into clusters defined by the vector ids */
void findClusters(vector<vector<int> > &sequences,vector<unsigned int> &ids,vector<vector<unsigned int> > &assignments){
    for(unsigned int j = 0; j < ids.size(); j++){
	assignments.push_back(vector<unsigned int>(0));
    }
    for(unsigned int i = 0; i < sequences.size(); i++){
	//find the closest medoid to the ith sequence
	float maxsim = 0;unsigned int maxsimid = 0;
	for(unsigned int j = 0; j < ids.size(); j++){
	  float sim = pairwise[i][ids[j]];
	    if(sim > maxsim){
		maxsim = sim;
		maxsimid = j;
	    }
	}
	assignments[maxsimid].push_back(i);
    }
}

/* This function finds a medoid for a given cluster */
void findMedoids(vector<vector<int> > &sequences,vector<unsigned int> &ids,vector<vector<unsigned int> > &assignments){
  clock_t s,e;
    //do it for every cluster
    for(unsigned int i = 0; i < assignments.size(); i++){
	//for every sequence in the cluster find its average similarity to all other sequences
	//use the sampling here to make code run faster
	float maxsim = 0; unsigned int medoidId = 0;
	vector<unsigned int> v;
	unifrnd(assignments[i].size(),(int) rint(sampling*assignments[i].size()),v);
	//	unifrnd(assignments[i].size(),(int) round(sampling*assignments[i].size()),v);
	s = clock();
	for(unsigned int j = 0; j < assignments[i].size(); j++){
	    float sum = 0.0;
 	    for(unsigned int k = 0; k < v.size(); k++)
	      sum += pairwise[assignments[i][j]][assignments[i][v[k]]];
	    if(sum > maxsim){
		medoidId = assignments[i][j];
		maxsim = sum;
	    }
	}
	e = clock();
	ids.push_back(medoidId);
    }
}

/* this function evaluates the clustering as the sum of intercluster similarities */
float objFun(vector<vector<int> > &sequences,vector<unsigned int> &ids,vector<vector<unsigned int> > &assignments){
    float sum = 0;
    //do it for every cluster
    for(unsigned int i = 0; i < assignments.size(); i++){
	//for every sequence in the cluster find its similarity to the medoid
	for(unsigned int j = 0; j < assignments[i].size(); j++)
	  sum += pairwise[assignments[i][j]][ids[i]];
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

void printClusters(vector<vector<int> > &sequences,vector<unsigned int> &ids,vector<vector<unsigned int> > &assignments,ofstream &out){
    for(unsigned int i = 0; i < assignments.size(); i++){
	out<<"Cluster "<<i<<"\n";
	out<<"Medoid :";
	printSequence(sequences[ids[i]],out);
	out<<endl;
	for(unsigned int j = 0; j < assignments[i].size(); j++){
	    printSequence(sequences[assignments[i][j]],out);
	    out<<endl;
	}
    }
}

void printMedoids(vector<vector<int> > &sequences,vector<unsigned int> &ids,vector<vector<unsigned int> > &assignments,ofstream &out){
  for(unsigned int i = 0; i < assignments.size(); i++){
    printSequence(sequences[ids[i]],out);
    out<<endl;
  }
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
    }else if(strcmp(argv[argi],"-s")==0){
      if(++argi<argc) sampling=atof(argv[argi]);
      else usage("-s requires sampling ratio argument\n");
    }else if(strcmp(argv[argi],"-e")==0){
      if(++argi<argc) eps=atof(argv[argi]);
      else usage("-e requires eps argument\n");
    }else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) measure=atoi(argv[argi]);
      else usage("-m requires measure argument\n");
    }else if(strcmp(argv[argi],"-f")==0){
      if(++argi<argc) simfile=strdup(argv[argi]);
      else usage("-f requires similarity matrix argument\n");
    }else if(strcmp(argv[argi],"-v")==0){
      if(++argi<argc) verbose=atoi(argv[argi]);
      else usage("-v requires verbosity argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!sequence_filename) usage("No input sequence file specified\n");
  if(!output_filename)   usage("No output file specified\n");
  if(numclusters==0)     usage("Number of clusters not specified\n");
  if((measure < 1) || (measure > 5)) usage("Incorrect similarity measure specified.\n");
  if((measure == 2) && (!simfile)) usage("No alphabet map file specified\n");
  if((measure == 5) && (!simfile)) usage("No similarity matrix specified\n");
  if((sampling<=0.0) || (sampling > 1.0)) usage("Incorrect value for sampling ratio.\n");
}

int main(int argc, char *argv[]){
  clock_t start,end;
  parseArgs(argc,argv);
  //read sequences
  vector<vector<int> > sequences;
  //select numclusters random points
  ifstream in(sequence_filename);
  if(in.fail()) {cerr << "Error: sequence file does not exist\n"; exit(0);}
  readSequences(sequences,in);
  in.close();

  map<pair<int,int>,float> alphaMap;
  if(measure == 2){
    ifstream mapin(simfile);
    if(mapin.fail()) {cerr << "Error: alphabet map file does not exist\n"; exit(0);}
    loadMap(alphaMap,mapin);
    mapin.close();
  }
  //calculate max_symbol and max_occ_any_symbol
  int max_symbol = -1;int max_occ_any_symbol = -1;
  if(measure == 4){
    findMax(sequences,&max_symbol,&max_occ_any_symbol);
  }
  if(measure != 5){
    //find pairwise similarity between all sequences
    for(unsigned int i = 0; i < sequences.size();i++){
      vector<float> tmp(sequences.size(),0.0);
      pairwise.push_back(tmp);
    }
    for(unsigned int i = 0; i < sequences.size();i++){
      for(unsigned int j = 0; j < sequences.size();j++){
	if(j < i){
	  pairwise[i][j] = pairwise[j][i];
	}else{
	  switch(measure){
	  case 1:
	    pairwise[i][j] = SMC(sequences[i],sequences[j]);
	    break;
	  case 2:
	    pairwise[i][j] = WSMC(sequences[i],sequences[j],alphaMap);
	    break;
	  case 3:
	    pairwise[i][j] = LCS_DP(sequences[i],sequences[j]);
	    break;
	  case 4:
	    pairwise[i][j] = LCS_HY(sequences[i],sequences[j],max_symbol,max_occ_any_symbol);
	    break;
	  default:
	    cerr<<"Unsupported measure "<<measure<<"\n";
	    exit(0);
	    break;
	  }
	}
      }
    }
  }else{
    ifstream simin(simfile);
    if(simin.fail()) {cerr << "similarity matrix file does not exist\n"; exit(0);}
    readSequences(pairwise,simin);
    simin.close();
  }

  //find unique sequences
  map<vector<int>,int> unique;vector<unsigned int> uniqueids;
  for(unsigned int i = 0; i < sequences.size(); i++){
    if(unique.find(sequences[i]) == unique.end()){
      unique[sequences[i]] = i;
      uniqueids.push_back(i);
    }
  }
  if(unique.size() < (unsigned int) numclusters){cerr<<"Data has fewer unique sequences than number of clusters required.\n";exit(0);}
  //select random cluster IDs
  vector<unsigned int> ids,ids1;
  ids.resize(0);ids1.resize(0);
  unifrnd(uniqueids.size(),numclusters,ids1);
  for(unsigned int i = 0; i < ids1.size(); i++){
    ids.push_back(uniqueids[ids1[i]]);
  }
  //start clustering
  int numiters = 0;
  queue<float> objvals;
  bool stoppingFlag = false;

  while((numiters != MAXITERS) && !(stoppingFlag)){
    vector<vector<unsigned int> > assignments;
    assignments.resize(0);
    findClusters(sequences,ids,assignments);
    ids.resize(0);
    start=clock();
    findMedoids(sequences,ids,assignments);
    end = clock();
    if(verbose) cerr<<"Medoids found in "<<(end-start)/CLOCKS_PER_SEC<<" seconds\n";
    start = clock();
    float objval = objFun(sequences,ids,assignments);
    end = clock();
    if(verbose) cerr<<"Objective function computed in  "<<(end-start)/CLOCKS_PER_SEC<<" seconds\n";
    if(objvals.size() == MAXQUEUE) objvals.pop();
    objvals.push(objval);
    if(objvals.size() >= MAXQUEUE) stoppingFlag = stoppingCondition(objvals);
    numiters++;
    
    if((numiters == MAXITERS) || (stoppingFlag)){
      ofstream out(output_filename);
      printMedoids(sequences,ids,assignments,out);
      out.close();
      if(verbose) cerr<<objval<<"\n";
    }
  }
}
