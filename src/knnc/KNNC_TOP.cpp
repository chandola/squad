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
char *train_filename = 0;
char *test_filename = 0;
char *output_filename = 0;
int nn=0;
unsigned int x=0;
vector<vector<float> > pairwise;
#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n\n" \
"./knn -i test_filename -t train_filename -o output_filename -k nearest_neighbors [-x top_anomalies]\n\n"\
"-i     test file name\n"\
"-t     train file name\n"\
"-o     output file for anomaly scores\n"\
"-k     number of nearest neighbors\n"\
"-x     top anomalies to be scored (default all)\n"\
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
    } else if(strcmp(argv[argi],"-x")==0){
      if(++argi<argc) x=atoi(argv[argi]);
      else usage("-x requires number of top anomalies argument\n");
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
}
class compare_A{
public:
  bool operator() (pair<float,unsigned int> a, pair<float,unsigned int>  b){
    return a.first > b.first;
  }
};

int main(int argc, char *argv[]){
  parseArgs(argc,argv);
  ifstream trainin(train_filename);
  ifstream testin(test_filename);
  if(trainin.fail()) {cerr << "training file does not exist\n"; exit(0);}
  if(testin.fail()) {cerr << "testing file does not exist\n"; exit(0);}
  vector<vector<float> > trainSequences;
  vector<vector<float> > testSequences;
  readSequences(trainSequences,trainin);
  readSequences(testSequences,testin);
  trainin.close();testin.close();
  if(x == 0){ x = testSequences.size();}
  priority_queue<pair<float,unsigned int>,vector<pair<float,unsigned int> >, compare_A > q;
  for(unsigned int i = 0; i < testSequences.size(); i++){
    unsigned int index = 0;
      
    if(q.size() < x){
      float d = queryNN(testSequences[i],trainSequences,nn,&index);
      q.push(pair<float,unsigned int>(d,i));
    }else{
      float d = queryNN(testSequences[i],trainSequences,nn,&index);
      if(d > q.top().first){
 	q.pop();
	q.push(pair<float,unsigned int>(d,i));
      }
    }
  }
  map<unsigned int,float> m;
  while(!q.empty()){
    m[q.top().second] = q.top().first;
    q.pop();
  }  
  vector<float> scores(testSequences.size(),0.0);
  for(map<unsigned int,float>::iterator it = m.begin(); it != m.end(); it++)
    scores[it->first] = it->second;

  ofstream out;out.open(output_filename);
  for(unsigned int i = 0; i < testSequences.size(); i++)
    out<<scores[i]<<"\n";
  out.close();
}
