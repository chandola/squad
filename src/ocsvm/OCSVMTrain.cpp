#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "svm.h"
#include "sax.h"
#include "tools.h"

using namespace std;
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

char *infilename;
char *modelfilename;
int len = 0;
double nu = 0.1;
double dgamma = -1.0;
int hop = 1;
int segment = 1;

struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;

#define usage(msg, args...) \
{\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n\n" \
"./SVRTrain -i sequence_filename -m model_filename -l length [-g gamma -n nu -s segment -x hop]\n\n" \
"   -n [0..1]  -> nu parameter for one class SVM (default 0.5)\n"\
"   -l [2..]  -> window length\n" \
"   -g        -> gamma parameter for rbf kernel\n"\
"   -x        -> number of hops (default 1)\n"\
"   -s        -> segment length (default 1)\n"\
"\n");\
  exit(1); \
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) infilename=strdup(argv[argi]);
      else usage("-i requires input file argument\n");
    } else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) modelfilename=strdup(argv[argi]);
      else usage("-m requires model file argument\n");
    }else if(strcmp(argv[argi],"-l")==0){
      if(++argi<argc) len=atoi(argv[argi]);
      else usage("-l requires measure argument\n");
    }else if(strcmp(argv[argi],"-n")==0){
      if(++argi<argc) nu=atof(argv[argi]);
      else usage("-n requires nu argument\n");
    }else if(strcmp(argv[argi],"-g")==0){
      if(++argi<argc) dgamma=atof(argv[argi]);
      else usage("-g requires gamma parameter argument\n");
    }else if(strcmp(argv[argi],"-x")==0){
      if(++argi<argc) hop=atoi(argv[argi]);
      else usage("-x requires hop argument\n");
    }else if(strcmp(argv[argi],"-s")==0){
      if(++argi<argc) segment=atoi(argv[argi]);
      else usage("-s requires segment argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!infilename) usage("No input sequence filename specified\n");
  if(!modelfilename)   usage("No model filename specified\n");
  if(len==0)     usage("Window length not specified\n");
  if(dgamma == -1.0) dgamma = (float)1/len;
}

int main (int argc, char *argv[]){ 
  parseArgs(argc,argv);
  
  //prepare the parameters
  param.svm_type = ONE_CLASS;
  param.degree = 1;
  param.coef0 = 0;
  param.cache_size = 100;
  param.C = 1.0;
  param.eps = 1e-3;
  param.p = param.eps;
  param.shrinking = 1;
  param.probability = 0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;

  param.kernel_type = RBF;
  param.gamma = dgamma;
  param.nu = nu;
  // Read training set
  vector<vector<float> > ts,ts_segs,windows;
  ifstream in(infilename);
  readSequences(ts,in);
  in.close();
  if(segment > 1) TStoPAA(ts,ts_segs,segment);
  else ts_segs = ts;
  vector<unsigned int> lens(0);
  getWindows_hop(ts_segs,windows,len,hop,lens);
  if(windows.size() == 0){
    cerr<<"Error: No windows created \n";
    exit(0);
  }
  ts.clear();ts_segs.clear();
  int elements = 0;
  for(unsigned int i = 0; i < lens.size(); i++)
    elements += lens[i];
  elements *= len;
  //check parameters
  const char *error_msg = svm_check_parameter(&prob,&param);
  if(error_msg)	{
    fprintf(stderr,"Error: %s\n",error_msg);
    exit(1);
  }
  //prepare data in svm format
  prob.l = (int) windows.size();
  prob.y = Malloc(double,prob.l);
  prob.x = Malloc(struct svm_node *,prob.l);
  x_space = Malloc(struct svm_node,elements+windows.size());
  unsigned int max_index = windows[0].size();
  int j=0;
  for(unsigned int i=0;i<windows.size();i++){
    prob.x[i] = &x_space[j];
    //fill data in x_space and prob.y[i]
    for(int k = 0; k < (int) max_index; k++){
      x_space[j].index =  k+1;
      x_space[j].value = (double) windows[i][k];
      j++;
    }
    prob.y[i] = 1.00;
    x_space[j++].index = -1;
  }
  windows.clear();
  //run svm-learn
  model = svm_train(&prob,&param);

  svm_save_model(modelfilename,model);
  svm_destroy_model(model);
  free(prob.y);
  free(prob.x);
  free(x_space);
}
