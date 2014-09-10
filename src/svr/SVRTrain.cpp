#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "tools.h"
#include "svm.h"
#include "sax.h"

using namespace std;
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

char *infilename;
char *modelfilename;
int len = 0;
double eps = 0.1;
int kernel = 0;
int degree = -1;
double dgamma = -1.0;
int hop = 1;
double C = -1.0;
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
"./SVRTrain -i sequence_filename -m model_filename -l length [-e epsilon -k kernel -d degree -g gamma -s segment -x hop -C C]\n\n" \
"   -e [0..]  -> epsilon width of tube for regression (default 0.1)\n"\
"   -l [2..]  -> window length\n" \
"   -k [0..2] -> kernel 0-linear (default) 1-polynomial 2-rbf\n"\
"   -d        -> degree parameter (required for polynomial kernels)\n"\
"   -g        -> gamma parameter (required for polynomial and rbf kernels)\n"\
"   -x        -> number of hops (default 1)\n"\
"   -s        -> segment length (default 1)\n"\
"   -C        -> trade off between training error and margin (default avg((x*x)^-1)) "\
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
    }else if(strcmp(argv[argi],"-e")==0){
      if(++argi<argc) eps=atof(argv[argi]);
      else usage("-e requires epsilon argument\n");
    }else if(strcmp(argv[argi],"-k")==0){
      if(++argi<argc) kernel=atoi(argv[argi]);
      else usage("-k requires kernel argument\n");
    }else if(strcmp(argv[argi],"-d")==0){
      if(++argi<argc) degree=atoi(argv[argi]);
      else usage("-d requires degree parameter argument\n");
    }else if(strcmp(argv[argi],"-g")==0){
      if(++argi<argc) dgamma=atof(argv[argi]);
      else usage("-g requires gamma parameter argument\n");
    }else if(strcmp(argv[argi],"-C")==0){
      if(++argi<argc) C=atof(argv[argi]);
      else usage("-C requires C parameter argument\n");
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
  if((kernel < 0) || (kernel > 2)) usage("Kernel argument should be 0,1, or 2.\n");
  if(kernel == 1){ //POLY Kernel
    if((degree <= 0) || (dgamma == -1.0)){
      usage("POLY kernel requires degree and gamma parameter.\n");
    }
  }
  if(kernel == 2){ //RBF Kernel
    if(dgamma == -1.0){
      usage("RBF kernel requires gamma parameter.\n");
    }
  }
}

int main (int argc, char *argv[]){ 
  parseArgs(argc,argv);
  
  //prepare the parameters
  param.svm_type = EPSILON_SVR;
  param.degree = 3;
  param.gamma = 0; 
  param.coef0 = 0;
  param.nu = 0.5;
  param.cache_size = 100;
  param.C = C;
  param.eps = 1e-3;
  param.p = eps;
  param.shrinking = 1;
  param.probability = 0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;
  switch(kernel){
  case 0:
    param.kernel_type = LINEAR;
    break;
  case 1:
    param.kernel_type = POLY;
    param.degree = degree;
    param.gamma = dgamma;
    break;
  case 2:
    param.kernel_type = RBF;
    param.gamma = dgamma;
    break;
  default:
    param.kernel_type = LINEAR;
    break;
  }
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
  //estimate a reliable estimate for parameter C if not specified
  if(C == -1.0){
    double avgC = 0;
    for(unsigned int i = 0; i < windows.size(); i++){
      double dotsum = 0;
      for(unsigned int j = 0; j < windows[i].size(); j++) 
	dotsum += windows[i][j]*windows[i][j];
      avgC += dotsum;
    }
    avgC /= windows.size();
    param.C = 1/avgC;
  }
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
  x_space = Malloc(struct svm_node,elements);
  unsigned int max_index = windows[0].size();
  int j=0;
  for(unsigned int i=0;i<windows.size();i++){
    prob.x[i] = &x_space[j];
    //fill data in x_space and prob.y[i]
    for(int k = 0; k < (int) max_index-1; k++){
      x_space[j].index =  k+1;
      x_space[j].value = (double) windows[i][k];
      j++;
    }
    prob.y[i] = windows[i][max_index-1];
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
