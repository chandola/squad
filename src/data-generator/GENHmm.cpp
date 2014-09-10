#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

int states = 0;
int sc = 0;
char *outputfileprefix = 0;
float delta = 0.0;
float delta2 = 0.0;
float mu = 0.0;

#define usage(msg, args...) \
{\
 fprintf(stderr,"Error: ");\
  fprintf(stderr, msg, ##args);\
  fprintf(stderr, \
"Usage:\n" \
"./GENHmm -n numstates -o hmmfileprefix -s shortcircuit -d delta -d2 delta2 -m mu\n"\
"   -n  [4.. ] -> number of states in HMMs\n"\
"   -o         -> prefix for hmm model files.\n"\
"   -s  [1..n] -> shortcircuit parameter\n"\
"   -d  [0...1]-> delta transition probability within a cluster\n"\
"   -d2 [0..1] -> delta transition probability from one cluster to another\n"\
"   -m         -> mu or observation probability" \
"    \n"); \
  exit(1); \
}

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) outputfileprefix=strdup(argv[argi]);
      else usage("-o requires output file prefix argument\n");
    } else if(strcmp(argv[argi],"-n")==0){
      if(++argi<argc) states=atoi(argv[argi]);
      else usage("-n requires number of states argument\n");
    }else if(strcmp(argv[argi],"-s")==0){
      if(++argi<argc) sc=atoi(argv[argi]);
      else usage("-s requires short circuit argument\n");
    }else if(strcmp(argv[argi],"-d")==0){
      if(++argi<argc) delta=atof(argv[argi]);
      else usage("-d requires delta argument\n");
    }else if(strcmp(argv[argi],"-d2")==0){
      if(++argi<argc) delta2=atof(argv[argi]);
      else usage("-d2 requires delta2 argument\n");
    }else if(strcmp(argv[argi],"-m")==0){
      if(++argi<argc) mu=atof(argv[argi]);
      else usage("-m requires mu argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!outputfileprefix) usage("No output file prefix specified.");
  if(delta == 0.0) usage("No delta parameter specified.");
  if(delta2 == 0.0) usage("No delta2 parameter specified.");
  if(delta > 0.05){
    cerr<<"Warning: Delta parameter should be below 0.05. Reducing delta\n";
    delta = 0.05;
  }
  if(states == 0) usage("No shortcircuit parameter specified.");
  if(states%2 == 1){
    cerr<<"Warning: Number of states should be even. Decrementing states to "<<states-1<<"\n";
    states = states -1;
  }
  if(states < 4){
    cerr<<"Error: Number of states should be at least 4\n";
    exit(0);
  }
  if((sc > states/2) || (sc <= 0)){
    cerr<<"Error: Shortcircuit parameter should be between 1 and "<<states/2<<"\n";
    exit(0);
  }
  if(mu == 0.0) usage("No mu parameter specified.");
}

//this program generates an HMM that can be used to generate artificial sequence data.

void printmat(float **mat, int nrows, int ncols, std::ostream &out){
  for(int i = 0; i < nrows;i++){
    for(int j = 0; j < ncols-1; j++){
      out<<mat[i][j]<<" ";
    }
    out<<mat[i][ncols-1]<<"\n";
  } 
}

int main(int argc, char* argv[]){
  using namespace std;
  parseArgs(argc,argv);
  float delta1 = 10*delta;
  int alphas = states/2;
  float onestate = 1-(mu*(alphas-1));
  if(onestate < 0.5){
    cerr<<"The maximum observation frequency is too low. Reduce mu.\n";
    exit(0);
  }
  ofstream out1,out2,out3;
  char filename[32],filename1[32],filename2[32];
  sprintf(filename,"%s.norm.hmm",outputfileprefix);out1.open(filename);
  sprintf(filename1,"%s.anom1.hmm",outputfileprefix);out2.open(filename1);
  sprintf(filename2,"%s.anom2.hmm",outputfileprefix);out3.open(filename2);
  out1<<"M= "<<alphas<<"\n";out2<<"M= "<<alphas<<"\n";out3<<"M= "<<alphas<<"\n";
  out1<<"N= "<<states<<"\n";out2<<"N= "<<states<<"\n";out3<<"N= "<<states<<"\n";
  out1<<"A:"<<"\n";out2<<"A:"<<"\n";out3<<"A:\n";
  float **statetrans1,**statetrans2,**statetrans3;
  //allocate memory
  statetrans1 = (float **) malloc(states*sizeof(float *));
  statetrans2 = (float **) malloc(states*sizeof(float *));
  statetrans3 = (float **) malloc(states*sizeof(float *));
  for(int i = 0; i < states;i++){
    statetrans1[i] = (float *) malloc(states*sizeof(float));
    statetrans2[i] = (float *) malloc(states*sizeof(float));
    statetrans3[i] = (float *) malloc(states*sizeof(float));
  }
  float firsttrans = 1 - ((sc-1)*delta);
  for(int i = 0; i < states;i++){
    for(int j = 0; j < states; j++){
      statetrans1[i][j] = 0.0;
      statetrans2[i][j] = 0.0;
      statetrans3[i][j] = 0.0;
    }
    if(i<states/2){
      statetrans1[i][(i+1)%(states/2)] = firsttrans;
      statetrans2[i][(i+1)%(states/2)] = firsttrans;
      statetrans3[i][(i+1)%(states/2)] = firsttrans;
      for(int j = i+2; j < i+sc+1; j++){
	statetrans1[i][j%(states/2)] = delta;
	statetrans2[i][j%(states/2)] = delta;
	statetrans3[i][j%(states/2)] = delta;
      }
    }else{
      statetrans1[i][(states/2)+(i+1)%(states/2)] = firsttrans;
      statetrans2[i][(states/2)+(i+1)%(states/2)] = firsttrans;
      statetrans3[i][(states/2)+(i+1)%(states/2)] = firsttrans;
      for(int j = i+2; j < i+sc+1; j++){
	statetrans1[i][(states/2)+j%(states/2)] = delta;
	statetrans2[i][(states/2)+j%(states/2)] = delta;
	statetrans3[i][(states/2)+j%(states/2)] = delta;
      }
    }
  }
  statetrans1[1][states/2+1] = delta2;statetrans1[1][2%(states/2)] -= delta2;
  statetrans1[states/2][0] = statetrans1[states/2][(states/2)+1]-delta2;statetrans1[states/2][(states/2)+1] = delta2;
  printmat(statetrans1,states,states,out1);

  statetrans2[(states/2)+1][1] = delta;statetrans2[(states/2)+1][(states/2)+((states/2+2)%(states/2))] -= delta;
  statetrans2[0][states/2] = statetrans2[0][1]-delta;statetrans2[0][1] = delta;
  printmat(statetrans2,states,states,out2);

  statetrans3[1][states/2+1] = delta1;statetrans3[1][2%(states/2)] -= delta1;
  statetrans3[states/2][0] = statetrans3[states/2][(states/2)+1]-delta1;statetrans3[states/2][(states/2)+1] = delta1;
  printmat(statetrans3,states,states,out3);
  //output observation matrix
  out1<<"B:"<<"\n";out2<<"B:"<<"\n";out3<<"B:\n";
  float **obs;float **obs1;
  //allocate memory
  obs = (float **) malloc((states/2)*sizeof(float *));
  obs1 = (float **) malloc((states/2)*sizeof(float *));
  for(int i = 0; i < states/2;i++){
    obs[i] = (float *) malloc(alphas*sizeof(float));    
    obs1[i] = (float *) malloc(alphas*sizeof(float));
    for(int j = 0; j < alphas; j++){
      obs[i][j] = mu;
      obs1[i][j] = mu;
    }
  }
  //fill up the matrix
  for(int i = 0; i < states/2;i++) obs[i][i] = onestate;
  for(int i = 0; i < states/2;i++) obs1[i][alphas-i-1] = onestate;
  printmat(obs,states/2,alphas,out1);
  printmat(obs1,states/2,alphas,out1);
  printmat(obs,states/2,alphas,out2);
  printmat(obs1,states/2,alphas,out2);
  printmat(obs,states/2,alphas,out3);
  printmat(obs1,states/2,alphas,out3);
  //print out the initialization vector
  out1<<"pi:\n";  out2<<"pi:\n";  out3<<"pi:\n";
  float perstate = (float) 2/states;
  for(int i = 0; i < states/2;i++){
    out1<<perstate<<" ";
    out2<<"0 ";
    out3<<perstate<<" ";
  }
  for(int i = 0; i < (states/2)-1;i++){
    out1<<"0 ";
    out2<<perstate<<" ";
    out3<<"0 ";
  }
  out1<<"0\n";
  out2<<perstate<<"\n";
  out3<<"0\n";
  out1.close();
  out2.close();
  out3.close();
}
