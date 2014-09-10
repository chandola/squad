#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <math.h>

#include "sax.h"
using namespace std;
#define MININT -MAXINT;

void TStoPAA(vector<vector<float> > &sequences,vector<vector<float> > &paasequences,int segment){
  if(segment <= 0){cerr<<"Error: Segment length should be greater than 0\n";return;}
  paasequences.resize(0);
  for(unsigned int i=0; i < sequences.size(); i++){
    vector<float> paa(0);
    int numsegs = (int) ceil((float) sequences[i].size()/segment);
    for(int j = 0; j < numsegs; j++){
      unsigned int beg = j*segment;
      unsigned int end = beg+segment;
      if(end > sequences[i].size()) end = sequences[i].size();
      float sum = 0;
      for(unsigned int k = beg; k < end; k++) sum += sequences[i][k];
      sum = sum/(end-beg);
      paa.push_back(sum);
    }
    paasequences.push_back(paa);
  }
}

void PAAtoSAX(vector<vector<float> > &paasequences,vector<vector<int> > &symsequences,int alphasize){
  vector<float> breakpoints;getBreakpoints(breakpoints,alphasize);
  for(unsigned int i = 0; i < paasequences.size(); i++){
    vector<int> vec(0);
    for(unsigned int j = 0; j < paasequences[i].size(); j++){
      vec.push_back(getalphabet(paasequences[i][j],breakpoints));
    }
    symsequences.push_back(vec);
  }
}

int getalphabet(float val, vector<float> &breakpoints){
  unsigned int ind = breakpoints.size();
  for(unsigned int i = 0; i < breakpoints.size(); i++){
    if(val <= breakpoints[i]){ind = i;break;}
  }
  return 1+ind;
}

void saveAlphaMap(int alphasize,ofstream &out){
  vector<float> breakpoints;
  getBreakpoints(breakpoints,alphasize);
  if(breakpoints.size() != 0){
    for(int i = 0; i < alphasize; i++){
      for(int j = 0; j < alphasize; j++){
	float sim = 0.0;
	if(abs(i-j) <= 1){
	  sim = 0.0;
	}else{
	  if(i > j) sim = breakpoints[i-1] - breakpoints[j];
	  else sim = breakpoints[j-1] - breakpoints[i];
	}
	sim = 1/(1+sim);
	int a = 1 + i;int b = 1 + j;
	out<<a<<","<<b<<","<<sim<<"\n";
      }
    }
  }
}

void getBreakpoints(vector<float> &breakpoints, int alphasize){
  breakpoints.resize(0);
  switch (alphasize){
  case 3:
    breakpoints.push_back(-0.43);breakpoints.push_back(0.43);
    break;
  case 4:
    breakpoints.push_back(-0.67);breakpoints.push_back(0);breakpoints.push_back(0.67);
    break;
  case 5:
    breakpoints.push_back(-0.84);breakpoints.push_back(-0.25);breakpoints.push_back(0.25);breakpoints.push_back(0.84);
    break;
  case 6:
    breakpoints.push_back(-0.97);breakpoints.push_back(-0.43);breakpoints.push_back(0);breakpoints.push_back(0.43);breakpoints.push_back(0.97);
    break;
  case 7:
    breakpoints.push_back(-1.07);breakpoints.push_back(-0.57);breakpoints.push_back(-0.18);breakpoints.push_back(0.18);breakpoints.push_back(0.57);breakpoints.push_back(1.07);
    break;
  case 8:
    breakpoints.push_back(-1.15);breakpoints.push_back(-0.67);breakpoints.push_back(-0.32);breakpoints.push_back(0);breakpoints.push_back(0.32);breakpoints.push_back(0.67);breakpoints.push_back(1.15);
    break;
  default:
    cerr<<"Error: Value "<<alphasize<<" not supported.\n";
    break;
  }
}
