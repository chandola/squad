#include<iostream>
#include<fstream>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "sax.h"
using namespace std;
int main(int argc, char *argv[]){
  if(argc < 5){
    cerr<<"USAGE: ./SAXConv infile outfile segmentlength alphasize\n";
    exit(0);
  }
  char* filename = strdup(argv[1]);
  char* outfilename = strdup(argv[2]);
  int segment = atoi(argv[3]);
  int alphasize = atoi(argv[4]);
  if((alphasize < 3) || (alphasize > 8)) {
    cerr<<"Alphabet size can be between 3 and 8 only.\n";
    exit(0);
  }
  ifstream in;
  in.open(filename);
  ofstream out(outfilename);
  vector<vector<float> > sequences(0);
  readSequences(sequences,in);
  in.close();
  vector<vector<float> > paasequences(0); 
  TStoPAA(sequences,paasequences,segment);
  vector<vector<int> > symsequences(0);
  PAAtoSAX(paasequences,symsequences,alphasize);
  printSequences(symsequences, out);
  out.close();
}
