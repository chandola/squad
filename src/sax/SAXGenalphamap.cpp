#include<iostream>
#include<fstream>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "sax.h"
using namespace std;
int main(int argc, char *argv[]){
  if(argc < 3){
    cerr<<"USAGE: ./SAXGenalphamap outfile alphasize\n";
    exit(0);
  }
  char* outfilename = strdup(argv[1]);
  int alphasize = atoi(argv[2]);
  ofstream out(outfilename);
  saveAlphaMap(alphasize,out);
  out.close();
}
