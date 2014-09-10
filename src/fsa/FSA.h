//fsa.h - Varun Chandola

#include <vector>
#include <map>
#include <fstream>

using namespace std;

class node{
public:
  int id;
  map<pair<int,int>,int> transitions;
  int numedges;

  node(){}
  node(int _id){
    id = _id;
    transitions.clear();
    numedges = 0;
  }
};

void readSequences(vector<vector<vector<int>  > > &windows,
		   int max_alphabet,
		   char *seqfile,
		   int windowsize);

void trainSequences(map<int,node> &nodesMap,
		    map<vector<int>  ,int> &str2intnodemap,
		    map<int,vector<int>  > &int2strnodemap,
		    map<vector<int>  ,int> &str2intedgemap,
		    map<int,vector<int>  > &int2stredgemap,
		    vector<vector<vector<int>  > > &windows,
		    int n,int l, int *totaledges);

void predictSequences(map<int,node> &nodesMap,
		      map<vector<int>  ,int> &str2intnodemap,
		      map<int,vector<int>  > &int2strnodemap,
		      map<vector<int>  ,int> &str2intedgemap,
		      map<int,vector<int>  > &int2stredgemap,
		      vector<vector<vector<int>  > > &windows,
		      ofstream &out,
		      int n,int l,int totaledges,int z);

void predictSequencesStide(map<int,node> &nodesMap,
		      map<vector<int>  ,int> &str2intnodemap,
		      map<int,vector<int>  > &int2strnodemap,
		      map<vector<int>  ,int> &str2intedgemap,
		      map<int,vector<int>  > &int2stredgemap,
		      vector<vector<vector<int>  > > &windows,
		      ofstream &out,
		      int n,int l,int totaledges);

void printFSA(map<int,node> &nodesMap,
	      ofstream &out,
	      map<vector<int>  ,int> &str2intnodemap,
	      map<int,vector<int>  > &int2strnodemap,
	      map<vector<int>  ,int> &str2intedgemap,
	      map<int,vector<int>  > &int2stredgemap,
	      int n, int l,int totaledges);

void loadFSA(map<int,node> &nodesMap,
	     ifstream &in,
	     map<vector<int>,int> &str2intnodemap,
	     map<int,vector<int> > &int2strnodemap,
	     map<vector<int>,int> &str2intedgemap,
	     map<int,vector<int> > &int2stredgemap,
	     int *n, int *l,int *totaledges);
