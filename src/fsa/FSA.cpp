#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include "FSA.h"
#include "tools.h"

void readSequences(vector<vector<vector<int> > > &windows,
		   int max_alphabet,
		   char *seqfile, 
		   int windowsize){
  //start reading the sequences
  ifstream in;in.open(seqfile);
  if(in.fail()) {cerr << "sequence file does not exist\n"; exit(0);}
  string bufstr;
  int lineid=0;
  while(!in.fail()){
    lineid++;
    getline(in,bufstr);
    if((bufstr.length() == 0) || (bufstr[0] == '#')) continue;
    vector<int> vals;
    Tokenize(bufstr,vals);
    //check if the length of sequence is less than or equal to the specified window size
    if(vals.size() <= (unsigned int) windowsize){
      cerr<<"ERROR: Sequence "<<lineid<<" is shorter than the windowsize.\n";
      exit(0);
    }
    //check if any symbol is out of range, i.e., [1,max_alphabet]
    for(unsigned int i = 0; i < vals.size();i++){
      if((vals[i] < 1) || (vals[i] > max_alphabet)){
	cerr<<"ERROR: Sequence "<<lineid<<" contains an invalid symbol "<<vals[i]<<".\n";
	exit(0);
      }
    }
    //run a sliding window
    vector<vector<int> > tmp(0);
    vector<int>::iterator it,it1;
    for(it=vals.begin(),it1=vals.begin()+windowsize-1;it1 != vals.end();it++,it1++) tmp.push_back(vector<int>(it,it1));
    windows.push_back(tmp);
  }
  in.close();
}

void trainSequences(map<int,node> &nodesMap,
		    map<vector<int> ,int> &str2intnodemap,
		    map<int,vector<int> > &int2strnodemap,
		    map<vector<int> ,int> &str2intedgemap,
		    map<int,vector<int> > &int2stredgemap,
		    vector<vector<vector<int> > > &windows,
		    int n, int l,int *totaledges){
  if(windows.size() == 0){
    cerr<<"No windows. No training possible. Exiting ... \n";return;
  }
  int currEdgeID = 0;
  int currNodeID = 0;
  //find the largest id in int2strnodemap and int2stredgemap
  map<int,vector<int> >::iterator it2;
  for(it2 = int2strnodemap.begin();it2 != int2strnodemap.end();it2++){
    if(it2->first > currNodeID) currNodeID = it2->first;
  }
  for(it2 = int2stredgemap.begin();it2 != int2stredgemap.end();it2++){
    if(it2->first > currEdgeID) currEdgeID = it2->first;
  }
  for(unsigned int j = 0; j < windows.size(); j++){
    map<int,node>::iterator it;
    map<pair<int,int>,int>::iterator it1;
    vector<int> buf,buf1,buf2;
    int c_node,c_edge,n_node,n_edge;
    
    //read the first window
    buf = windows[j][0];
    buf1 = vector<int>(buf.begin(),buf.begin()+n-1);
    buf2 = vector<int>(buf.begin()+n-1,buf.end());
    if(str2intnodemap.find(buf1) != str2intnodemap.end()){
      c_node = str2intnodemap[buf1];
    }else{
      c_node = currNodeID;
      currNodeID++;
      str2intnodemap[buf1] = c_node;
      int2strnodemap[c_node] = buf1;
    }
    if(str2intedgemap.find(buf2) != str2intedgemap.end()){
      c_edge = str2intedgemap[buf2];
    }else{
      c_edge = currEdgeID;
      currEdgeID++;
      str2intedgemap[buf2] = c_edge;
      //check this
      int2stredgemap[c_edge] = buf2;
    }
    *totaledges += 1;
    
    if(nodesMap.find(c_node) == nodesMap.end()){
      //create a new node in the FSA
      node firstNode(c_node);
      nodesMap[c_node] = firstNode;
    }
    
     for(unsigned int i = 1; i < windows[j].size(); i++){
       buf = windows[j][i];
       buf1 = vector<int>(buf.begin(),buf.begin()+n-1);
       buf2 = vector<int>(buf.begin()+n-1,buf.end());
       if(str2intnodemap.find(buf1) != str2intnodemap.end()){
	 n_node = str2intnodemap[buf1];
       }else{
	 n_node = currNodeID;
	 currNodeID++;
	 str2intnodemap[buf1] = n_node;
	 int2strnodemap[n_node] = buf1;
       }
       if(str2intedgemap.find(buf2) != str2intedgemap.end()){
	 n_edge = str2intedgemap[buf2];
       }else{
	 n_edge = currEdgeID;
	 currEdgeID++;
	 str2intedgemap[buf2] = n_edge;
	 int2stredgemap[n_edge] = buf2;
       }
       *totaledges += 1;
      
       it = nodesMap.find(c_node);
       if(it == nodesMap.end()){cerr<<"Error.\n";exit(0);}
       it->second.numedges += 1;
       if(nodesMap.find(n_node) != nodesMap.end()){
 	it1 = (it->second).transitions.find(pair<int,int>(n_node,c_edge));
 	if(it1 != (it->second).transitions.end()){
 	  //scenario 1: n_node exists and a transition from c_node to n_node
 	  //with edge c_edge exists
 	  it1->second += 1;//increment the count
 	}else{
 	  //scenario 2: n_node exists and a transition from c_node to n_node
 	  //with edge c_edge does not exist
 	  (it->second).transitions[pair<int,int>(n_node,c_edge)] = 1;
 	}
       }else{
 	//scenario 3: n_node does not exist
 	node one_node(n_node);
 	nodesMap[n_node] = one_node;
 	(it->second).transitions[pair<int,int>(n_node,c_edge)] = 1;
       }
       c_node = n_node;
       c_edge = n_edge;
     }
 }
}

void predictSequences(map<int,node> &nodesMap,
		      map<vector<int> ,int> &str2intnodemap,
		      map<int,vector<int> > &int2strnodemap,
		      map<vector<int> ,int> &str2intedgemap,
		      map<int,vector<int> > &int2stredgemap,
		      vector<vector<vector<int> > > &windows,
		      ofstream &out,
		      int n,int l,int totaledges,int z){

  for(unsigned int j = 0; j < windows.size(); j++){
    vector<float> transitionProbabilities;
    
    vector<int> buf,buf1,buf2;
    int c_node,c_edge,n_node,n_edge;
    int lines = 0;
    buf = windows[j][0];
    buf1 = vector<int>(buf.begin(),buf.begin()+n-1);
    buf2 = vector<int>(buf.begin()+n-1,buf.end());
    
    if(str2intnodemap.find(buf1) != str2intnodemap.end())	c_node = str2intnodemap[buf1];
    else c_node = -1;
    if(str2intedgemap.find(buf2) != str2intedgemap.end())	c_edge = str2intedgemap[buf2];
    else c_edge = -1;
    lines++;
    
    for(unsigned int i = 1; i < windows[j].size(); i++){
      buf = windows[j][i];
      buf1 = vector<int>(buf.begin(),buf.begin()+n-1);
      buf2 = vector<int>(buf.begin()+n-1,buf.end());
      if(str2intnodemap.find(buf1) != str2intnodemap.end())      n_node = str2intnodemap[buf1];
      else    n_node = -1;
      
      if(str2intedgemap.find(buf2) != str2intedgemap.end())      n_edge = str2intedgemap[buf2];
      else    n_edge = -1;
      if(c_node != -1){
	//a sanity check
	if(nodesMap.find(c_node) == nodesMap.end()){cerr<<"The node is not in the map\n";exit(0);}
	if(nodesMap[c_node].transitions.find(pair<int,int>(n_node,c_edge)) != nodesMap[c_node].transitions.end()){
	  //scenario 1: c_node exists and transition to <n_node,c_edge> exists
	  int freq = nodesMap[c_node].transitions[pair<int,int>(n_node,c_edge)];
	  float tprob = (float) freq/nodesMap[c_node].numedges;
	  transitionProbabilities.push_back(tprob);
	}else{
	  //scenario 2: c_node exists but transition to <n_node,c_edge> does not exist
	  transitionProbabilities.push_back(0.0);
	}             
      }else{
	//scenario 3: c_node does not exist
	transitionProbabilities.push_back((float) z);
      }
      c_node = n_node;
      c_edge = n_edge;
    }
    
    //write out the predictions
    for(unsigned int i = 0; i < transitionProbabilities.size()-1; i++) out<<transitionProbabilities[i]<<" ";
    out<<transitionProbabilities[transitionProbabilities.size()-1]<<"\n";
  }
}

void predictSequencesStide(map<int,node> &nodesMap,
			   map<vector<int> ,int> &str2intnodemap,
			   map<int,vector<int> > &int2strnodemap,
			   map<vector<int> ,int> &str2intedgemap,
			   map<int,vector<int> > &int2stredgemap,
			   vector<vector<vector<int> > > &windows,
			   ofstream &out,
			   int n,int l,int totaledges){

  for(unsigned int j = 0; j < windows.size(); j++){
    vector<float> nodeProbabilities;
    
    vector<int> buf,buf1,buf2;
    int c_node,c_edge,n_node,n_edge;
    int lines = 0;
    //start reading the sequences
    
    buf = windows[j][0];
    buf1 = vector<int>(buf.begin(),buf.begin()+n-1);
    buf2 = vector<int>(buf.begin()+n-1,buf.end());
    
    if(str2intnodemap.find(buf1) != str2intnodemap.end())	c_node = str2intnodemap[buf1];
    else c_node = -1;
    if(str2intedgemap.find(buf2) != str2intedgemap.end())	c_edge = str2intedgemap[buf2];
    else c_edge = -1;
    lines++;
    
    for(unsigned int i = 1; i < windows[j].size(); i++){
      buf = windows[j][i];
      buf1 = vector<int>(buf.begin(),buf.begin()+n-1);
      buf2 = vector<int>(buf.begin()+n-1,buf.end());
      if(str2intnodemap.find(buf1) != str2intnodemap.end())      n_node = str2intnodemap[buf1];
      else    n_node = -1;
      
      if(str2intedgemap.find(buf2) != str2intedgemap.end())      n_edge = str2intedgemap[buf2];
      else    n_edge = -1;
      if(c_node != -1){
	//a sanity check
	if(nodesMap.find(c_node) == nodesMap.end()){cerr<<"The node is not in the map\n";exit(0);}
	if(nodesMap[c_node].transitions.find(pair<int,int>(n_node,c_edge)) != nodesMap[c_node].transitions.end()){
	  //scenario 1: c_node exists and transition to <n_node,c_edge> exists
	  int freq = nodesMap[c_node].transitions[pair<int,int>(n_node,c_edge)];
	  float nprob = (float) freq;///totaledges;
	  nodeProbabilities.push_back(nprob);
	}else{
	  //scenario 2: c_node exists but transition to <n_node,c_edge> does not exist
	  nodeProbabilities.push_back(0.0);
	}             
      }else{
	//scenario 3: c_node does not exist
	nodeProbabilities.push_back(0.0);
      }
      c_node = n_node;
      c_edge = n_edge;
    }
    
    //write out the predictions
    for(unsigned int i = 0; i < nodeProbabilities.size()-1; i++) out<<nodeProbabilities[i]<<" ";
    out<<nodeProbabilities[nodeProbabilities.size()-1]<<"\n";
  }
}

void printFSA(map<int,node> &nodesMap,
	      ofstream &out,
	      map<vector<int>,int> &str2intnodemap,
	      map<int,vector<int> > &int2strnodemap,
	      map<vector<int> ,int> &str2intedgemap,
	      map<int,vector<int> > &int2stredgemap,
	      int n, int l, int totaledges){
  map<int,node>::iterator it;
  map<pair<int,int>,int>::iterator it1;
  out<<n<<" "<<l<<" "<<totaledges<<"\n";
  for(it = nodesMap.begin();it != nodesMap.end(); it++){
    printSequence(int2strnodemap[it->first],out);
    out<<endl;
    out<<(it->second).transitions.size()<<":"<<it->second.numedges<<"\n";
    for(it1 = (it->second).transitions.begin(); it1 != (it->second).transitions.end(); it1++){
      
      printSequence(int2stredgemap[(it1->first).second],out);
      out<<":";
      printSequence(int2strnodemap[(it1->first).first],out);
      out<<":"<<it1->second<<"+";
    }
    out<<"\n";
  }
}

void loadFSA(map<int,node> &nodesMap,
	     ifstream &in,
	     map<vector<int> ,int> &str2intnodemap,
	     map<int,vector<int> > &int2strnodemap,
	     map<vector<int> ,int> &str2intedgemap,
	     map<int,vector<int> > &int2stredgemap,
	     int *n, int *l, int *totaledges){
  //load an existing FSA from a model file

  int currEdgeID = 0;int currNodeID = 0;
  //find the largest id in int2strnodemap and int2stredgemap
  map<int,vector<int> >::iterator it2;
  for(it2 = int2strnodemap.begin();it2 != int2strnodemap.end();it2++){if(it2->first > currNodeID) currNodeID = it2->first;}
  for(it2 = int2stredgemap.begin();it2 != int2stredgemap.end();it2++){if(it2->first > currEdgeID) currEdgeID = it2->first;}
  map<int,node>::iterator it;
  string bufstr;getline(in,bufstr);
  vector<int> vals;
  Tokenize(bufstr,vals);
  *n = vals[0];
  *l = vals[1];
  *totaledges += vals[2];
  while(1){
    string bufstr;
    //first line: name of the node
    getline(in,bufstr);
    if(in.fail()) break;
    vector<int> vals1;
    Tokenize(bufstr,vals1);
    if(str2intnodemap.find(vals1) == str2intnodemap.end()){
      //create a new node
      node one_node(currNodeID);
      nodesMap[currNodeID] = one_node;
      str2intnodemap[vals1] = currNodeID;
      int2strnodemap[currNodeID] = vals1;
      currNodeID++;
    }
    int c_node = str2intnodemap[vals1];
    getline(in,bufstr);
    string::size_type ind = bufstr.find(":",0);
    if(ind == string::npos){
      cerr<<"Code 0 : Incorrect format of model file. Exiting ... \n";
      exit(0);
    }
    int numdistinctedges = atoi(bufstr.substr(0,ind).c_str());
    int numedges = atoi(bufstr.substr(ind+1).c_str());
    it = nodesMap.find(c_node);
    it->second.numedges = numedges;

    getline(in,bufstr);
    int currind = 0;
    for(int i = 0; i < numdistinctedges-1; i++){
      int newind = bufstr.find("+",currind);
      if(newind == -1){	cerr<<"Code 1 : Incorrect format of model file. Exiting ... \n";	exit(0);}
      string str = bufstr.substr(currind,newind-currind);
      string ed;string nd;vector<int> edvals;vector<int> ndvals;int ed_i;int nd_i;int ct = 0;
      int newind1 = str.find(":",0);
      if(newind1 == -1){	cerr<<"Code 2 : Incorrect format of model file. Exiting ... \n";	exit(0);}
      ed = str.substr(0,newind1);Tokenize(ed,edvals);
      int newind2 = str.find(":",newind1+1);
      if(newind2 == -1){	cerr<<"Code 3 : Incorrect format of model file. Exiting ... \n";	exit(0);}
      nd = str.substr(newind1+1,newind2-newind1-1);Tokenize(nd,ndvals);
      ct = atoi(str.substr(newind2+1).c_str());
      //find the index corresponding to ed and nd
      if(str2intedgemap.find(edvals) == str2intedgemap.end()) {
	ed_i = currEdgeID;str2intedgemap[edvals] = ed_i;int2stredgemap[ed_i] = edvals;currEdgeID++;
      }else{
	ed_i = str2intedgemap[edvals];
      }
      if(str2intnodemap.find(ndvals) == str2intnodemap.end()) {
	nd_i = currNodeID;str2intnodemap[ndvals] = nd_i;int2strnodemap[nd_i] = ndvals;currNodeID++;
	node one_node(nd_i);  nodesMap[nd_i] = one_node;
      }else{
	nd_i = str2intnodemap[ndvals];
      }

      pair<int,int> tmp(nd_i,ed_i);
      if((it->second).transitions.find(tmp) == (it->second).transitions.end()){
	(it->second).transitions[tmp] = ct;
      }else{
	(it->second).transitions[tmp] += ct;
      }
      currind = newind+1;
    }
    
    if(numedges != 0){
      string str = bufstr.substr(currind);
      string ed;string nd;vector<int> edvals;vector<int> ndvals;int ed_i;int nd_i;int ct = 0;
      int newind1 = str.find(":",0);
       if(newind1 == -1){	cerr<<"Code 4 : Incorrect format of model file. Exiting ... \n";exit(0);}
       ed = str.substr(0,newind1);Tokenize(ed,edvals);
       int newind2 = str.find(":",newind1+1);
       if(newind2 == -1){	cerr<<"Code 5 : Incorrect format of model file. Exiting ... \n";exit(0);}
       nd = str.substr(newind1+1,newind2-newind1-1);Tokenize(nd,ndvals);
       ct = atoi(str.substr(newind2+1).c_str());
       //find the index corresponding to ed and nd
       if(str2intedgemap.find(edvals) == str2intedgemap.end()) {
 	ed_i = currEdgeID;str2intedgemap[edvals] = ed_i;int2stredgemap[ed_i] = edvals;currEdgeID++;
       }else{
 	ed_i = str2intedgemap[edvals];
       }
       if(str2intnodemap.find(ndvals) == str2intnodemap.end()) {
	 nd_i = currNodeID;str2intnodemap[ndvals] = nd_i;int2strnodemap[nd_i] = ndvals;currNodeID++;
	 node one_node(nd_i);  nodesMap[nd_i] = one_node;
       }else{
	 nd_i = str2intnodemap[ndvals];
       }
 
       pair<int,int> tmp(nd_i,ed_i);
       if((it->second).transitions.find(tmp) == (it->second).transitions.end()){
 	(it->second).transitions[tmp] = ct;
       }else{
 	(it->second).transitions[tmp] += ct;
       }
    }
  }
}
