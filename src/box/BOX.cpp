#include<iostream>
#include<fstream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <queue>
#include <map>

#include "tools.h"

#define INF 10E6
using namespace std;


class Space{
public:
  float x;
  float y;
  float z;
};

class Box{
public:
  vector<Space> elements;
  vector<int> points;
  float volume;
  float x[2], y[2],z[2];
  void computeVolume(){
    x[0] = y[0] = z[0] = 0; // MAX
    x[1] = y[1] = z[1] = INF; // MIN
		  
    for(int i = 0; i < elements.size(); i++){
      if(elements[i].x > x[0])
	x[0] = elements[i].x;		
      if(elements[i].x < x[1])
	x[1] = elements[i].x;
      if(elements[i].y > y[0])
	y[0] = elements[i].y;		
      if(elements[i].y < y[1])
	y[1] = elements[i].y;
      if(elements[i].z > z[0])
	z[0] = elements[i].z;		
      if(elements[i].z < z[1])
	z[1] = elements[i].z;
    }
    volume = (x[0]-x[1])*(y[0]-y[1])*(z[0]-z[1]);
  }

  void remove(){
    elements.resize(0);
    points.resize(0);
    volume = 0;
    x[0] = y[0] = z[0] = x[1] = y[1] = z[1] = 0;
  }
};


class compare_A{
public:
  bool operator() (pair<int,float> a, pair<int,float>  b){
    return a.second > b.second;
  }
};

float calcDistance_center(Space feature, Box box){
  float x1, y1, z1, dist;
  x1 = ((box.x[0]+box.x[1])/2-feature.x);
  y1 = ((box.y[0]+box.y[1])/2-feature.y);
  z1 = ((box.z[0]+box.z[1])/2-feature.z);
  dist = x1*x1 + y1*y1 + z1*z1;
  return dist;
}

float calcDistance(Space feature, Box box){
  // 0 : MAX ; 1 : MIN
  float x1=0.0, y1=0.0, z1=0.0, dist=0.0;
  if(feature.x <= box.x[0] && feature.x >= box.x[1])
    if(feature.y <= box.y[0] && feature.y >= box.y[1])
      if(feature.z <= box.z[0] && feature.z >= box.z[1])
	return 0;
  
  if(feature.z >= box.z[0])
    z1 = -box.z[0] + feature.z;
  else if(feature.z <= box.z[1])
    z1 = box.z[1] - feature.z;
	  
  if(feature.y >= box.y[0])
    y1 = -box.y[0] + feature.y;
  else if(feature.y <= box.y[1])
    y1 = box.y[1] - feature.y;
	  
  if(feature.x >= box.x[0])
    x1 = -box.x[0] + feature.x;
  else if(feature.x <= box.x[1])
    x1 = box.x[1] - feature.x;
	 
  dist = x1*x1 + y1*y1 + z1*z1;
  return dist;
}

pair<float, int> calcDistance(Space feature, Box box1, Box box2){

  float dist1 = calcDistance(feature, box1);
  float dist2 = calcDistance(feature, box2);
  pair<float, int> result;
  if(dist1 < dist2){
    result.first=dist1;
    result.second=1;	
    return result;
  }
  result.first=dist2;
  result.second=2;	
  return result;
}

bool isInBox(Space feature, Box box1){
  if(calcDistance(feature,box1) == 0)
    return true;
  return false;
}

Box merge(Box box1, Box box2){
	
  Box box_new;
  for(int i = 0; i < box1.elements.size(); i++){
    box_new.elements.push_back(box1.elements[i]);
    box_new.points.push_back(box1.points[i]);
  }
  for(int i = 0; i < box2.elements.size(); i++){
    box_new.elements.push_back(box2.elements[i]);
    box_new.points.push_back(box2.points[i]);
  }
  box_new.computeVolume();
  return box_new;
}

float computeVolumeDifference(Box box1, Box box2){

  Box box_new = merge(box1, box2);
  return (box_new.volume - box1.volume - box2.volume);
}

char  *trainfile = 0;
char   *testfile = 0;
char *outputfile = 0;
int  rand_starts = 0;
int            k = 0;
#define usage(msg, args...)			\
  {									\
    fprintf(stderr, msg, ##args);					\
    fprintf(stderr,							\
	    "Usage:\n\n"						\
	    "./BOX -i -trainfilename -t testfilename -k numboxes -r numrandomrestarts -o outputfilename\n\n" \
	    "-i      input training file name\n"			\
	    "-t      input testing file name\n"				\
	    "-o      output file\n"					\
	    "-k      number of boxes \n"				\
	    "-r      number of random restarts\n"			\
	    "\n");							\
    exit(1);								\
  }

void parseArgs(int argc, char** argv){
  int argi;
  for(argi=1; argi<argc; argi++){
    if(strcmp(argv[argi],"-i")==0){
      if(++argi<argc) trainfile=strdup(argv[argi]);
      else usage("-i requires input train file argument\n");
    } else if(strcmp(argv[argi],"-t")==0){
      if(++argi<argc) testfile=strdup(argv[argi]);
      else usage("-t requires input test file argument\n");
    } else if(strcmp(argv[argi],"-o")==0){
      if(++argi<argc) outputfile=strdup(argv[argi]);
      else usage("-o requires output file argument\n");
    } else if(strcmp(argv[argi],"-k")==0){
      if(++argi<argc) k=atoi(argv[argi]);
      else usage("-k requires number of boxes argument\n");
    } else if(strcmp(argv[argi],"-r")==0){
      if(++argi<argc) rand_starts=atoi(argv[argi]);
      else usage("-r requires number of random starts argument\n");
    } else if(strcmp(argv[argi], "-h")==0 || strcmp(argv[argi], "-?")==0){
      usage("Usage:");
    }
    else
      usage("Invalid switch '%s'\n", argv[argi]);
  }
  if(!trainfile) usage("No input training file specified\n");
  if(!testfile) usage("No input testing file specified\n");
  if(!outputfile) usage("No output file specified\n");
  if(k == 0) usage("Number of boxes not specified\n");
  if(rand_starts == 0) usage("Number of random restarts not specified\n");
}

int main(int argc, char *argv[]){
  
  parseArgs(argc,argv);
  if(argc < 6){
    exit(0);
  }
  vector<vector<float> > train, test;
  ifstream in;
  in.open(strdup(trainfile));
  readSequences(train,in);
  in.close();
  in.open(strdup(testfile));
  readSequences(test,in);
  in.close();

  int trainSize = train.size();
  int testSize = test.size();
  int dimensions = train[0].size();

  for(int i = 0; i < train.size(); i++)
    if(train[i].size() != dimensions){
      cout<<"Train Dimension Mismatch "<<dimensions<<endl;
      return 0;
    }
  for(int i =0; i < test.size(); i++)
    if(test[i].size() != dimensions){
      cout<<"Test Dimension Mismatch "<<dimensions<<endl;
      return 0;
    }  
  //ofstream fout2("features_train");
  //Creating 3D feature space
  float sum_x=0.0, sum_y=0.0, sum_z=0.0, std_x=0.0, std_y=0.0, std_z=0.0;
  vector<vector<Space> > features;
  for(int i = 0; i < trainSize; i++){
    vector<Space> tmp;
    for(int j = 0; j < dimensions; j++){
      Space oneFeature;
      if(j > 0){
	oneFeature.x = train[i][j];	
	oneFeature.y = train[i][j] - train[i][j-1];
	oneFeature.z = oneFeature.y - tmp[j-1].y;
	sum_x+=oneFeature.x;
	sum_y+=oneFeature.y;
	sum_z+=oneFeature.z;
      }
      else{
	sum_x = train[i][j];
	sum_y = train[i][j];
	sum_z = train[i][j]; 
	oneFeature.x = train[i][j];
	oneFeature.y = train[i][j];
	oneFeature.z = train[i][j];
      }
      tmp.push_back(oneFeature);
    }
    int n=dimensions;
    for(int j = 0; j < dimensions; j++){
      std_x+=(tmp[j].x-sum_x/n)*(tmp[j].x-sum_x/n);
      std_y+=(tmp[j].y-sum_y/n)*(tmp[j].y-sum_y/n);
      std_z+=(tmp[j].z-sum_z/n)*(tmp[j].z-sum_z/n);
    }
    std_x = sqrt(std_x/n);
    std_y = sqrt(std_y/n);
    std_z = sqrt(std_z/n);
    for(int j = 0; j < dimensions; j++){
      tmp[j].x = (tmp[j].x-sum_x/n)/std_x;
      tmp[j].y = (tmp[j].y-sum_y/n)/std_y;
      tmp[j].z = (tmp[j].z-sum_z/n)/std_z;
    }
    features.push_back(tmp);
    sum_x=sum_y=sum_z=0;
    std_x=std_y=std_z=0;
  }
  train.clear();
  vector<vector<Space> > test_features;

  for(unsigned int i = 0; i < test.size(); i++){
    vector<Space> tmp;
    for(int j = 0; j < dimensions; j++){
      Space oneFeature;
      if(j > 0){	
	oneFeature.x = test[i][j];	
	oneFeature.y = test[i][j] - test[i][j-1];
	oneFeature.z = oneFeature.y - tmp[j-1].y;
	sum_x+=oneFeature.x;
	sum_y+=oneFeature.y;
	sum_z+=oneFeature.z;
      }
      else{
	sum_x = test[i][j];
	sum_y = test[i][j];
	sum_z = test[i][j];
	oneFeature.x = test[i][j];
	oneFeature.y = test[i][j];
	oneFeature.z = test[i][j];
      }
      tmp.push_back(oneFeature);
    }
    int n=dimensions;
    for(int j = 0; j < dimensions; j++){
      std_x+=(tmp[j].x-sum_x/n)*(tmp[j].x-sum_x/n);
      std_y+=(tmp[j].y-sum_y/n)*(tmp[j].y-sum_y/n);
      std_z+=(tmp[j].z-sum_z/n)*(tmp[j].z-sum_z/n);
    }
    std_x = sqrt(std_x/n);
    std_y = sqrt(std_y/n);
    std_z = sqrt(std_z/n);
    for(int j = 0; j < dimensions; j++){
      tmp[j].x = (tmp[j].x-sum_x/n)/std_x;
      tmp[j].y = (tmp[j].y-sum_y/n)/std_y;
      tmp[j].z = (tmp[j].z-sum_z/n)/std_z;
    }
    test_features.push_back(tmp);
    sum_x=sum_y=sum_z=0;
    std_x=std_y=std_z=0;
  }
  test.clear();
  //cout<<"Features Calculated"<<endl;
  //Smoothing function : optional
  int seq[rand_starts];
  for(int i = 0; i < rand_starts; i++){
    int p = rand();
    float x = (float) p/RAND_MAX;
    //cout<<x<<" "<<trainSize<<endl;
    seq[i] = (int) (x*trainSize);    
    //cout<<seq[i]<<endl;
  }

  float finalScore[testSize];  
  for(int i=0; i<testSize; i++)
    finalScore[i]=0;

  for(int s = 0; s < rand_starts; s++){
    
    //First step : Creating N-1 boxes with consecutive elements merged.
    int noBoxes = dimensions - 1;		
    Box* boxes = new Box[noBoxes];
    //cout<<"Start sequence: "<<seq[s]<<endl;
    for(int j = 0; j < dimensions-1; j++){
      boxes[j].elements.push_back(features[seq[s]][j]);
      boxes[j].elements.push_back(features[seq[s]][j+1]);
      boxes[j].points.push_back(j);
      boxes[j].points.push_back(j+1);
      boxes[j].computeVolume();
    }
    // cout<<"Number of boxes created : "<<noBoxes<<endl; 	 

    // GreedySplit. First fill the V_diff queue by calculating the Volume difference when consecutive boxes are merged.
    pair<int , float> pr;	
    priority_queue<pair<int , float> , vector<pair<int, float> >, compare_A> V_diff;
    float diff, diff1, diff2;
    for(int j = 0; j < noBoxes-1; j++){
      diff = computeVolumeDifference(boxes[j] , boxes[j+1]);
      //cout<<" Diff "<<diff<<endl;
      if(V_diff.size() < noBoxes){
	pr.first = j;
	pr.second = diff;
	V_diff.push(pr);
      }
    }

    // Reduce the noBoxes to k by merging the closest consecutive ones and making updates to V_diff queue. 
    int mergeIndex;
    while( noBoxes > k ){
      mergeIndex = V_diff.top().first;
      boxes[mergeIndex] = merge(boxes[mergeIndex], boxes[mergeIndex+1]); //update the box contents with all the entries 
      // Delete box[mergeIndex+1] and update the next boxes indices so that total number of boxes decrease by one.
      for(int m = mergeIndex+1; m < noBoxes-1; m++)
	boxes[m] = boxes[m+1];
      boxes[noBoxes-1].remove();
      V_diff.pop();
      noBoxes--;

      /* Update 2 old entries of V_diff which are diff(mergeIndex,mergeIndex-1) and diff(mergeIndex +1 , mergeIndex +2) which will become diff(mergeIndex , mergeIndex +1) after updating entries in deleteBox function.*/
      if(mergeIndex > 0)
	diff1 = computeVolumeDifference(boxes[mergeIndex - 1], boxes[mergeIndex]);
      if(mergeIndex < noBoxes-1){
	diff2 = computeVolumeDifference(boxes[mergeIndex], boxes[mergeIndex+1]);
	pr.first = mergeIndex;
	pr.second = diff2;
	V_diff.push(pr);
      }

      // Update the V_diff here
      vector<pair<int, float> > prs;
      while(V_diff.size()!=0){	
	pr = V_diff.top();
	if(pr.first == mergeIndex -1)
	  pr.second = diff1;
	if(pr.first==mergeIndex+1){
	  V_diff.pop();
	  continue;
	}
	if(pr.first>mergeIndex)
	  pr.first--;
	prs.push_back(pr);
	V_diff.pop();
      }
      //cout<<" V_diff size "<<prs.size()<<endl;
      for(int l = 0; l < prs.size(); l++)
	V_diff.push(prs[l]);
      prs.clear();
    }
    
    //cout<<" Final number of Boxes "<<noBoxes<<endl;
    for(int l = 0; l < k; l++){
      int count =0;
      for(int t = 0; t < boxes[l].points.size(); t++){
	//cout<<boxes[l].points[t]<<" ";
	count++;
      }
    }

    // Box Expand : Order Dependent 
    for(int i = 0 ; i < trainSize; i++){
      int index[dimensions];
      if(i == seq[s]) continue;
      for(int j = 0; j < dimensions; j++){
	float minDist = INF;
	//cout<<j<<" : "<<endl;
	for(int l = 0; l < k; l++){
	  float dist = calcDistance_center(features[i][j], boxes[l]);
	  if(dist < minDist){
	    minDist = dist;
	    index[j] = l;
	  }
	  //cout<<endl;
	}		
      }
      for(int j = 0; j < dimensions; j++){
	//cout<<j<<" "<<index[j]<<endl;
	boxes[index[j]].elements.push_back(features[i][j]);
      }
       
      for(int l = 0; l < k; l++){
	boxes[l].computeVolume();
      }
    }
    // Testing :
    float individualAnomalyScore[testSize];
    for(int dim=0; dim<testSize; dim++)
      individualAnomalyScore[dim]=0; 
    for(int t = 0; t < testSize; t++){
      int curr_box;
      float anomalyScore[dimensions];
      for(int p = 0; p < dimensions; p++){
	bool flag = false;
	anomalyScore[p] = 1;	  
	// Testing the box in which point 0 fits as all the other boxes depend on the previous box locations.	  
	if(p == 0){			
	  for(int l = 0; l < k; l++){
	    if(isInBox(test_features[t][p],boxes[l])){	
	      flag = true;  
	      curr_box = l;
	      anomalyScore[p] = 0;
	      break;
	    }
	  }
	  if(!flag){
	    float dist, min=INF;
	    for(int l = 0; l < k; l++){   
	      dist = calcDistance(test_features[t][p], boxes[l]);
	      if(dist < min){
		min = dist;
		curr_box = l;
	      }
	    }
	    anomalyScore[p] = dist;
	    //cout<<anomalyScore[p]<<endl;
	  }
	}
	else{
	  pair<float, int> pr;
	  //Testing if the point lies in current box or the next box, then its anomalyScore is 0;      
	  if(isInBox(test_features[t][p], boxes[curr_box]))
	    anomalyScore[p] = 0;
	  else if(curr_box+1<k) {
	    if(isInBox(test_features[t][p], boxes[curr_box+1])){
	      anomalyScore[p] = 0;
	      curr_box=curr_box+1;
	    }
	    else{
	      pr = calcDistance(test_features[t][p], boxes[curr_box], boxes[curr_box+1]);
	      if(pr.second==2)
		curr_box=curr_box+1;
	      anomalyScore[p]=pr.first;
	    }
	  }
	  else{
	    if(isInBox(test_features[t][p], boxes[0])){
	      anomalyScore[p] = 0;
	      curr_box=0;
	    }
	    else{
	      pr = calcDistance(test_features[t][p], boxes[curr_box], boxes[0]);
	      if(pr.second==2)
		curr_box=0;
	      anomalyScore[p]=pr.first;
	    }
	  }
	} 	     
	individualAnomalyScore[t] += anomalyScore[p];	
	finalScore[t] += anomalyScore[p];
      }
    }
    delete [] boxes;
  }
  ofstream fout(outputfile);
  for(int n = 0; n < testSize; n++){
    fout<<finalScore[n]<<endl;
  }
  fout.close();
  cout<<"Done"<<endl;
  return 0;
}
