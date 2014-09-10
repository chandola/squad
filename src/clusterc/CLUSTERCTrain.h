#include <fstream>
#include <string>
#include <vector>
#include <queue>

using namespace std;
void unifrnd(int N, unsigned int k, vector<unsigned int> &s);
void findClusters(vector<vector<float> > &sequences,vector<vector<float> > &centroids,vector<vector<unsigned int> > &assignments);
void findCentroids(vector<vector<float> > &sequences,vector<vector<float> > &centroids,vector<vector<unsigned int> > &assignments);
float objFun(vector<vector<float> > &sequences,vector<vector<float> > &centroids,vector<vector<unsigned int> > &assignments);
bool stoppingCondition(queue<float> &objvals);
