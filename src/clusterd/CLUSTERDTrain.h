#include <fstream>
#include <string>
#include <vector>
#include <queue>

using namespace std;
void unifrnd(int N, unsigned int k, vector<unsigned int> &s);
void findClusters(vector<vector<int> > &sequences,vector<unsigned int> &ids,vector<vector<unsigned int> > &assignments);
void findMedoids(vector<vector<int> > &sequences,vector<unsigned int> &ids,vector<vector<unsigned int> > &assignments);
void printClusters(vector<vector<int> > &sequences,vector<unsigned int> &ids,vector<vector<unsigned int> > &assignments,ofstream &out);
float objFun(vector<vector<int> > &sequences,vector<unsigned int> &ids,vector<vector<unsigned int> > &assignments);
bool stoppingCondition(queue<float> &objvals);
bool verifyInitMedoids(vector<vector<int> > &sequences,vector<unsigned int> &ids);
