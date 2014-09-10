#include<vector>
#include<fstream>
using namespace std;
void TStoPAA(vector<vector<float> > &sequences,vector<vector<float> > &paasequences,int segment);
void PAAtoSAX(vector<vector<float> > &paasequences,vector<vector<int> > &symsequences,int alphasize);
int getalphabet(float val, vector<float> &breakpoints);
void getBreakpoints(vector<float> &breakpoints, int alphasize);
void saveAlphaMap(int alphasize,ofstream &out);
