#define MAXENTROPY      0
#define LEASTSQUARES    1
#define TRUE  1
#define FALSE 0

int AutoRegression(double *,int,int,double *,int);
int ARMaxEntropy(double *,int,int,double **,double *,double *,double *,double *);
int ARLeastSquare(double *,int,int,double *);
int SolveLE(double **,double *,unsigned int);

