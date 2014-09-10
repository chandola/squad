/*
**      Author: Varun Chandola, chandola@cs.umn.edu
**      Date:   4 January 2009 
**      File:   esthmmall.c
**      Purpose: Wrapper for esthmm.c
**      Organization: University of Minnesota
**
*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"
#include <sys/types.h>
#include <unistd.h>

/* Modify the MAXLENGTH to handle sequences longer than 2000.
   Make sure to compile with -D__LONG flag.
   Default compilation has maximum length set to 500.
*/
#ifdef __LONG
#define MAXLENGTH 2000
#else
#define MAXLENGTH 500
#endif

void Usage(char *name);

int main (int argc, char **argv)
{
  int 	  T,N,M,*O, nflg=0, mflg=0, errflg =0,iflg = 0,oflg=0,c,seed,niter;
  HMM  	  hmm;
  double  **alpha, **beta, **gamma,logprobinit,logprobfinal; 
  char	  *inputfile,*hmmoutputfile;
  FILE	*fp;
  extern char *optarg;
  extern int optind, opterr, optopt;

  while ((c= getopt(argc, argv, "h:i:o:n:m:")) != EOF)
    switch (c) {
    case 'h': 
      Usage(argv[0]);
      exit(1);
      break;
    case 'i':
      iflg++;  
      inputfile = argv[optind-1];
      break;   
    case 'o':
      oflg++;  
      hmmoutputfile  = argv[optind-1];
      break;   
    case 'n':
      nflg++;  
      sscanf(optarg, "%d", &N);
      break;   
    case 'm':  
      mflg++;  
      sscanf(optarg, "%d", &M);
      break;   
    case '?':
      errflg++;
    }
  if ((!nflg) || (!mflg) || (!iflg) || (!oflg)) errflg++;     /* both N and M should be specified */
  if (errflg) {
    Usage(argv[0]);
    exit (1);
  }
  
  /* open the observed sequences file */
  fp = fopen(inputfile, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found \n", inputfile);
    exit (1);
  }
  
  /* initialize the hmm model */
  seed = hmmgetseed();
  InitHMM(&hmm, N, M, seed);

  while(1){
    if(ReadSequence1(fp, &T, &O, MAXLENGTH) == -1) break;
    /* allocate memory */
    alpha = dmatrix(1, T, 1, hmm.N);
    beta = dmatrix(1, T, 1, hmm.N);
    gamma = dmatrix(1, T, 1, hmm.N);
    /* call Baum Welch */
    BaumWelch(&hmm,T,O,alpha,beta,gamma,&niter,&logprobinit,&logprobfinal);
    /* free memory */
    free_ivector(O, 1, T);
    free_dmatrix(alpha, 1, T, 1, hmm.N);
    free_dmatrix(beta, 1, T, 1, hmm.N);
    free_dmatrix(gamma, 1, T, 1, hmm.N);
  }
  fclose(fp);
  /* save the HMM */
  fp = fopen(hmmoutputfile,"w");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found \n", hmmoutputfile);
    exit (1);
  }
  PrintHMM(fp, &hmm);
  fclose(fp);
  FreeHMM(&hmm);
  return 0;
}

void Usage(char *name)
{
  printf("Usage error. \n");
  printf("Usage: %s -i <input sequences> -o <output hmm model file> -n <num_states> -m <max_alphabet>\n", 
	 name);
  printf("  i - file containing the input seqences\n");
  printf("  o - file for saving hmm\n");
  printf("  n - number of states\n");
  printf("  m - maximum alphabet\n");
}
