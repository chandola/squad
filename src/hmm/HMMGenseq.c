/*
**      Author: Varun Chandola, chandola@cs.umn.edu
**      Date:   06 January 2009
**      File:   HMMGenseq.c
**      Purpose: driver for testing the Viterbi code.
**      Organization: University of Minnesota
*/

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"
#include <sys/types.h>
#include <unistd.h> 

void Usage(char *name);

int main (int argc, char **argv)
{
	HMM  	hmm; 	/* the HMM */
	int  	T; 	/* length of observation sequence */
	int	*O;	/* the observation sequence O[1..T]*/
	int	*q; 	/* the state sequence q[1..T] */
	FILE 	*fp,*fp1;
	char    *inputfile,*outputfile;
	int	iflg = 0,oflg = 0,rflg=0, tflg=0,nflg =0, errflg = 0;
	int	randomize=0,seed=1;	/* random number seed */
	int	n,c,i;
	extern char *optarg;
	extern int optind, opterr, optopt;
	while ((c= getopt(argc, argv, "h:i:o:r:t:n:")) != EOF)
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
	    outputfile  = argv[optind-1];
	    break;   
	  case 'r':
	    rflg++;  
	    sscanf(optarg, "%d", &randomize);
	    break;   
	  case 't':  
	    tflg++;  
	    sscanf(optarg, "%d", &T);
	    break;   
	  case 'n':  
	    nflg++;  
	    sscanf(optarg, "%d", &n);
	    break;   
	  case '?':
	    errflg++;
	  }
	if (errflg || !tflg || !iflg || !oflg||!nflg) {
	  Usage(argv[0]);
	  exit(1);
	}
	/* read HMM file */
	fp = fopen(inputfile,"r");		
	if (fp == NULL) {
	  fprintf(stderr, "Error: File %s not found \n", argv[optind]);
	  exit (1);
	}
	ReadHMM(fp, &hmm);
	fclose(fp);
	fp1 = fopen(outputfile,"w");
	int maxdeviation = (int) T/10;
	int randnum = 0;
	for(i=0;i<n;i++){
	  int T1 = T;
	  seed = (unsigned)time( NULL )+i;
	  if(randomize) {
	    srand(seed);
	    unsigned long r = rand(); r <<= 15; r += rand();
	    int randnum = (r%(2*maxdeviation))-maxdeviation;	    
	    T1 = T + randnum;
	  }
	  O = ivector(1,T1); /* alloc space for observation sequence O */
	  q = ivector(1,T1); /* alloc space for state sequence q */
	  GenSequenceArray(&hmm, seed, T1, O, q);
	  PrintSequence1(fp1, T1, O);
	  free_ivector(O, 1, T1);
	  free_ivector(q, 1, T1);
	}
	fclose(fp1);
	FreeHMM(&hmm);
	return 1;
}

void Usage(char *name)
{
  printf("Usage error \n");
  printf("Usage: %s -i <mod.hmm> -o <outputfile> -n <numsequences> -t <sequence length> [-r <randomize>] \n", name);
  printf("  t = length of sequence (will vary if random bit is set)\n");
  printf("  r =  1 - randomize 0 - not randomize\n");
  printf("  i = mod.hmm is a file with HMM parameters\n");  
  printf("  o = outputfile is the file to which the output will be written\n");
  printf("  n = numsequences is the number of sequences to be generated\n");
}
