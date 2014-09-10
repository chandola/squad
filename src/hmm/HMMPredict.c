/*
**      Author: Varun Chandola, chandola@cs.umn.edu
**      Date:   06 January 2009
**      File:   HMMPredict.c
**      Purpose: driver for testing the Viterbi code.
**      Organization: University of Minnesota
*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"

/* Modify the MAXLENGTH to handle sequences longer than 2000.
   Make sure to compile with -D__LONG flag.
   Default compilation has maximum length set to 500.
*/
#ifdef __LONG
#define MAXLENGTH 2000
#else
#define MAXLENGTH 500
#endif

int main (int argc, char **argv)
{
	int 	T; 
	HMM  	hmm;
	int	*O;	/* observation sequence O[1..T] */
	int	*q;	/* state sequence q[1..T] */
	double **delta;
	int	**psi;
	double 	proba; 
	FILE	*fp,*fp1;

	if (argc != 4) {
		printf("Usage error \n");
		printf("Usage: testvitall <model.hmm> <test.seq> <outputfile>\n");
		exit (1);
	}
	
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found\n", argv[1]);
		exit (1);
	}
	ReadHMM(fp, &hmm);
	fclose(fp);

	fp = fopen(argv[2], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found\n", argv[2]);
		exit (1);
	}

	fp1 = fopen(argv[3], "w");
	if (fp1 == NULL) {
		fprintf(stderr, "Error: File %s could not be opened\n", argv[3]);
		exit (1);
	}

	while(1){
	  if(ReadSequence1(fp, &T, &O, MAXLENGTH) == -1) break;
	  q = ivector(1,T);
	  delta = dmatrix(1, T, 1, hmm.N);
	  psi = imatrix(1, T, 1, hmm.N);
	  Viterbi(&hmm, T, O, delta, psi, q, &proba); 
	  //fprintf(stdout, "Viterbi  MLE log prob = %E\n", log(proba));
	  //fprintf(stdout, "Optimal state sequence:\n");
	  //	  PrintSequence1(stdout, T, q);
	  PrintTransProbs(fp1,T,q,&hmm);
	  //printf("------------------------------------\n");
	  free_ivector(q, 1, T);
	  free_ivector(O, 1, T);
	  free_imatrix(psi, 1, T, 1, hmm.N);
	  free_dmatrix(delta, 1, T, 1, hmm.N);
	}
	fclose(fp);fclose(fp1);
	FreeHMM(&hmm);
	return 0;
}
