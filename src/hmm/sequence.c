/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   22 February 1998 
**      File:   sequence.c
**      Purpose: Routines for generating, reading and writing sequence of
**		observation symbols.
**      Organization: University of Maryland
**	
**	Update:
**	Author: Tapas Kanungo
**	Purpose: To make calls to generic random number generators
**		and to change the seed everytime the software is executed.
**
**      $Id: sequence.c,v 1.2 1998/02/23 06:19:41 kanungo Exp kanungo $
*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"

void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q)
{
        int     t = 1;
	hmmsetseed(seed); 
 
        q[1] = GenInitalState(phmm);
        O[1] = GenSymbol(phmm, q[1]);
 
        for (t = 2; t <= T; t++) {
                q[t] = GenNextState(phmm, q[t-1]);
                O[t] = GenSymbol(phmm, q[t]);
        }
}

int GenInitalState(HMM *phmm)
{
        double val, accum;
        int i, q_t;
 
        val = hmmgetrand();
        accum = 0.0;
        q_t = phmm->N;
        for (i = 1; i <= phmm->N; i++) {
                if (val < phmm->pi[i] + accum) {
                        q_t = i;
                        break;
                }
                else {
                                accum += phmm->pi[i];
                }
        }
 
        return q_t;
}

int GenNextState(HMM *phmm, int q_t)
{
        double val, accum;
        int j, q_next;
 
        val = hmmgetrand();
        accum = 0.0;
        q_next = phmm->N;
        for (j = 1; j <= phmm->N; j++) {
                if ( val < phmm->A[q_t][j] + accum ) {
                        q_next = j;
                        break;
                }
                else
                        accum += phmm->A[q_t][j];
        }
 
        return q_next;
}
int GenSymbol(HMM *phmm, int q_t)
{
        double val, accum;
        int j, o_t;
 
        val = hmmgetrand();
        accum = 0.0;
        o_t = phmm->M;
        for (j = 1; j <= phmm->M; j++) {
                if ( val < phmm->B[q_t][j] + accum ) {
                       o_t = j;
                       break;
                }
                else
                        accum += phmm->B[q_t][j];
        }
 
        return o_t;
}
 
void ReadSequence(FILE *fp, int *pT, int **pO)
{
        int *O;
        int i;
 
        fscanf(fp, "T= %d\n", pT);
        O = ivector(1,*pT);
        for (i=1; i <= *pT; i++)
                fscanf(fp,"%d", &O[i]);
        *pO = O;
}
 
int ReadSequence1(FILE *fp, int *pT, int **pO, int maxlength)
{
  int *O,*O1;
  int i;
  char buf[2*maxlength];
  if(fgets(buf,2*maxlength,fp) == NULL) return -1;
  int len = strlen(buf);
  buf[len-1] = '\0';
  O1 = (int *) malloc(sizeof(int)*maxlength);
  char *buf1;
  buf1 = strtok(buf," ");
  int ind = 0;
  while(buf1 != NULL){
    O1[ind] = atoi(buf1);
    buf1 = strtok(NULL," ");
    ind++;
  }
  *pT = ind;
  O = ivector(1,*pT);
  for (i=1; i <= *pT; i++) O[i] = O1[i-1];
  *pO = O;
  free(O1);
  return 1;
}
 
void PrintSequence(FILE *fp, int T, int *O){
        int i;
 
        fprintf(fp, "T= %d\n", T);
        for (i=1; i <= T; i++) 
                fprintf(fp,"%d ", O[i]);
	fprintf(fp,"\n"); 
}

void PrintSequence1(FILE *fp, int T, int *O){
        int i;
        for (i=1; i <= T; i++) fprintf(fp,"%d ", O[i]);
	fprintf(fp,"\n");
}

void PrintTransProbs(FILE *fp, int T, int *O, HMM *phmm){
  int i;
  for (i=2; i < T; i++)
    fprintf(fp,"%.6f ",phmm->A[O[i-1]][O[i]]);
  fprintf(fp,"%.6f\n",phmm->A[O[T-1]][O[T]]);
}
