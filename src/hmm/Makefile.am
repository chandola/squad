AM_LDFLAGS = -L$(libdir) 
INCLUDES = -I$(includedir) 
LDADD = -ltools -lm
bin_PROGRAMS = HMMGenseq HMMPredict HMMTrain
EXTRA_DIST = hmm.h nrutil.h

HMMGenseq_SOURCES = HMMGenseq.c sequence.c nrutil.c hmmrand.c hmmutils.c

HMMPredict_SOURCES = HMMPredict.c viterbi.c nrutil.c sequence.c hmmutils.c  hmmrand.c

HMMTrain_SOURCES = HMMTrain.c baum.c nrutil.c sequence.c hmmutils.c forward.c backward.c hmmrand.c

dist_doc_DATA = hmm_README