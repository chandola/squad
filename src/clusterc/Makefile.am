LDADD = -ltools -lfftw3
AM_CPPFLAGS = -I$(includedir) 

bin_PROGRAMS = CLUSTERCTrain CLUSTERCPredict
EXTRA_DIST = CLUSTERCTrain.h

CLUSTERCTrain_SOURCES = CLUSTERCTrain.cpp
CLUSTERCTrain_LDFLAGS = -L$(libdir) 

CLUSTERCPredict_SOURCES = CLUSTERCPredict.cpp
CLUSTERCPredict_LDFLAGS = -L$(libdir) 

dist_doc_DATA = clusterc_README