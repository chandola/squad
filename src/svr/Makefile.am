AM_LDFLAGS = -L$(libdir) 
INCLUDES = -I$(includedir) 
LDADD = -ltools
bin_PROGRAMS = SVRTrain SVRPredict
EXTRA_DIST = svm.h sax.h
SVRTrain_SOURCES = SVRTrain.cpp svm.cpp sax.cpp
SVRPredict_SOURCES = SVRPredict.cpp svm.cpp sax.cpp
dist_doc_DATA = svr_README
