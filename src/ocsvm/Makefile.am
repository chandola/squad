AM_LDFLAGS = -L$(libdir) 
INCLUDES = -I$(includedir) 
LDADD = -ltools
bin_PROGRAMS = OCSVMTrain OCSVMPredict
EXTRA_DIST = svm.h sax.h
OCSVMTrain_SOURCES = OCSVMTrain.cpp svm.cpp sax.cpp
OCSVMPredict_SOURCES = OCSVMPredict.cpp svm.cpp sax.cpp
dist_doc_DATA = ocsvm_README
