LDADD = -ltools -lfftw3
AM_CPPFLAGS = -I$(includedir) 
bin_PROGRAMS = KNNC findsimc KNNC_TOP

KNNC_SOURCES = KNNC.cpp
KNNC_LDFLAGS = -L$(libdir) 

KNNC_TOP_SOURCES = KNNC_TOP.cpp
KNNC_TOP_LDFLAGS = -L$(libdir) 

findsimc_SOURCES = findsimc.cpp
findsimc_LDFLAGS = -L$(libdir) 

dist_doc_DATA = knnc_README