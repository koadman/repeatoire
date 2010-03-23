
if DEBUG
D_CXXFLAGS = -Wall -g -D__GNDEBUG__ 
else
endif
OPTIMIZATION = -O3 -Wall -funroll-loops -fomit-frame-pointer -ftree-vectorize -funsafe-math-optimizations
AM_CXXFLAGS = -DCOMMAND_LINE $(OPTIMIZATION) $(D_CXXFLAGS) 
AM_LDFLAGS = @STATIC_FLAG@ 

bin_PROGRAMS = repeatoire 
repeatoire_SOURCES = \
repeatoire.cpp

LIBRARY_CL = $(DEPS_LIBS) 

repeatoire_LDADD = $(LIBRARY_CL)

INCLUDES = @DEPS_CFLAGS@ 


