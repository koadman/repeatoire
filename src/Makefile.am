
if DEBUG
D_CXXFLAGS = -Wall -g -D__GNDEBUG__ 
else
endif
OPTIMIZATION = -O2 -funroll-loops -fomit-frame-pointer -ftree-vectorize -funsafe-math-optimizations
AM_CXXFLAGS = -DCOMMAND_LINE $(OPTIMIZATION) $(D_CXXFLAGS) 
AM_LDFLAGS = @STATIC_FLAG@ 

TUIUIU_H = tuiuiu/tuiglobals.h tuiuiu/index_str.h tuiuiu/display.h tuiuiu/itree.h tuiuiu/tuilist.h  tuiuiu/util.h tuiuiu/emptyBlock.h tuiuiu/mono_and_multi_commons.h
REPSEEK_H = repseek/sequence.h repseek/readfst.h repseek/help.h repseek/families.h repseek/KMRK.h repseek/KMRK_Seeds.h repseek/KMRK_merge_seeds.h repseek/KMRK_mask.h repseek/read_seeds.h repseek/filter.h repseek/align.h repseek/memory.h repseek/smin.h repseek/lmin.h repseek/sort.h repseek/memory_align.h repseek/repseek_types.h repseek/output.h

TUIUIU_SRC = \
tuiuiu/tuiglobals.c tuiuiu/index_str.c tuiuiu/display.c  tuiuiu/itree.c tuiuiu/tuilist.c  tuiuiu/mono_and_multi_commons.c tuiuiu/util_excelent_mono.c

REPSEEK_SRC = \
repseek/KMRK.c repseek/KMRK_Seeds.c repseek/KMRK_merge_seeds.c repseek/KMRK_mask.c repseek/align_blast2like.c repseek/align_di.c repseek/align_matrix.c repseek/align_seeds.c repseek/families.c repseek/families_2seqs.c repseek/filter.c repseek/filter2seq.c repseek/help.c repseek/lmin.c repseek/memory.c repseek/memory_align.c repseek/output.c repseek/read_seeds.c repseek/readfst.c repseek/sequence.c repseek/smin.c repseek/sort.c

bin_PROGRAMS = repeatoire 
repeatoire_SOURCES = \
repeatoire.cpp $(TUIUIU_SRC) $(REPSEEK_SRC) 

LIBRARY_CL = $(DEPS_LIBS) 

repeatoire_LDADD = $(LIBRARY_CL)

INCLUDES = @DEPS_CFLAGS@ 


