#
# Design for Repseek
# compile all src and output one binary
# <amikezor>
# April 2004
# 

# MACHINE = SUN
# MACHINE = SGI
# MACHINE = LINUX
# MACHINE = OSF

MACHINE = MACOSX


# ???: MALLOC=-lmalloc
# others: MALLOC =  
# macosx debug: MALLOC =  -lMallocDebug
MALLOC= 


# SGI: PROTO= PROTO=1
# others: PROTO= PROTO=0
PROTO= PROTO=0

# SGI: RANLIB= touch
# others: RANLIB= ranlib

RANLIB= ranlib 

# Users can choose:

#CC= cc
CC= gcc


CFLAGS= -O4 -Wall
#CFLAGS= -O2 
#CFLAGS= -pg
#CFLAGS= -g

INSTALLDIR = $$HOME/bin

##### defined

SHELL = bash

LDFLAGS =  -lm $(MALLOC)
#LDFLAGS =  -g -lm $(MALLOC)

SRC  = sort.c\
       help.c\
       output.c\
       filter2seq.c \
       filter.c \
       families.c \
       families_2seqs.c\
       memory.c\
       memory_align.c\
       sequence.c\
       readfst.c\
       lmin.c \
       smin.c\
       KMRK.c \
       KMRK_Seeds.c \
       KMRK_merge_seeds.c\
       KMRK_mask.c\
       align_matrix.c\
       align_seeds.c\
       align_di.c\
       align_blast2like.c\
       read_seeds.c


OBJ  = $(SRC:.c=.o)


## Rules

default:
	@echo "++ Repeats Search Engines ++"
	@echo "edit Makefile and set MACHINE, RANLIB, PROTO, CC, CFLAGS and eventually MALLOC"
	@echo " "
	@echo "To compile: make repseek"
	@echo "To clean:   make clean"


%.o: %.c
	$(CC) $(CFLAGS) -D$(MACHINE) -c -o $@ $<;


repseek:  $(OBJ) main_repseek.c
	$(CC) $(CFLAGS)  -D$(MACHINE) -o $@ $(OBJ) main_repseek.c $(LDFLAGS);


install: repseek
	cp  repseek $(INSTALLDIR)


clean:
	+rm -f *.o repseek

archive:
	\rm -rf RepSeek; mkdir RepSeek
	cp -r $(SRC) *.h Makefile main_repseek.c *.txt RepSeek;
	tmp=`date | awk '{print $$3$$2$$NF}'`;tar czvf RepSeek.$$tmp.tgz  RepSeek;
	\rm -rf RepSeek

