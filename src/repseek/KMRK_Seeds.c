/*
	Copyright (C) 2006 G achaz, F boyer, E rocha, A viari and E coissac.

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public License
	as published by the Free Software Foundation; either version 2.1
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 
 	for more information, please contact guillaume achaz <achaz@abi.snv.jussieu.fr>
*/
#include "KMRK_Seeds.h"
#include "memory.h"
#include <stdlib.h>
#include <string.h>
#include "sequence.h"


static void SetMultipleLenInvSeeds(
                   SmallSeed_type* seeds,           /* an array of a structure with pos1 and pos2 */
				   int32_t         nseeds,          /* number of small seeds */
				   int32_t         wordSize,        /* lmin */
				   int8_t          same,            /* 1 if one sequence, 0 otherwise */
				   AllSeeds_type *PtrAllSeeds);     /* the final seeds with length */


/*
	Concatenate
	DirSeq\@InvSeq\0
*/
static char* makeDirInvSeq(char* seq, int32_t size)
{
  char *SeqInv;

	seq = (char *)MyRealloc( (void *)seq, (size*2+2)*sizeof(char),
                            (size+1)* sizeof(char),  "makeDirInvSeq: Cannot allocate space for reverse sequence");
	
	SeqInv= seq + size + 1;
	seq[size]=0;
	invseq(seq, SeqInv);
	seq[size]='@';
	SeqInv[size]=0;

  return seq;
}

/*
	Merge the seq1 with seq2
*/
static char *merge2seq(char* seq1, char* seq2,
		       int32_t size1, int32_t size2)
{
  char * dest;

  seq1 = (char *)MyRealloc((void *)seq1,  (size1+size2+2) *sizeof(char),
           (size1+1)*sizeof(char), "merg2seq: Cannot allocate space for reverse sequence");

  dest = seq1 + size1 + 1;
  seq1[size1]='@';
  memcpy(dest,seq2,size2);
  dest[size2]=0;

  return seq1;
}

static int32_t dirDelta(SmallSeed_type* seed)
{
    return seed->pos2 - seed->pos1;
}

void KMRK_SetMultipleLenDirSeeds(SmallSeed_type* seeds,
                              int32_t         nseeds,
                              int32_t         wordSize,
                              AllSeeds_type   *PtrAllSeeds)
{
    
    int32_t i,j;                                /* dummy counters j=kept seeds ; i = current seed */
    int32_t curLen=wordSize;                    /* Length of the current seed */  
    int32_t delta;
    SmallSeed_type *mainSeed;
    SmallSeed_type *curSeed;
    int32_t add;

/*    fprintf(stderr,"New Version\n");*/

    KMRK_sortSeeds(seeds,nseeds,KMRK_cmpDeltaSeedsPos);

    for(j=0,mainSeed=seeds ; 
        j < nseeds; 
        j=i, mainSeed=curSeed)
    {    
        /* foreach seed */
        delta = dirDelta(mainSeed);
        curLen=wordSize;
        
        for (i=j+1,curSeed=mainSeed+1; 
             i < nseeds && 
             dirDelta(curSeed)==delta &&
             (curSeed->pos1 - mainSeed->pos1) <= curLen;
             i++,curSeed++)
        {
            add=wordSize - mainSeed->pos1 - curLen + curSeed->pos1;
            if (add < 0) add = 0;
            curLen+=add;
        }   
        
        KMRK_pushSeed(PtrAllSeeds,
                      seeds[j].pos1,seeds[j].pos2,
                      curLen,
                      1);
/*	
	if(seeds[j].pos1 == seeds[j].pos2)
		fprintf(stderr, "KMRK_SetMultipleLenDirSeeds: %d %d %d, bye\n", seeds[j].pos1, seeds[j].pos2, curLen), exit(0);
*/	
	
    }
}

static int32_t invDelta(SmallSeed_type* seed)
{
    return seed->pos2 + seed->pos1;
}

static void SetMultipleLenInvSeeds(
                   SmallSeed_type* seeds,           /* an array of a structure with pos1 and pos2 */
				   int32_t         nseeds,		    /* number of small seeds */
				   int32_t         wordSize,	    /* lmin */
				   int8_t          same,		    /* 1 if one sequence, 0 otherwise */
				   AllSeeds_type *PtrAllSeeds)	    /* the final seeds with length */
{

    int32_t i,j;                                /* dummy counters j=kept seeds ; i = current seed */
    int32_t curLen=wordSize;                    /* Length of the current seed */  
    int32_t delta;                              /* the pos1+pos2 - constant if seeds should be merged */
    int32_t pos2;                               /* after merge, new pos2 */
    SmallSeed_type *mainSeed;                   /* the first one of a series */
    SmallSeed_type *curSeed;                    /* the current one */
    /*int32_t add;                                 what is added to the Lenght for the merged form */



	KMRK_sortSeeds(seeds,nseeds,KMRK_cmpDeltaInvSeedsPos);   /* qsort() seeds by delta, then pos1 */
 	

     for(j=0,mainSeed=seeds ;  j < nseeds;   j=i, mainSeed=curSeed)
    {
        
        delta = invDelta(mainSeed);        /* pos1+pos2 */
        curLen=wordSize;                   /* at start curLen is equal to lmin */
        
        for (i=j+1, curSeed=mainSeed+1; 
             i < nseeds &&
			 invDelta(curSeed)==delta &&                         /* Constant if they should merge */
			 (curSeed->pos1 - (curSeed-1)->pos1) <= wordSize;    /* if they overlap */
             i++, curSeed++)
		{;}
		
		
		curLen += (  (curSeed-1)->pos1 - mainSeed->pos1);

        if ( same && ( seeds[j].pos1+curLen-1 > seeds[i-1].pos2) ){
           pos2 = seeds[j].pos1;
           curLen += seeds[i-1].pos2 - seeds[j].pos1;
        }
		else
		   pos2 = seeds[i-1].pos2;
	

 /*
 		Old Code - there was a bug for palindrome
 
        for (i=j+1, curSeed=mainSeed+1; 
             i < nseeds && invDelta(curSeed)==delta && (curSeed->pos1 - mainSeed->pos1) <= curLen;
             i++, curSeed++)
		{
        {
			add = wordSize-curLen - (mainSeed->pos1 - curSeed->pos1);
            if (add < 0) add = 0;             
            curLen+=add;
        }
			
        if ( same && (seeds[j].pos1+curLen>= seeds[i-1].pos2) )
        {
            curLen = (seeds[i-1].pos2 - seeds[j].pos1 + 1) * 2;
            pos2 = seeds[j].pos1;
        }
        else
             pos2 = seeds[i-1].pos2;
*/
		
		
        KMRK_pushSeed(PtrAllSeeds,
                      seeds[j].pos1,pos2,
                      curLen,
                      0);


    }
}

AllSeeds_type *KMRK_allocSeeds(AllSeeds_type *AllSeeds, 
			       int32_t size, 
			       int8_t opt_dir, 
			       int8_t opt_inv)
{

	AllSeeds_type *reponse;

	if(opt_inv != 1 && opt_dir != 1)
	{
		fprintf(stderr,"AllocSeeds: requiere at least "
		               "one of opt_dir and opt_inv to be 1\n");
		exit(4);
	}

	reponse = AllSeeds;

	if (!reponse)
		reponse = MyCalloc( 1, sizeof(AllSeeds_type),"KMRK_allocSeeds: cannot allocate new data structure");

	if(opt_dir)
	{
		if (reponse->dirSeeds==NULL)
			reponse->dirSeeds = (Seed_type *)MyCalloc( size,sizeof(Seed_type),"KMRK_allocSeeds: cannot allocate Direct Seeds");
		else
		{
			if(size)
				reponse->dirSeeds = (Seed_type *)MyRealloc( (void *)reponse->dirSeeds,
				                     size*sizeof(Seed_type),  reponse->cDirSeeds*sizeof(Seed_type),
											"allocSeeds: cannot reallocate Direct Seeds" );
			else
	  		{
				MyFree( reponse->dirSeeds, reponse->cDirSeeds*sizeof(Seed_type) );
	    			reponse->dirSeeds=NULL;
	    			reponse->cDirSeeds=0;
			}
		}

		reponse->cDirSeeds=size;
	}
	
	
	if(opt_inv)
	{
		if (reponse->invSeeds==NULL)
			reponse->invSeeds = (Seed_type *)MyCalloc( size, sizeof(Seed_type),"KMRK_allocSeeds: cannot allocate Inverted Seeds");
		else
		{
			if(size)
				reponse->invSeeds = (Seed_type *)MyRealloc( (void *)reponse->invSeeds,
				                     size*sizeof(Seed_type),  reponse->cInvSeeds*sizeof(Seed_type),
											"allocSeeds: cannot reallocate Inverted Seeds" );
			else
	  		{
				MyFree( reponse->invSeeds, reponse->cInvSeeds*sizeof(Seed_type) );
	    			reponse->invSeeds=NULL;
	    			reponse->cInvSeeds=0;
			}
		}

		reponse->cInvSeeds=size;
	}
	
      
  return reponse;
}

void KMRK_freeSeeds(AllSeeds_type *AllSeeds)
{
  if (!AllSeeds)
    return;

  if (AllSeeds->dirSeeds)
    MyFree(AllSeeds->dirSeeds, AllSeeds->cDirSeeds*sizeof(Seed_type) );
	AllSeeds->dirSeeds=NULL;

  if (AllSeeds->invSeeds)
    MyFree(AllSeeds->invSeeds, AllSeeds->cInvSeeds*sizeof(Seed_type) );
	AllSeeds->dirSeeds=NULL;

  MyFree(AllSeeds, 1*sizeof( AllSeeds_type ) );
}

void KMRK_compactSeeds(AllSeeds_type *AllSeeds)
{
  if (AllSeeds)
    {
      if (AllSeeds->dirSeeds)
        KMRK_allocSeeds(AllSeeds, AllSeeds->nDirSeeds, 1, 0);

      if (AllSeeds->invSeeds)
         KMRK_allocSeeds(AllSeeds, AllSeeds->nInvSeeds, 0, 1);

    }
}


void KMRK_pushSeed(AllSeeds_type *AllSeeds,
		   int32_t       pos1,
		   int32_t       pos2,
		   int32_t       length,
		   int8_t        dir)
{
  Seed_type* stack;
  int32_t    maxcount;
  int32_t    index;

  if (dir)
    {
      dir      = 1;
      stack    = AllSeeds->dirSeeds;
      maxcount = AllSeeds->cDirSeeds;
      index    = AllSeeds->nDirSeeds;
    }
  else
    {
      dir      = 0;
      stack    = AllSeeds->invSeeds;
      maxcount = AllSeeds->cInvSeeds;
      index    = AllSeeds->nInvSeeds;
    }

  if (index == maxcount)
    {
      (void) KMRK_allocSeeds(AllSeeds,
			     maxcount * 2,
			     dir,
			     !dir);
      
      if (dir)
	stack    = AllSeeds->dirSeeds;
      else
	stack    = AllSeeds->invSeeds;
    }

  stack+=index;

  stack->pos1   = pos1;
  stack->pos2   = pos2;
  stack->length = length;
 
  if (dir)
    AllSeeds->nDirSeeds++;
  else
    AllSeeds->nInvSeeds++;

}




AllSeeds_type* KMRK_enumerateDirectCouple(AllSeeds_type* Seeds,
					  int32_t expected,
					  int32_t wordSize,
					  vn_type* stack,
					  int32_t seq)
{

  int32_t xmin;
  int32_t xmax;
  int32_t i;
  int32_t j;
  int32_t k;
  int32_t next;
  int32_t* n;
  int32_t nseeds;

  SmallSeed_type *list;

  list = (SmallSeed_type *)MyCalloc( expected, sizeof(SmallSeed_type) ,
	      "KMRK_enumerateDirectCouple: cannot allocate DirectCouple memory");

  nseeds = 0;

  n    = GETVECTOR(stack,n);
  
  if (seq)
    xmin = stack->limits[seq-1];
  else
    xmin = 0;

  xmax = stack->limits[seq];

  for (i=1; i <= xmax; i++)
	  if (IS_MARKED(n,i))     /* Check begining of chained list */
      {
				/* Look for begining of sequence of interest */   
		  for( ;(i <= xmin) && (i != GET(n,i));
			   i=GET(n,i));
		  
				/* for each factor in  sequence of interest */   
		  for (j=i; 
			   (j != GET(n,j)) && (j <= xmax);
			   j = GET(n,j))
		  {
			  next = GET(n,j);
			  if (next <= xmax)
				  do
				  {
					  k = next;
					  next = GET(n,k);
					  list[nseeds].pos1 = j-1;
					  list[nseeds].pos2 = k-1;
					  nseeds++;
				  } while ((k!=next) && (next <= xmax));
		  }
      } ;
  
  Seeds = KMRK_allocSeeds(Seeds,
			  expected/20+1,
			  1,0);

/*  fprintf(stderr,"Expected direct couple : %d\n",expected);*/
  
  KMRK_SetMultipleLenDirSeeds(list,nseeds,wordSize,Seeds);
  MyFree(list, expected*sizeof(SmallSeed_type) );
  KMRK_compactSeeds(Seeds);
	  
  return Seeds;
}


/*
	From KMR-K Stacks to SmallSeeds
*/
AllSeeds_type* KMRK_enumerateInvertedCouple(
                        AllSeeds_type* Seeds,
					    int32_t expected,              /* the expected number of couple */
					    int32_t wordSize,              /* the Lmin */
					    vn_type* stack)
{
  int32_t xmax;
  int32_t invmax;
  int32_t posinv;
  int32_t i;
  int32_t j;
  int32_t k;
  int32_t memk;
  int32_t* n;
  int32_t next;

  int32_t nseeds;        /* number of seeds */

  SmallSeed_type *list;   /* seed list with only pos1 and pos2 --simple output from kmrk */


  list = (SmallSeed_type *)MyCalloc( expected, sizeof(SmallSeed_type) ,
	      "KMRK_enumerateInvertedCouple: cannot allocate InvertedCouple memory");

  nseeds = 0;

  n    = GETVECTOR(stack,n);
  
  xmax = stack->limits[0];
  invmax = stack->limits[1];

  for (i=1; i <= xmax; i++)
	  if (IS_MARKED(n,i))
      {
		  for(memk=i ;
			  (memk < xmax) && 
			  memk != GET(n,memk) &&
			  (memk <= invmax);
			  memk=GET(n,memk));
		  
		  if ((memk > xmax) && (memk <= invmax))
			  for (j=i;
				   (j <= xmax) && (j != GET(n,j));
				   j=GET(n,j))
			  {
				  next = memk;
				  do
				  {
					  k = next;
					  next = GET(n,k);
					  posinv = 2 * xmax - k -wordSize + 3;
					  if (j <= posinv)
					  {
						  list[nseeds].pos1=j-1;
						  list[nseeds].pos2=posinv-1;
						  nseeds++;
					  }
				  } while ((next <= invmax) && (k != next));
			  }
      }
		  
		  Seeds = KMRK_allocSeeds(Seeds, expected/20+1, 0,1);

/*  fprintf(stderr,"Expected inverted couple : %d\n",expected);*/

  SetMultipleLenInvSeeds(list,nseeds,wordSize,1,Seeds);   /* turn the Small seeds into merged seeds */
  MyFree(list, expected*sizeof(SmallSeed_type) );
  KMRK_compactSeeds(Seeds);


  return Seeds;
}

AllSeeds_type* KMRK_enumerateInterCouple(AllSeeds_type* Seeds,
					 int32_t seq1,
					 int32_t seq2,
					 int32_t expected,
					 int32_t wordSize,
					 vn_type* stack)
{
  int32_t xmin;
  int32_t xmax;
  int32_t ymax;
  int32_t ymin;
  int32_t pos1;
  int32_t pos2;
  int32_t i;
  int32_t j;
  int32_t k;
  int32_t memj;
  int32_t memk;
  int32_t* n;
  int32_t next;

  int32_t nseeds;


  SmallSeed_type *list;
  nseeds=0;

  list = (SmallSeed_type *)MyCalloc( expected, sizeof(SmallSeed_type) ,
	      "KMRK_enumerateInterCouple: cannot allocate InterCouple memory");

  n    = GETVECTOR(stack,n);
  
  if (seq1==0)
    xmin=0;
  else
    xmin = stack->limits[seq1-1];

  xmax = stack->limits[seq1];
  ymin = stack->limits[seq2-1];
  ymax = stack->limits[seq2];

  for (i=1; i <= xmax; i++)
	  if (IS_MARKED(n,i))
      {
		  for(memj=i ;
			  (memj < xmin) && 
			  memj != GET(n,memj);
			  memj=GET(n,memj));
		  
		  if ((memj > xmin) && (memj <= xmax))
		  {
			  for(memk=memj ;
				  (memk < ymin) && 
				  memk != GET(n,memk);
				  memk=GET(n,memk));
			  
			  if ((memk > ymin) && (memk <= ymax))
				  for (j=memj;
					   (j <= xmax) && (j != GET(n,j));
					   j=GET(n,j))
				  {
					  next = memk;
					  do
					  {
						  k = next;
						  next = GET(n,k);
						  if (seq1 > 0)
							  pos1 = j - xmin - 1;
						  else
							  pos1 = j;
						  pos2 = k - ymin - 1;
						  list[nseeds].pos1 = pos1 - 1;
						  list[nseeds].pos2 = pos2 - 1;
						  nseeds++;
					  } while ((next <= ymax) && (k != next));
				  }
		  }
      }

	  
	  Seeds = KMRK_allocSeeds(Seeds,
							  expected/20+1,
							  1,0);
	  
	/*  fprintf(stderr,"Expected inter-direct couple : %d\n",expected);*/
	  KMRK_SetMultipleLenDirSeeds(list,nseeds,wordSize,Seeds);
	  
	  MyFree(list, expected*sizeof(SmallSeed_type) );
	  
	  KMRK_compactSeeds(Seeds);
	  
	  return Seeds;
}

AllSeeds_type* KMRK_enumerateInterInvertedCouple(AllSeeds_type* Seeds,
						 int32_t seq2,
						 int32_t expected,
						 int32_t wordSize,
						 vn_type* stack)

{
  int32_t xmin;
  int32_t xmax;
  int32_t ymax;
  int32_t ymin;
  int32_t posinv;
  int32_t pos2;
  int32_t i;
  int32_t j;
  int32_t k;
  int32_t memj;
  int32_t memk;
  int32_t* n;
  int32_t next;

  int32_t nseeds;

  SmallSeed_type *list;

  list = (SmallSeed_type *)MyCalloc( expected, sizeof(SmallSeed_type) ,
	      "KMRK_enumerateInterCouple: cannot allocate InterCouple memory");


  nseeds = 0;
  n    = GETVECTOR(stack,n);

  if (seq2 < 2)
    {
      fprintf(stderr,"enumerateInterInvertedCouple: seq2 must be differente to 0");
      exit(4);
    }
    
  
  xmin = stack->limits[0];
  xmax = stack->limits[1];

  ymin = stack->limits[seq2-1];
  ymax = stack->limits[seq2];

  Seeds = KMRK_allocSeeds(Seeds,
			  expected,
			  0,1);

  for (i=1; i <= xmax; i++)
	  if (IS_MARKED(n,i))
      {
		  for(memj=i ;
			  (memj < xmin) && 
			  memj != GET(n,memj);
			  memj=GET(n,memj));
		  
		  if ((memj > xmin) && (memj <= xmax))
		  {
			  for(memk=memj ;
				  (memk < ymin) && 
				  memk != GET(n,memk);
				  memk=GET(n,memk));
			  
			  if ((memk > ymin) && (memk <= ymax))
				  for (j=memj;
					   (j <= xmax) && (j != GET(n,j));
					   j=GET(n,j))
				  {
					  next = memk;
					  do
					  {
						  k = next;
						  next = GET(n,k);
						  posinv = 2 * xmin - j -wordSize + 3;
						  pos2 = k - ymin - 1;
						  list[nseeds].pos1=posinv-1;
						  list[nseeds].pos2=pos2-1;
						  nseeds++;
					  } while ((next <= ymax) && (k != next));
				  }
		  }
      }

		  Seeds = KMRK_allocSeeds(Seeds,
			  expected/20+1,
			  0,1);

/*  fprintf(stderr,"Expected inter-inverted couple : %d\n",expected);*/
  SetMultipleLenInvSeeds(list,nseeds,wordSize,0,Seeds);
  MyFree(list, expected* sizeof(SmallSeed_type) );
  KMRK_compactSeeds(Seeds);

		  
  return Seeds;
}

  
int32_t KMRK_cmpSeedsPos(SmallSeed_type *s1, SmallSeed_type *s2)
{
  if (s1->pos1==s2->pos1)
    return s1->pos2 - s2->pos2;
  else
    return s1->pos1 - s2->pos1;
}

int32_t KMRK_cmpDeltaSeedsPos(SmallSeed_type *s1, SmallSeed_type *s2)
{
    int32_t delta1 =  s1->pos2-s1->pos1;
    int32_t delta2 =  s2->pos2-s2->pos1;
    
    if (delta1==delta2)
        return s1->pos1 - s2->pos1;
    else
        return delta1 - delta2;
}

/*
	Sort by deltas, then pos1
*/
int32_t KMRK_cmpDeltaInvSeedsPos(SmallSeed_type *s1, SmallSeed_type *s2)
{
    int32_t delta1 =  s1->pos2+s1->pos1;
    int32_t delta2 =  s2->pos2+s2->pos1;
    
    if (delta1==delta2)
        return s1->pos1 - s2->pos1;        /* then sort by pos1 */
    else
        return delta1 - delta2;            /* sort by deltas */
}


void KMRK_sortSeeds(SmallSeed_type* seeds,
		    int32_t         nseeds,
		    KMRK_SORT_SEEDS_FUNC_PTR(compare))
{

  qsort(seeds, 
	nseeds, 
	sizeof(SmallSeed_type),
	  (int (*)(const void *, const void *))compare);
}

AllSeeds_type* KMRK_get_seeds(char **seq,
                              int32_t SimpleSeqLen,
                              int16_t Lmin,
                              int8_t opt_dir,
                              int8_t opt_inv, 
                              int8_t opt_verbose,
                              masked_area_table_t *mask) 
{
  AllSeeds_type* AllSeeds;
  char *ConcatSeq = *seq;
  vn_type * stacks;
  int32_t dirExpect=0;
  int32_t invExpect=0;
  KMRK_QUORUM_FUNC_PTR(quorum);
  
  if(opt_inv != 1 && opt_dir != 1)
      fprintf(stderr, "get_seeds: requiered at least opt_dir or opt_inv to be 1\n"), exit(4);

		
  if(opt_inv)
      ConcatSeq = makeDirInvSeq(ConcatSeq,SimpleSeqLen);   /* create a sequence with "DirSeq\@InvSeq\0" */

	if (opt_inv){                                     /* Are we interested in dir, inv or both ? */
		if (opt_dir)
			quorum = KMRK_DirInvQuorum;
		else
			quorum = KMRK_InvQuorum;
	}
	else
		quorum = KMRK_CoupleQuorum;

  stacks = KMRK_RunTo(ConcatSeq,             /* concatenation of both sequence */
		      Lmin,                  /* minimum length of the repeated elements = word size to construct */
		      opt_inv,               /* shall we search for inverted repeats as well ?*/
		      2,                     /* the count -- a mysterious value passed to the quorum function (unused ??) */
		      1,                     /* how many functions */
		      quorum,                /* a pointer to the selected quorum function */
		      KMRK_HashOneSequence,  /* a pointer to the init function */
              mask);

  invExpect =0;

  KMRK_markStart(stacks);

  if (opt_inv)
    {
      ConcatSeq = (char *)MyRealloc( (void *)ConcatSeq,  (SimpleSeqLen+1)*sizeof(char),       
		                            (2*SimpleSeqLen+2)*sizeof(char) , "KRMK_get_seeds: Cannot shrink memory");   /* reset mem to a sigle sequence */
      ConcatSeq[SimpleSeqLen]=0;
    }

    if(opt_inv)
        invExpect = KMRK_upperInvertedCount(stacks,Lmin);

    if(opt_dir)
        dirExpect = KMRK_upperCoupleCount(stacks);
 

  AllSeeds = NULL;

	
  MyFree(stacks->v, (stacks->vectorsize) * sizeof(int32_t));
  stacks->v=NULL;


  if (opt_dir)
    AllSeeds = KMRK_enumerateDirectCouple(AllSeeds,dirExpect,Lmin ,stacks,0);

  if (opt_inv)
    AllSeeds = KMRK_enumerateInvertedCouple(AllSeeds,invExpect,Lmin,stacks);

  KMRK_FreeVN(stacks);

  *seq = ConcatSeq;
				 
  return AllSeeds;
}  

AllSeeds_type* KMRK_get_seeds_2seqs(char **seq1,
                                    char **seq2,
                                    int32_t size1,
                                    int32_t size2,
                                    int16_t Lmin,
                                    int8_t opt_dir,
                                    int8_t opt_inv, 
                                    int8_t opt_verbose,
                                    masked_area_table_t *mask) 
{

  AllSeeds_type* AllSeeds;
  char *sequence1 = *seq1;
  char *sequence2 = *seq2;
  vn_type * stacks;
  int32_t dirExpect=0;
  int32_t invExpect=0;
  KMRK_QUORUM_FUNC_PTR(quorum);
  int32_t sizef;                  /* =size1 if only dir OR =2*size1+size2 if inv */

  if(opt_inv != 1 && 
     opt_dir != 1)
    {
      fprintf(stderr,
	      "get_seeds_2seqs: requiered at least "
	      "opt_dir or opt_inv to be 1\n");
      exit(4);
    }

  sizef = size1;

  if(opt_inv)
    {
      sequence1 = makeDirInvSeq(sequence1,size1);
      sizef+=(1+size1);
    }

  sequence1 = merge2seq(sequence1,sequence2,sizef,size2);

  if (opt_inv)
    if (opt_dir)
      quorum = KMRK_IntDirInv12Quorum;
    else
      quorum = KMRK_IntInv12Quorum;
  else
    quorum = KMRK_Int12Quorum;

  stacks = KMRK_RunTo(sequence1,
		      Lmin,
		      opt_inv,
		      2,
		      2,
		      quorum,
		      KMRK_HashOneSequence,
              mask);

  KMRK_markStart(stacks);

  sequence1= (char *)MyRealloc( (void *)sequence1, (size1+1)*sizeof(char),
	      (sizef+size2+2)*sizeof(char), "KMRK_get_seeds_2seqs: shrink memory back to size1");

  sequence1[size1]=0;

  if (opt_dir){
    if (opt_inv)
      dirExpect = KMRK_upperInterCount(stacks,0,2,Lmin);
    else
      dirExpect = KMRK_upperInterCount(stacks,0,1,Lmin);
	}
  if (opt_inv)
    invExpect = KMRK_upperInterCount(stacks,1,2,Lmin);

  AllSeeds = NULL;
  MyFree(stacks->v, stacks->vectorsize*sizeof(int32_t) );    /*  free vector v from the stacks */
  stacks->v=NULL;

  if (opt_dir){
    if (opt_inv)
      AllSeeds = KMRK_enumerateInterCouple(AllSeeds,
					   0,2,
					   dirExpect,
					   Lmin ,
					   stacks);
    else
      AllSeeds = KMRK_enumerateInterCouple(AllSeeds,
					   0,1,
					   dirExpect,
					   Lmin ,
					   stacks);

	}
  if (opt_inv)
    AllSeeds = KMRK_enumerateInterInvertedCouple(AllSeeds,
						 2,
						 invExpect,
						 Lmin ,
						 stacks);


  KMRK_FreeVN(stacks);

  *seq1 = sequence1;
				 
  return AllSeeds;
}

#define SIGN(x) (((x)<0) ? -1:(((x)>0) ? 1:0))

static int32_t compareSeedsByPos(Seed_type* s1,Seed_type* s2)
{
  if (s1->pos1 == s2->pos1)
    return SIGN(s1->pos2 - s2->pos2);
  else
    return SIGN(s1->pos1 - s2->pos1);
}

void KMRK_sortSeedsByPos(Seed_type* seeds, int32_t count)
{
  qsort(seeds,
	count,
	sizeof(Seed_type),
	(int (*)(const void *, const void *))compareSeedsByPos);
};

#undef SIGN

