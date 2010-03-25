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
#include "KMRK.h"
#include "repseek_types.h"
#include "memory.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


/********************************************
 ********************************************
 **
 ** Declaration of static functions
 **
 ********************************************
 ********************************************/


/** 
 * @brief Allocate a couple of vector v and n of same size.
 * Never return if memory space is not available.
 * 
 * @param size size of the new allocated vector couple
 * @param count number of sequence encodedable in the VN
 *              structure allocated 
 * 
 * @return a pointer to a vn_type structure
 */

static vn_type* new_vn(int32_t size,int32_t count);





/** 
 * Set to 0 all elements of the n vector of a
 * vn structure.
 * 
 * @param x a pointer to a vn_type structure
 */

static void clearNVector(vn_type* x);

  
/** 
 * Set to 0 all elements of the v vector of a
 * vn structure.
 * 
 * @param x a pointer to a vn_type structure
 */

static void clearVVector(vn_type* x);


/* static int32_t *copyNVector(vn_type *x, int32_t *dest); */


/** 
 * @brief Encode a nucleotide in a numeric format.
 * 
 * @param nucleic a char value containing A,C,G or T. 
 *                All others values induce 0 as result.
 * 
 * @return an interger code corresponding to the nucleotide
 */

static int32_t code(char nucleic);


/** 
 * @brief Construct linked list of similar symboles present in
 * v part of x.
 *
 * Equivalent to filling of P stack in standard implementation
 *
 * Result is stored in n part of x.
 * 
 * @param x a pointer to a vn_type structure
 */

static void phase1(vn_type* x);


/* static void unroll1(vn_type* x); */


static void phase2(vn_type* x, int32_t k);


static void unroll2(vn_type* x);


static void unstack(vn_type* x, int32_t k);

static vn_type *hashSequence(char *sequence, int32_t lmax, int32_t *k,
                             masked_area_table_t *mask);

static void clearNMark(vn_type *x);

static vn_type* encodeSequence(char *sequence,
                               masked_area_table_t *mask);

static void applyQuorum(vn_type *x, 
			int32_t count, 
			int32_t countSeq,
			KMRK_QUORUM_FUNC_PTR(quorum));

/********************************************
 ********************************************
 **
 ** Defintion of static functions
 **
 ********************************************
 ********************************************/



static int32_t code(char nucleic)
{

	if(nucleic == 'A')
		return 1;
	if(nucleic == 'C')
		return 2;
	if(nucleic == 'G')
		return 3;
	if(nucleic == 'T')
		return 4;
	
	return 0;
}



static void clearNVector(vn_type* x)
{
  memset(x->n,0,sizeof(int32_t)*x->size);
};

static void clearVVector(vn_type* x)
{
  memset(x->v,0,sizeof(int32_t)*x->size);
};

static void clearNMark(vn_type *x)
{
  int32_t i;
  int32_t *n;
  int32_t size;

  n = x->n;
  size = x->size;

  for (i=0; i < size; i++)
    if (n[i] < 0) n[i] = -n[i];

}

static vn_type* new_vn(int32_t size,int32_t count)
{
	vn_type* reponse;

	reponse = MyCalloc( 1, sizeof(vn_type)+sizeof(int32_t) * (count -1) , "vn_type: first calloc" );
    
	reponse -> v = (int32_t *)MyCalloc(size, sizeof(int32_t), "new_vn: I cannot not allocate a new VN structure");
	reponse->n = (int32_t *)MyCalloc(size, sizeof(int32_t),   "new_vn: I cannot not allocate a new VN structure");


  reponse->size = size;
  reponse->vectorsize = size;
  reponse->seqCount= count;

  return reponse;
}


/*
	Fill-up the n vector with a "chained list" of positions
*/
static void phase1(vn_type* x)
{

  int32_t *v;
  int32_t *n;
  int32_t size;
  int32_t i;
  int32_t symbole;
  int32_t oldlast;

  v    = GETVECTOR(x,v);
  n    = GETVECTOR(x,n);
  size = x->size;

  clearNVector(x);

  for(i = 1; i <= size; i++)
    {
      symbole = GET(v,i);
      if (symbole == i)

	SET(n,i,i);

      else if (symbole)
	{
	  oldlast = GET(n,symbole);
	  SET(n,symbole,i);
	  SET(n,i,GET(n,oldlast));
	  SET(n,oldlast,i);
	}
    }

}



static void phase2(vn_type* x, int32_t k)
{

  int32_t *v;        /* vector v of the stack 'x' */
  int32_t *n;        /* vector n of the stack 'x' */
  int32_t size;      /* size of the vector v of stack 'x' */
  int32_t i;
  int32_t j;
  int32_t next;
  int32_t jump;
  int32_t pi;
  int32_t nipi;
  int32_t symbmax;

  v    = GETVECTOR(x,v);
  n    = GETVECTOR(x,n);
  size = x->size - k;
  symbmax = x->symbmax;

  for (j=1; j <= symbmax; j++)
    if (SYMBOLE(v,j) == j)
      {				                 /* start the red chain */
	next = j;
	do
	  {
	    i = next;
	    next = GET(n,i);

	    jump = SYMBOLE(v,i+k);

	    if (jump)
	      {
		pi = GET(v,jump);
                nipi = GET(n,pi);
		if ((! IS_MARKED(n,pi)) ||
		    ( SYMBOLE(v,i) != SYMBOLE(v,nipi)) ||
		    ( jump != SYMBOLE(v,nipi+k)))
		  {
				/* beginning of a green list */
		    SETMARKED(v,jump,i);
		    //MARK(v,jump);

		    SETMARKED(n,i,i);
		    //MARK(n,i);
		  }
		else
		  {
		    SET(n,pi,i);              /* le nouveau dernier  */
		    SET(n,i,GET(n,nipi));  /*  */
		    SET(n,nipi,i);
		    MARK(n,pi);
		  }
	      }

	  } while (( i != next ) && 
		   (next <= size));

      }
}


static void unroll2(vn_type* x)
{
  int32_t *v;
  int32_t *n;
  int32_t size;
  int32_t i;
  int32_t second;
  int32_t last;

  v    = GETVECTOR(x,v);
  n    = GETVECTOR(x,n);
  size = x->size;

  for(i = 1; i <= size; i++)
    if (IS_MARKED(n,i))
      {
	last   = GET(n,i);
	second = GET(n,last);
	SET(n,last,last);
	SETMARKED(n,i,second);
      }
}



static void unstack(vn_type* x, int32_t k)
{
  int32_t *v;
  int32_t *n;
  int32_t size;
  int32_t i;
  int32_t j;
  int32_t classe=0;
  int32_t next;  

  v    = GETVECTOR(x,v);
  n    = GETVECTOR(x,n);
  size = x->size;

  clearVVector(x);

  for(j = 1; j <= size; j++)
    if (IS_MARKED(n,j))
    {

      classe=j;
      next  =j;
    
      do
	{

	  i = next;
	  next = GET(n,i);
	  SET(v,i,classe);
	  
	} while (i != next);
    }

  x->symbmax=classe;
  x->size-=k;
 }


static vn_type *hashSequence(
                 char *sequence,   /* a pointer to the Concateneted sequece */
			     int32_t lmax, 
			     int32_t *k,
                 masked_area_table_t *mask)
{
  vn_type *reponse;                 /* vntype is an object containing two arrays and some values on them */
  int32_t ws;                       /* let see, word size ? */
  int32_t clef;
  int32_t *n;
  int32_t *v;
  int32_t size;
  int32_t symbole;
  int32_t word;
  int32_t i;
  int32_t maskZero;
  int32_t maskClef;
  int32_t zero;
  int32_t count;                    /* Number of Seq\@ within the sequence */
  char *index;
  int32_t symbmax;
  int32_t posinseq;
  
  posinseq=0;

  size=0;
  count=1;
  symbmax = 0;
  for (index=sequence; *index; index++)
    {
      size++;
      if (*index == '@')
          count++;
    }

  reponse = new_vn(size,count);     /* get memory for two vectors and count-1 objects */

  n    = GETVECTOR(reponse,n);      /* init local n and v */
  v    = reponse->v;

                    /* estimate hash word size */
  
  ws    = (int32_t)floor(log((double)size+1)/log(4.0));
  if (ws > lmax)
    ws=lmax;

  if (ws <= 1)
    {
      *k=1;
      return NULL;
    }

  maskClef=0;
  maskZero=0;
  clef    =0;
  zero    =0;
  count   =0;
  
  maskZero= (1 << ws) - 1;
  maskClef= (1 << (ws*2)) - 1;


  for (i=0; i < ws-1; i++)
    {
      if (sequence[i]=='@')
        {
          reponse->limits[count]=i;
          count++;
          posinseq=0;
        }
      symbole = code(sequence[i]);

      if (KMRK_isMasked(mask,count,posinseq))
          symbole=0;

      clef  <<= 2;
      clef   |= symbole-1;

      zero  <<= 1;
      zero   |= (symbole) ? 0:1;
      posinseq++;
    }

  clef   &= maskClef;
  zero   &= maskZero;

  for (i=ws-1; i < size; i++)
    {
      if (sequence[i]=='@')
        {
          reponse->limits[count]=i;
          count++;
          posinseq=0;
        }
      symbole = code(sequence[i]);

      if (KMRK_isMasked(mask,count,posinseq))
          symbole=0;
      
      clef  <<= 2;
      clef   |= symbole-1;
      clef   &= maskClef;


      zero  <<= 1;
      zero   |= (symbole) ? 0:1;
      zero   &= maskZero;

      word = (zero) ? 0:(clef+1);

      if (word)
        {
          symbole    = GET(n,word);
          if (!symbole)
            {
              symbole = i+2-ws;
              SET(n,word,symbole);
            }
        }
      else
        symbole=0;
      
      if (symbole > symbmax)
          symbmax = symbole;
      
      v[i-ws+1]=symbole;
      posinseq++;
    }

  reponse->limits[count]= size;
  reponse->symbmax = symbmax;

  for (i=reponse->size; i < size; i++)
    v[i]=0;

  *k=ws;

  return reponse;
}

static vn_type* encodeSequence( char *sequence,  masked_area_table_t *mask )
{
  int32_t size;
  vn_type *reponse;
  int32_t i;
  int32_t nucleotide;
  int32_t symbole;
  int32_t *n;
  int32_t *v;
  int32_t count;
  char*   index;
  int32_t symbmax;
  int32_t posinseq;
  
  size=0;
  count=1;
  posinseq=0;
  
  for (index=sequence; *index; index++)
    {
      size++;
      if (*index == '@')
          count++;
    }

  reponse = new_vn(size,count);

  n = GETVECTOR(reponse,n);
  v = reponse->v;
  count = 0;
  symbmax = 0;

  for (i=0; i < size; i++)
    {
      if (sequence[i]=='@')
        {
          reponse->limits[count]=i;
          count++;
          posinseq=0;
        }
      nucleotide = code(sequence[i]);
      
      if (KMRK_isMasked(mask,count,posinseq))
          nucleotide=0;
      
      
      if (nucleotide)
        {
          symbole    = GET(n,nucleotide);
          if (!symbole)
            {
              symbole = i + 1;
              SET(n,nucleotide,symbole);
            }
        }
      else
        symbole=0;

      if (symbole > symbmax)
        symbmax = symbole;

      v[i]=symbole;
      posinseq++;
    }

  reponse->limits[count]= size;
  reponse->symbmax = symbmax;
  return reponse;
}


static void applyQuorum(vn_type *x, 
			int32_t count, 
			int32_t countSeq,
			KMRK_QUORUM_FUNC_PTR(quorum))
{
  int32_t *n;
  int32_t *v;
  int32_t size;
  int32_t i;

  n    = GETVECTOR(x,n);
  v    = GETVECTOR(x,v);
  size = x->size;

  for (i=1; i < size; i++)
    if ((IS_MARKED(n,i)) &&
	(! quorum(x,i,count,countSeq)))
      UNMARK(n,i);
}




/********************************************
 ********************************************
 **
 ** Defintion of public functions
 **
 ********************************************
 ********************************************/

/*
	Count the number of direct pairs in the KMR-k result
	input 'x' is the input stack -- see KMRK.h
	vector n is a chained list where each n[i] is the next occurence of i
*/
int32_t KMRK_upperCoupleCount(vn_type *x)
{
  int32_t i;                /* i and j, dummy counters */
  int32_t j;
  int32_t *n;					 /* the n vector of stack 'x' */
  int32_t size;				 /* its size */
  int32_t reponse;			 /* number of pairs  */
  int32_t lllength;			 /* number of occurence */
  int32_t next;				 /* the next occurence in vector 'n' of stack 'x' */
  int32_t xmax;             /* limit of the first sequence: from 1 to xmax */

	xmax = x->limits[0];      

	reponse = 0;

	n    = GETVECTOR(x,n);      /* set n as the vector x->n */
	size = x->size;             /* set its size */

	for (j=1; j < xmax; j++)    /* for all pos of the first sequence */
		if (IS_MARKED(n,j))       /* check if it is marked */
      {
			next = j;
			lllength = 0;
			do 
			{
				i = next;
				next = GET(n,i);   /* get the next occurence position */
				lllength++;
								
			} while ((next <= xmax) && (i != next));   /* while the next occurence are in direct seq *
			                                            * and not themselves */


			reponse += lllength * (lllength - 1) / 2;  /* lllength choose 2 */
			
      }

  return reponse;
}


/*
	Count the number of inverted pairs in the KMR-k result
	input 'x' is the input stack -- see KMRK.h
	vector n is a chained list where each n[i] is the next occurence of i
*/
int32_t KMRK_upperInvertedCount(vn_type* x,int32_t wordsize)
{
  int32_t i;           /* i and j, dummy counters */
  int32_t j;
  int32_t *n;          /* the n vector of stack 'x' */
  int32_t size;        /* its size */
  int32_t reponse;     /* number of pairs with at least one inverted occurence */
  int32_t direct;      /* occurence in the direct sequence */
  int32_t inverted;    /* ccurence in the inverted sequence */
  int32_t next;        /* the next occurence in vector 'n' of stack 'x' */
  int32_t posinv;
  int32_t posdeb;

	/*
		if there was no inverted sequence, there is no inverted seeds
	*/
  if ((x->seqCount == 1) ||   
      (! x->complement))
    return 0;

  reponse = 0;

  n    = GETVECTOR(x,n);                 /* set n as the vector n of stack 'x' */
  size = x->size;                        /* set its size */ 

	for (j=1; j <= x->limits[0]; j++)     /* for all position of the first (direct) sequence */
		if( IS_MARKED(n,j) )               /* if it is tagged (negative)  */
		{
			posdeb   = j;
			next     = j;
			direct   = 0;
			inverted = 0;
			
			do
			{
				i = next;                     /* set i as the position to check */
				next = GET(n,i);              /* next is the symbol n[i] -- pos of the next occurence */

				if ((i > x->limits[0]) &&
				    (i < x->limits[1])    )   /* if it position is in seq inverted */
 					inverted++;             
				else
					direct++;
					
			} while (i != next);             /* when the next occurence is the same pos, stop */

			posinv = 2 * x->limits[0] - i -wordsize + 3;
			if (inverted && (posinv == posdeb))
                           reponse += direct;
			reponse += inverted * direct;    /* the number of occurence for this seed. *
			                                  * NB: if inv is 0, there is no inv occurence */ 
      }

/*  return (reponse+1)/2;   .... weird ?? */

	return reponse / 2;
}

int32_t KMRK_upperInterCount(vn_type* x,int32_t seq1,int32_t seq2,int32_t wordsize)
{
  int32_t i;
  int32_t j;
  int32_t *n;
  int32_t size;
  int32_t reponse;
  int32_t seqcount1;
  int32_t seqcount2;
  int32_t next;

  int min1,min2,max1,max2;

  if (seq1 > seq2)
    {
      i=seq1;
      seq1=seq2;
      seq2=i;
    }
  
  if (seq1==seq2 && seq1 == 0)
    return KMRK_upperCoupleCount(x);

  if (seq1==0 && seq2==1 && x->complement)
    return KMRK_upperInvertedCount(x,wordsize);

  reponse = 0;

  n    = GETVECTOR(x,n);
  size = x->size;

  if (seq1==0)
    min1=0;
  else
    min1 = x->limits[seq1-1];

  max1 = x->limits[seq1];

  min2 =x->limits[seq2-1];
  max2 =x->limits[seq2];

  for (j=1; j <= max1; j++)
	  if (IS_MARKED(n,j))
      {
		  next      = j;
		  seqcount1 = 0;
		  seqcount2 = 0;
		  do 
		  {
			  i = next;
			  next = GET(n,i);
			  if ((i > min1) &&
				  (i <=max1))
				  seqcount1++;
			  else if ((i > min2) &&
					   (i <=max2))
				  seqcount2++;
		  } while ((i != next) && (next <= max2));
		  
		  reponse += seqcount1 * seqcount2;
      }
		  
  /*fprintf(stderr,"\ncouple count : %i\n",reponse);*/
  return reponse;
}


void KMRK_markStart(vn_type* x)
{
	int32_t *v;
	int32_t *n;
	int32_t size;
	int32_t i;

 v    = GETVECTOR(x,v);
 n    = GETVECTOR(x,n);
 size = x->size;

 for(i = 1; i <= size; i++)
    if (SYMBOLE(v,i) == i)
	    MARK(n,i);
    else
	    UNMARK(n,i);

						  
}




/*
	All Quorums functions
	those functions do some check on a given position
*/

int32_t KMRK_CoupleQuorum(              /* return TRUE if array_n[pos] is not linked to itself */
              vn_type* x,
			  int32_t pos,
			  int32_t count, 
			  int32_t countSeq)
{
  return (GET(GETVECTOR(x,n),pos) != pos) ? 1:0;
}

int32_t KMRK_DirInvQuorum(              /* return TRUE if array_n[pos] from DirSeq is linked something else */
              vn_type* x,
			  int32_t pos,
			  int32_t count, 
			  int32_t countSeq)
{
  return  ((pos >= x->limits[0]) ||
           (GET(GETVECTOR(x,n),pos) == pos))  ?  0:1;
}

int32_t KMRK_InvQuorum(              /* return TRUE if array_n[pos] linked to the InvSeq section */
              vn_type* x, 
			  int32_t pos,
			  int32_t count, 
			  int32_t countSeq)
{
  return  ((pos < x->limits[0]) &&
           (GET(GETVECTOR(x,n),pos) >= x->limits[0])) ? 1:0;
}

int32_t KMRK_Int12Quorum(         /* return TRUE if array_n[pos] linked to the Seq2 section */
              vn_type* x, 
			  int32_t pos,
			  int32_t count, 
			  int32_t countSeq)
{
  return  ((pos < x->limits[0]) &&
           (GET(GETVECTOR(x,n),pos) >= x->limits[0])) ? 1:0;
}


int32_t KMRK_IntInv12Quorum(      /* return TRUE if array_n[pos] of DirSeq2 is linked to the InvSeq2 section */
                vn_type* x, 
			    int32_t pos,
			    int32_t count, 
			    int32_t countSeq)
{
  return  ((pos <  x->limits[1]) &&
           (pos >= x->limits[0]) &&
           (GET(GETVECTOR(x,n),pos) >= x->limits[1])) ? 1:0;
}

int32_t KMRK_IntDirInv12Quorum(      /* return TRUE if array_n[pos] is linked to the InvSeq2 section */
                   vn_type* x, 
			       int32_t pos,
			       int32_t count, 
			       int32_t countSeq)
{
  return  ((pos < x->limits[1]) &&
           (GET(GETVECTOR(x,n),pos) >= x->limits[1])) ? 1:0;
}








void KMRK_RunKMRK(vn_type *x,                     /* an object with all arrays and associated values */  
                  int32_t k,                      /* the current word size value */
                  int32_t count,                  /* number of occurence */
                  int32_t countSeq,               /* the number of sequences */
                  KMRK_QUORUM_FUNC_PTR(quorum))   /* a pointer to the quorum */
{
  clearNMark(x);
  phase2(x,k);
  applyQuorum(x,count,countSeq,quorum);
  unroll2(x);
  unstack(x,k);
}



/** 
 * Initialise a vn_type record from one sequence to run KMRK algorithm.
 * In more than KMRK_InitOneSequence, KMRK_HashOneSequence construct 
 * word of len k with an hash algorithm. k used is a function of
 * sequence size and alphabet size. If calculed k is superior to lmax
 * then k = lmax. 
 * 
 * @param sequence   pointer to a C string containing the sequence
 * @param complement != 0 means that seq one and two are the two strands of 
 *                   the same sequence.
 * @param count      parameter count passed to the quorun function
 * @param countSeq   parameter countSeq passed to the quorun function
 * @param k          maximum length of the created word (input)
 *                   length of the word represented by each symbole 
 *                   of v (output). 
 * @param quorum     pointer to a quorum function
 * 
 * @return a pointer to vn_type structure correctly initialized 
 *          to be used by KMRK_RunKMRK
 *
 * @see   KMRK_InitOneSequence
 */

vn_type *KMRK_HashOneSequence(char *sequence,                 /* the sequence either "DirSeq\0" or "DirSeq\@InvSeq\0" */          
                              int32_t complement,             /* set to 1 if the latter is true */
                              int32_t count,                  /* number of occurence, always 2 - unused */
                              int32_t countSeq,               /* Number of Sequence - should be 1 here !! */
                              int32_t *k,                     /* currenty word size - starts with Lmin value */
                              KMRK_QUORUM_FUNC_PTR(quorum),   /* a pointer to the quorum function */
                              masked_area_table_t *mask)      /* if mask option is set on */
{

  int32_t i;
  int32_t *n;
  int32_t *v;
  int32_t size;
  vn_type *reponse;                                       /* a structure with two arrays and associated values */
  int32_t lmax;

  lmax = *k;                                              /* do not hash larger than lmin */

  reponse = hashSequence(sequence,lmax,k,mask);           /* hash the Sequence ans fill up the vn_type structure */

  if (!reponse)                                           /* if HASH failed then encode the sequence */
    {
      reponse = encodeSequence(sequence,mask);
      *k=1;
    }

  reponse->complement = complement;                        /* does the sequence contains also InvSeq */

  phase1(reponse);

  n = GETVECTOR(reponse,n);
  v = GETVECTOR(reponse,v);
  size = reponse->size;

  for (i=1; i <= size; i++)
    if (SYMBOLE(v,i) == i)
      MARK(n,i);
  
  applyQuorum(reponse,count,countSeq,quorum);

  unroll2(reponse);
  unstack(reponse,(*k)-1);
  
  /*  fprintf(stderr, "Hash up to size %d\n",*k); */

  return reponse;
}


/** 
 * realises serveral run of KMR cycle to make from a sequence
 * a vn_type structure describing sequences of factors of a precise size.
 * 
 * @param sequence pointer to a C string containing the sequence
 * @param size word size to construct
 * @param count parameter count passed to the quorun function
 * @param countSeq parameter countSeq passed to the quorun function
 * @param quorum pointer to a quorum function
 * @param init pointer to a initialisation  function
 * 
 * @return a vn_type pointer to a structure containing sequences of factors
 */

vn_type *KMRK_RunTo(char *sequence,                      /* the concantenated seq ("DirSeq\0" OR "DirSeq\@InvSeq\0" */
                    int32_t size,                        /* lmin size */
                    int32_t complement,                  /* if 1, sequence is "DirSeq\@InvSeq\0" */ 
                    int32_t count,                       /* the number of occurence - always set to 2 - unused */
                    int32_t countSeq,                    /* number of sequence */
                    KMRK_QUORUM_FUNC_PTR(quorum),        /* which quorum */
                    KMRK_INIT_FUNC_PTR(init,quorum),     /* a pointer to the Hash function */
                    masked_area_table_t *mask)           /* if option mask is set on */
{
  
  
  int32_t k;               /* the current word size of merge */
  vn_type *reponse;        /* a pointer to the object with 2 arrays and associated values */

  k=size;                  /* set to final size - use for the first Hash */


  reponse = init(sequence, complement, count, countSeq, &k, quorum, mask);

  while(k <= size/2)
    {
      /*fprintf(stderr,"KMRK step %d\n",k);*/
      KMRK_RunKMRK(reponse,k,count,countSeq,quorum);
      k*=2;
    }

  if (k < size)
    KMRK_RunKMRK(reponse,size - k,count,countSeq,quorum);

  return reponse;
}


void KMRK_FreeVN(vn_type *x)
{

  if (x)
    {
      if (x->v) MyFree(x->v, (x->vectorsize)*sizeof(int32_t) );
      if (x->n) MyFree(x->n, (x->vectorsize)*sizeof(int32_t) );
      MyFree(x, 1* sizeof(vn_type)+sizeof(int32_t) * (x->seqCount-1)  );
    }
}



/*
	Old Unused Functions
*/


/*
static int32_t *copyNVector(vn_type *x, int32_t *dest)
{
  int32_t* reponse;

  if (dest)
    reponse=dest;
  else
    KMRK_MALLOC(reponse,int32_t,x->size * sizeof(int32_t),
		"copyNVector: I cannot not allocate a copy "
		"of n vector");

  memcpy(reponse,x->n,x->size*sizeof(int32_t));
  
  return reponse;
}

	
vn_type* KMRK_InitOneSequence(char *sequence,              
                              int32_t complement, 
                              int32_t count,
                              int32_t countSeq,
                              int32_t *k,
                              KMRK_QUORUM_FUNC_PTR(quorum),
                              masked_area_table_t *mask)
{
  
  vn_type *reponse;
  int32_t i;
  int32_t *n;
  int32_t *v;
  int32_t size;

  reponse = encodeSequence(sequence,mask);
  reponse->complement = complement;

  phase1(reponse);

  n = GETVECTOR(reponse,n);
  v = GETVECTOR(reponse,v);
  size = reponse->size;

  for (i=1; i <= size; i++)
    if (SYMBOLE(v,i) == i)
      MARK(n,i);

  applyQuorum(reponse,count,countSeq,quorum);

  unroll2(reponse);
  unstack(reponse,0);
  *k=1;
  return reponse;

}

static void displayVN(vn_type * x)
{
  int32_t i;
  
  for (i=0; i < x->size; i++)
    {
      printf("%-4i => %-4i ; %-4i\n",i+1,x->v[i],x->n[i]);
    }

  printf("\n");
}*/



