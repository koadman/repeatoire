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
/**
 * @file   repseek_types.h
 * @author Guillaume Achaz <gachaz@oeb.harvard.edu>
 * @date   April 2004
 * @modif  July 2004 turn scores into doubles
 * @brief  definition of general types and macros for repseek
 * 
 * 
 */



#ifndef _REPSEEK_TYPES_
#define  _REPSEEK_TYPES_

/*
	Version of the program
*/
#define REPSEEK_VERSION "6.5"
#define REPSEEK_DATE "Dec, 09, 2009"



/********** **********

        General Macros

 ********** **********/

/*
	Macros to compare 2 or three values
*/	
#define MAX2( A, B )     ( ((A)>(B))?(A):(B) )
#define MAX3( A, B, C )  ( ((MAX2(A,B))>(C))?(MAX2(A,B)):(C) )
#define MIN2( A, B )     ( ((A)<(B))?(A):(B) )

#define MAX( A, B )     MAX2( A, B )
#define MIN( A, B )     MIN2( A, B )

/*
	Absolute values
*/
#define ABS(x) (( (x)>=0 )? (x) : -(x))



/********** **********

   All types used in repseek

 ********** **********/

#include <stdio.h>     /* The type FILE * is defined here */

#ifdef OSF
typedef signed char          int8_t;
typedef short                int16_t; 
typedef long                 int32_t;
#else
#include <stdint.h>           /* all, the int??_t are defined in there for typical unix */
#endif

/**
 * Store informations about one STRICT repeat (seeds)
 *
 */
typedef struct {      /* the complete seed structure */

  int32_t pos1;       /**< position of the first copy */
  int32_t pos2;       /**< position of the second copy */

  int32_t length;     /**< length of the strict repeats */
  float   rmean;      /**< mean repeat leavel */
  
} Seed_type;

typedef struct {      /* Just after a KMRK length X, only the 2 pos matter */

  int32_t pos1;       /**< postion of the first copy */
  int32_t pos2;       /**< postion of the second copy */

} SmallSeed_type;   




/**
 * Store informations about all strict repeat (seeds)
 *
 */
typedef struct {

  int32_t      cDirSeeds;          /**< currently allocated space in dirSeeds array */
  int32_t      nDirSeeds;          /**< count of direct strict repeats */
  int32_t      nFilteredDirSeeds;  /**< ??? */

  Seed_type*   dirSeeds;           /**< array of direct repeats */

  int32_t      cInvSeeds;          /**< currently allocated space in invSeeds array */
  int32_t      nInvSeeds;          /**< count of inverted strict repeats */
  int32_t      nFilteredInvSeeds;  /**< ??? */

  Seed_type*   invSeeds;           /**< array of inverted repeats */

} AllSeeds_type;





/**
 * Store informations about one GENERIC repeat
 *
 */
typedef struct{              

  char type[20];               /* its name; i.e. Tandem, Palindrome, etc... */
  
  int32_t pos1, pos2,          /* both copies postions */
        len1, len2,            /* both copies length */
        seed_pos1,seed_pos2,   /* pos1 and pos2 of the originate seed */
        seed_len,              /* len  of  the seed	*/
		  
  	match, align_len;           /* number of match and length of alignment */
	
	double score;		          /* the alignment score */

	float seed_meanR;           /* the seed meanR */

	float meanR;                /* The mean R-level of the repeat */
	int32_t mainR;              /* its Mode R */
	float fraction_mainR;       /* the fraction of length containing the Mode R */


} Rep;



/**
 * Store informations about All GENERIC repeats
 *
 */
typedef struct {  
			  
	int32_t nDirRep;	        /* Total Number of Direct Repats in Mem */
	int32_t nDirBadRep;	        /* Direct repeats set to -1 -- filtered out and co. */
	Rep *DirRep;                /* The array of Dir Rep */

	int32_t nInvRep;	        /* Total Number of Inverted Repats in Mem */
	int32_t nInvBadRep;		    /* Inverted Repeats set to -1 -- filtered out and co. */
	Rep *InvRep;			    /* The array of Inverted Rep */

} Repeats;


#define MATRIX_SIZE 26

typedef struct {            /******* The scoring Matrix  ************/

	double matrix[MATRIX_SIZE][MATRIX_SIZE];  /* the matrix of match/missmatch */
	double gap_open;                          /* value of gap-open */
	double gap_ext;                           /* value of gap_ext */
	double expect;

} SCORING;



typedef struct {            /******* The Results of Alignement by Dynamik programming  ************/


	double *scores;                 /* the score strings (+/- 'matrix') */
	double *pscore;                 /* pointer to the current score */
	
	double *F;                      /* *F is a standard to describe the best score from top */

	int32_t *traces;                /* the path matrix - could take values 0 for diagonal, >0 for top and <0 for left */
	int32_t *Top;                   /* *Top is the number of steps from TOP -- needed for memorizing deletion in seq2 */
	int32_t Left;
	
	
	double *pBestScore;             /* pointer to it */
	int32_t BestScore_row;          /* its row and col */
	int32_t BestScore_col;   

	char *traceback;                /* all you need for bactracking */
	char *traceback_f;              /* memory for forward traceback and then check other seeds */
	char *traceback_b;              /* memory needed for backward traceback - to avoid erasing the forward one */

	int32_t alignment_len;          /* guess ?? */
	int32_t matches;                /* number of match (score>0 in scoring matrix) */

	int32_t nSegment;               /* number of segment */
	int32_t *Segment_begin;         /* begin and end of each segment */
	int32_t *Segment_end;
   
	int32_t max_scores;             /* size of the matrices ; only for memory purposes */
	int32_t max_col;
	int32_t max_row;
	int32_t max_alignlen;

} RESULTS;



#endif  /* _REPSEEK_TYPES_ */

