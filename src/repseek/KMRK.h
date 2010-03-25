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
#ifndef KMRK_h
#define KMRK_h


/********************************************
 ********************************************
 **
 ** Declaration of struct
 **
 ********************************************
 ********************************************/

#include "repseek_types.h"
#include "KMRK_mask.h"

/** 
 * Structure used to manipulate simultanously 
 * the v and n vector 
 * 
 */ 

typedef struct {

	int32_t  size;           /**< size of the current sequence to search for */
	int32_t  vectorsize;     /**< size of vectors */
	
	int32_t  seqCount;       /**< count of concatenated sequences */
	int32_t  complement;     /**< if seqCount > 1 and complement !=0
                              *   then second sequence is the inverted complement strand of first one */
	int32_t  symbmax;
	int32_t* v;               /**< sequence vector */
	int32_t* n;	              /**< linked list vector the special vector of KRM-K */
	int32_t  limits[1];       /**< array of end limits of concatenated sequences in v (array size is seqCount) */
	
} vn_type;



/********************************************
 ********************************************
 **
 ** Declaration of public macro
 **
 ********************************************
 ********************************************/

#define GETVECTOR(x,v)  (((x)->v) - 1)

#define IS_MARKED(x,i)  ((x)[i] < 0)
#define MARK(x,i)       ((x)[i]) = -ABS((x)[i])
#define UNMARK(x,i)     ((x)[i]) = ABS((x)[i])

#define SET(x,i,v)      ((x)[i]) = (v)
#define SETMARKED(x,i,v) ((x)[i]) = -(v)
#define GET(x,i)        ABS((x)[i])
#define SYMBOLE(x,i)    ((IS_MARKED((x),(i))) ? (i): (GET(x,i)))



/** 
 * Macro used to declare a pointer to a quorum function.
 * 
 * @param name name of the pointer
 * 
 */

#define KMRK_QUORUM_FUNC_PTR(name) int32_t (*name)( vn_type* x, int32_t pos, int32_t count, int32_t countSeq )




/** 
 * Macro used to declare a pointer to an initialisation function.
 * 
 * @param name name of the pointer
 * @param quorum name used for the quorum assiciated function
 *
 * @see KMRK_QUORUM_FUNC_PTR
 * 
 */

#define KMRK_INIT_FUNC_PTR(name,quorum)  vn_type* (*name)(char *sequence,              \
                                                          int32_t complement,          \
                                                          int32_t count,               \
                                                          int32_t countSeq,            \
                                                          int32_t *k,                  \
                                                          KMRK_QUORUM_FUNC_PTR(quorum),\
                                                          masked_area_table_t *mask)


/********************************************
 ********************************************
 **
 ** Declaration of public functions
 **
 ********************************************
 ********************************************/





/** 
 * Initialise a vn_type record from one sequence to run KMRK algorithm
 * 
 * @param sequence pointer to a C string containing the sequence
 * @param complement != 0 means that seq one and two are the two strands of 
 *                   the same sequence.
 * @param count parameter count passed to the quorun function
 * @param countSeq parametter countSeq passed to the quorun function
 * @param k length of the word represented by each symbole of v. 
 *        k is an output parametter
 * @param quorum pointer to a quorum function
 * 
 * @return a pointer to vn_type structure correctly initialized 
 *          to be used by KMRK_RunKMRK
 *
 * @see KMRK_HashOneSequence
 */

vn_type* KMRK_InitOneSequence(char *sequence,              
                              int32_t complement, 
                              int32_t count,
                              int32_t countSeq,
                              int32_t *k,
                              KMRK_QUORUM_FUNC_PTR(quorum),
                              masked_area_table_t *mask);





/** 
 * Initialise a vn_type record from one sequence to run KMRK algorithme.
 * In more than KMRK_InitOneSequence, KMRK_HashOneSequence construct 
 * word of len k with an hash algorithm. k used is a function of
 * sequence size and alphabet size. If calculed k is superior to lmax
 * then k = lmax. 
 * 
 * @param sequence   pointer to a C string containing the sequence
 * @param complement != 0 means that seq one and two are the two strands of 
 *                   the same sequence.
 * @param count      parametter count passed to the quorun function
 * @param countSeq   parametter countSeq passed to the quorun function
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

vn_type* KMRK_HashOneSequence(char *sequence,              
                              int32_t complement, 
                              int32_t count,
                              int32_t countSeq,
                              int32_t *k,
                              KMRK_QUORUM_FUNC_PTR(quorum),
                              masked_area_table_t *mask);


/** 
 * An example of quorum function testing than a factor is
 * present a least two times. Because of definition of this
 * quorum function count and countSeq parametter have no meanning
 * in this instance of quorum function
 * 
 * @param x   a pointer to vn_type structure to check
 * @param pos position in n of the beginning of the linked list to test
 * @param count minimal number of occurence of factor
 * @param countSeq minimal number of sequences concerned
 * 
 * @return 1 if quorum is ok 0 otherwise.
 */

int32_t KMRK_CoupleQuorum(vn_type* x, 
			  int32_t pos,
			  int32_t count, 
			  int32_t countSeq);


/** 
 * An example of quorum function testing than a factor is
 * present a least two times in the direct strand of a sequence or
 * at least one time in the direct strand and one time in the reverse
 * strand. Because of definition of this
 * quorum function count and countSeq parametter have no meanning
 * in this instance of quorum function
 * 
 * @param x   a pointer to vn_type structure to check
 * @param pos position in n of the beginning of the linked list to test
 * @param count minimal number of occurence of factor
 * @param countSeq minimal number of sequences concerned
 * 
 * @return 1 if quorum is ok 0 otherwise.
 */

int32_t KMRK_DirInvQuorum(vn_type* x, 
			  int32_t pos,
			  int32_t count, 
			  int32_t countSeq);

/** 
 * An example of quorum function testing than a factor is
 * present a least one time in the direct strand and one time in the reverse
 * strand. Because of definition of this
 * quorum function count and countSeq parametter have no meanning
 * in this instance of quorum function
 * 
 * @param x   a pointer to vn_type structure to check
 * @param pos position in n of the beginning of the linked list to test
 * @param count minimal number of occurence of factor
 * @param countSeq minimal number of sequences concerned
 * 
 * @return 1 if quorum is ok 0 otherwise.
 */

int32_t KMRK_InvQuorum(vn_type* x, 
		       int32_t pos,
		       int32_t count, 
		       int32_t countSeq);

int32_t KMRK_Int12Quorum(vn_type* x, 
			 int32_t pos,
			 int32_t count, 
			 int32_t countSeq);

int32_t KMRK_IntInv12Quorum(vn_type* x, 
			    int32_t pos,
			    int32_t count, 
			    int32_t countSeq);

int32_t KMRK_IntDirInv12Quorum(vn_type* x, 
			       int32_t pos,
			       int32_t count, 
			       int32_t countSeq);
/** 
 * realize one cycle of KMR.
 * 
 * @param x a pointer to vn_type created by an intialisation
 *           function or returned by this function.
 * @param k step used to join two words
 * @param count parametter count passed to the quorun function
 * @param countSeq parametter countSeq passed to the quorun function
 * @param KMRK_QUORUM_FUNC_PTR quorum pointer to a quorum function
 */

void KMRK_RunKMRK(vn_type *x, 
		  int32_t k, 
		  int32_t count,
		  int32_t countSeq,
		  KMRK_QUORUM_FUNC_PTR(quorum));




/** 
 * realises serveral run of KMR cycle to make from a sequence
 * a vn_type structure describing sequences of factors of a precise size.
 * 
 * @param sequence pointer to a C string containing the sequence
 * @param size word size to construct
 * @param count parametter count passed to the quorun function
 * @param countSeq parametter countSeq passed to the quorun function
 * @param quorum pointer to a quorum function
 * @param init pointer to a initialisation  function
 * 
 * @return a vn_type pointer to a structure containing sequences of factors
 */

vn_type *KMRK_RunTo(char *sequence,
                    int32_t size,
                    int32_t complement,
                    int32_t count,
                    int32_t countSeq,
                    KMRK_QUORUM_FUNC_PTR(quorum),
                    KMRK_INIT_FUNC_PTR(init,quorum),
                    masked_area_table_t *mask);



/** 
 * free memory associated to a vn_type pointer
 * 
 * @param x a pointer to vn_type structure
 */

void KMRK_FreeVN(vn_type *x);


int32_t KMRK_upperCoupleCount(vn_type *x);
int32_t KMRK_upperInvertedCount(vn_type* x,int32_t wordsize);
int32_t KMRK_upperInterCount(vn_type* x,int32_t seq1,int32_t seq2,int32_t wordsize);

void KMRK_markStart(vn_type* x);

#endif  /* KMRK_h */
