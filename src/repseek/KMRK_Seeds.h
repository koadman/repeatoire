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
#ifndef KMRK_Seeds_h
#define KMRK_Seeds_h

/********************************************
 ********************************************
 **
 ** Declaration of struct
 **
 ********************************************
 ********************************************/

#include <stdio.h>
#include "KMRK.h"
#include "repseek_types.h"




#define KMRK_SORT_SEEDS_FUNC_PTR(name) int32_t (*name)(SmallSeed_type*, \
						       SmallSeed_type*)

#define KMRK_DELTA_SEEDS_FUNC_PTR(name) int32_t (*name)(SmallSeed_type*)

/********************************************
 ********************************************
 **
 ** Declaration of public functions
 **
 ********************************************
 ********************************************/


AllSeeds_type *KMRK_allocSeeds(AllSeeds_type *AllSeeds, 
			       int32_t size, 
			       int8_t opt_dir, 
			       int8_t opt_inv);

void KMRK_SetMultipleLenDirSeeds(SmallSeed_type* seeds,
                                 int32_t         nseeds,
                                 int32_t         wordSize,
                                 AllSeeds_type   *PtrAllSeeds);


void KMRK_freeSeeds(AllSeeds_type *AllSeeds);

void KMRK_compactSeeds(AllSeeds_type *AllSeeds);

void KMRK_pushSeed(AllSeeds_type *AllSeeds,
		   int32_t       pos1,
		   int32_t       pos2,
		   int32_t       length,
		   int8_t        dir);

AllSeeds_type* KMRK_enumerateDirectCouple(AllSeeds_type* Seeds,
					  int32_t expected,
					  int32_t wordSize,
					  vn_type* stack,
					  int32_t seq);

AllSeeds_type* KMRK_enumerateInvertedCouple(AllSeeds_type* Seeds,
					    int32_t expected,
					    int32_t wordSize,
					    vn_type* stack);

AllSeeds_type* KMRK_enumerateInterCouple(AllSeeds_type* Seeds,
					 int32_t seq1,
					 int32_t seq2,
					 int32_t expected,
					 int32_t wordSize,
					 vn_type* stack);

AllSeeds_type* KMRK_enumerateInterInvertedCouple(AllSeeds_type* Seeds,
						 int32_t seq2,
						 int32_t expected,
						 int32_t wordSize,
						 vn_type* stack);


/** 
 * Compare two seeds and return an integer less than, equal to or greater 
 * than zero considering the relative order of the two seeds. This
 * version take into account only pos1 and pos2 of seeds without taking
 * account of the sequences or of the relative direction
 * 
 * @param s1 pointer to seed one
 * @param s2 pointer to seed two
 * 
 * @return a integer less than, equal to or greater than zero
 */

int32_t KMRK_cmpSeedsPos(SmallSeed_type *s1, SmallSeed_type *s2);
int32_t KMRK_cmpDeltaSeedsPos(SmallSeed_type *s1, SmallSeed_type *s2);
int32_t KMRK_cmpDeltaInvSeedsPos(SmallSeed_type *s1, SmallSeed_type *s2);

void KMRK_sortSeeds(SmallSeed_type* seeds,
		    int32_t         nseeds,
		    KMRK_SORT_SEEDS_FUNC_PTR(compare));


AllSeeds_type* KMRK_get_seeds(char **seq,
                              int32_t SimpleSeqLen,
                              int16_t Lmin,
                              int8_t opt_dir,
                              int8_t opt_inv, 
                              int8_t opt_verbose,
                              masked_area_table_t *mask);

AllSeeds_type* KMRK_get_seeds_2seqs(char **seq1,
                                    char **seq2,
                                    int32_t size1,
                                    int32_t size2,
                                    int16_t Lmin,
                                    int8_t opt_dir,
                                    int8_t opt_inv, 
                                    int8_t opt_verbose,
                                    masked_area_table_t *mask);

/** 
 * Order an array of seeds by pos1,pos2
 * 
 * @param seeds pointer to an array of Seed_type object to sort
 * @param count count of element in the array
 */
void KMRK_sortSeedsByPos(Seed_type* seeds, int32_t count);

#endif   /* KMRK_Seeds_h */
