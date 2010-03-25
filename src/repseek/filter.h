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
#ifndef KMRK_filter_h
#define KMRK_filter_h


#include "KMRK.h"
#include "KMRK_Seeds.h"
#include <stdio.h>

void KMRK_FamilySeeds(AllSeeds_type *AllSeeds, 
		      double min, 
		      double Max, 
		      int8_t opt_dir, 
		      int8_t opt_inv, 
		      int32_t size_chr);

int32_t *KMRK_SeedTableR(AllSeeds_type *AllSeeds, 
				int8_t opt_dir, 
				int8_t opt_inv, 
				int32_t size_chr);

void KMRK_SeedTableR2seq(AllSeeds_type *AllSeeds, 
			 int8_t opt_dir, 
			 int8_t opt_inv, 
			 int32_t size1,
			 int32_t size2,
			 int32_t **r1,
			 int32_t **r2);

void KMRK_FamilySeeds2seq(AllSeeds_type *AllSeeds, 
		      double min, 
		      double Max, 
		      int8_t opt_dir, 
		      int8_t opt_inv, 
		      int32_t size1,
			  int32_t size2);


void BuiltFilterFamily_seeds(AllSeeds_type *AllSeeds, 
			 int32_t *table_r, 
			 double min, 
			 double Max, 
			 int8_t opt_dir, 
			 int8_t opt_inv);
			 
void BuiltFilterFamily_seeds2seq(AllSeeds_type *AllSeeds, 
			     int32_t *table_r1, 
			     int32_t *table_r2,
			     double min, 
			     double Max, 
			     int8_t opt_dir, 
			     int8_t opt_inv);



#endif /* KMRK_filter_h */
