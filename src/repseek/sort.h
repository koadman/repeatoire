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
 * @file   sort.h
 * @author amikezor <gachaz@oeb.harvard.edu>
 * @date   april 2004
 * 
 * @brief   header fr sort.c
 * 
 * 
 */
 
 
#ifndef _SORT_H_
#define _SORT_H_

#include <stdio.h>
#include <stdlib.h>
#include "repseek_types.h"
#include "memory.h"


/*
	A simple qsort as described in numrec
*/
void double_Qsort(double *array1, int32_t n);


/***
	This function sort the first array, and for equal value, use the second one
	calling the qsort2() function
***/
void qsort2_1then2(int32_t *array1, int32_t *array2, int32_t n);


/***
	This function sort by pos1 and for equal value, use pos2
	calling the qsort2_repeats() function
***/

void qsort_1then2_repeats(Rep *AllReps, int32_t n);

#endif /* _SORT_H_*/
