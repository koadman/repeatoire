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
 * @file   align_memory.h
 * @author amikezor <gachaz@oeb.harvard.edu>
 * @date   April 2004
 * 
 * @brief  header for all memory alloc/dealloc for the alignemnt part.
 * 
 * 
 */
 
#ifndef _align_memory_h_
#define _align_memory_h_


void mem_Blast2Align(RESULTS *pResults);

void free_Blast2Align(RESULTS *pResults);

void ReallocScores(RESULTS *pResults, int32_t sig);

void ReallocTraceback(RESULTS *pResults, int32_t sig, char direction);

void ReallocAbove(RESULTS *pResults, int32_t sig);

void ReallocBegEnd(RESULTS *pResults, int32_t sig);


#endif /* _align_memory_h_ */
