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
 * @file   smin.h
 * @author Guillaume Achaz <gachaz@oeb.harvard.edu>
 * @date   April 2004
 * 
 * @brief  header for all shuffling/smin process
 * 
 * 
 */


#ifndef _smin_h_
#define _smin_h_


/*
	From a sequence, get length and dGC
	then compute the Smin using Waterman & Vingron
*/
double compute_smin(char *sequence, double pvalue, int32_t size, int32_t nX);


/*
	From two sequences, get average length and average dGC
	then compute the Smin using Waterman & Vingron
*/

double compute_smin_2seq(char *sequence1, char *sequence2, int32_t size1,  int32_t size2,   double pvalue, int32_t nX,  int32_t nX2);


/*
	tag all repeats with a score < smin
	set their positions and length at -1
*/
void UpdateRepeats(Repeats *AllRepeats, double Score_min);


#endif /* _smin_h_*/
