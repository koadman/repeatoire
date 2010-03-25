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
#ifndef _ALIGN_H_
#define _ALIGN_H_

/* 
	Macros to turn any char from A-Z to a symb 0-26 (for S&W scoring matrix)
*/
#define CHAR2SYMB( A ) ( (A) - 65 )
#define SYMB2CHAR( i ) ( (i) + 65 )

#define Match(A,B) pScoring->matrix[ CHAR2SYMB(seq1[ (A)-1 ]) ][ CHAR2SYMB(seq2[ (B)-1 ]) ]


/*
	The Seed2Repeat part
*/
Repeats align_seeds(AllSeeds_type *AllSeeds, char *sequence,  float Xg, float gap_open, float gap_ext,\
                    char matrix_type, int16_t Lmin, int8_t opt_dir, int8_t opt_inv, int8_t opt_overlap, float merge_repeats);

Repeats align_seeds_2seq(AllSeeds_type *AllSeeds, char *seq1, char *seq2,  float Xg, float gap_open, float gap_ext, char matrix_type,\
                    int16_t Lmin, int8_t opt_dir, int8_t opt_inv, float merge_repeats);


Rep alignd(int32_t pos1, int32_t pos2, int32_t len, char *sequence, int32_t sizeseq,
           float Xg, SCORING *pScoring, RESULTS *pResults , int8_t opt_overlap);

Rep aligni(int32_t pos1, int32_t pos2, int32_t len, char *sequence, char *invsequence,  int32_t sizeseq,\
           float Xg, SCORING *pScoring, RESULTS *pResults );

Rep alignd_2seq(int32_t pos1, int32_t pos2, int32_t len, char *seq1, char *seq2, \
              int32_t size1, int32_t size2, float Xg, SCORING *pScoring, RESULTS *pResults );

Rep aligni_2seq(int32_t pos1, int32_t pos2, int32_t len, char *seq1, char *invseq2, \
              int32_t size1, int32_t size2, float Xg, SCORING *pScoring, RESULTS *pResults );




/*
	The alignemnt itself
*/

void align_blast2like(char *seq1, char *seq2, float Xg, SCORING *pScoring,\
                      RESULTS *pResults, char direction, int32_t max_row, int32_t max_col, int32_t diff_max);
void log_matrix(char *sequence, float gap_open, float gap_ext, SCORING *pScoring );
void identity_matrix(char *sequence, float gap_open, float gap_ext, SCORING *pScoring);


#endif /*_ALIGN_H_*/
