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
/******
	file :     align_di.c
	function : get a seeds, align it forward and backward and return the aligned repeat
	
	created  :  Apr 11 2001
	modif    :  Jul 2003,Feb 2004, April 2004, July 2004: turn scores into doubles
	 
	author :   amikezor
*****/


#include "repseek_types.h"
#include "align.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


/*
	Get the score of a seed
	first and ast position are not taken into account since they are used in the alignments
*/
static double score_seed(int32_t pos, int32_t len, char *sequence, SCORING *pScoring){

	char *pseq;              /* pointer to the seqeunce */
	double score=0;          /* the score */
	int32_t index=0;         /* the symbol of the character in the reward matrix */
	
	pseq = sequence+pos;
	
	while(pseq < sequence+pos+len-1){
		
		index =  CHAR2SYMB(*pseq);
		score += pScoring->matrix[index][index];
		pseq++;
	}
	
	return score;
}


/*
	The repeat is a palindrome
*/
static Rep Palindome(int32_t pos1, int32_t pos2, int32_t len,  char *sequence, SCORING *pScoring){

	Rep aligned_seed;

	aligned_seed.seed_pos1 = pos1;
	aligned_seed.seed_pos2 = pos2;
	aligned_seed.seed_len = len;
	aligned_seed.pos1 = pos1;
	aligned_seed.len1 = len/2;
	aligned_seed.pos2 = pos1+aligned_seed.len1;
	aligned_seed.len2 = len/2;
  	aligned_seed.match = len/2;
	aligned_seed.align_len = len/2;            
	strcpy(aligned_seed.type, "Palindrome");
	aligned_seed.score = score_seed(pos1, len, sequence, pScoring);
	
	return aligned_seed;
	
}

/*
	Init the repeats value with a seed
	first and last position are removed since they are used in the alignements
*/
static void InitRep(int32_t pos1, int32_t pos2, int32_t len, double ScoreSeed, Rep *aligned_seed ){

	aligned_seed->pos1=pos1+1;
	aligned_seed->pos2=pos2+1;

	aligned_seed->len1= len-2;   /* len of seed + aligned_genomik */
	aligned_seed->len2= len-2;
	
	aligned_seed->seed_pos1 = pos1;
	aligned_seed->seed_pos2 = pos2;
	aligned_seed->seed_len = len;

	aligned_seed->match = len-2;
	aligned_seed->align_len = len-2;
	aligned_seed->score = ScoreSeed;
}

/*
	Update the repeats from we have on the right - forward-
	Only for direct repeats
*/
static void UpdateRepRight(Rep *aligned_seed , RESULTS *pResults){
	
	aligned_seed->len1 += pResults->BestScore_row+1;   /* len of seed + aligned_genomik */
	aligned_seed->len2 += pResults->BestScore_col+1;
	
	aligned_seed->match += pResults->matches;
	aligned_seed->align_len += pResults->alignment_len;
	aligned_seed->score += *pResults->pBestScore;
}

/*
	Update the repeats from we have on the left - backward -
	Only for direct repeats
*/
static void UpdateRepLeft(Rep *aligned_seed , RESULTS *pResults){

	aligned_seed->pos1 -= pResults->BestScore_row+1;   /* len of seed + aligned_genomik */
	aligned_seed->pos2 -= pResults->BestScore_col+1;

	aligned_seed->len1 += pResults->BestScore_row+1;   /* len of seed + aligned_genomik */
	aligned_seed->len2 += pResults->BestScore_col+1;
	
	aligned_seed->match += pResults->matches;
	aligned_seed->align_len += pResults->alignment_len;
	aligned_seed->score += *pResults->pBestScore;
}

/*
	Same functions but for inverted repeats
*/
static void UpdateRepRightInv(Rep *aligned_seed , RESULTS *pResults){

	aligned_seed->len1 += pResults->BestScore_row+1;   /* len of seed + aligned_genomik */

	aligned_seed->pos2 -= pResults->BestScore_col+1;
	aligned_seed->len2 += pResults->BestScore_col+1;
	
	aligned_seed->match += pResults->matches;
	aligned_seed->align_len += pResults->alignment_len;
	aligned_seed->score += *pResults->pBestScore;
}
static void UpdateRepLeftInv(Rep *aligned_seed , RESULTS *pResults){

	aligned_seed->pos1 -= pResults->BestScore_row+1;   /* len of seed + aligned_genomik */

	aligned_seed->len1 += pResults->BestScore_row+1;   /* len of seed + aligned_genomik */
	aligned_seed->len2 += pResults->BestScore_col+1;
	
	aligned_seed->match += pResults->matches;
	aligned_seed->align_len += pResults->alignment_len;
	aligned_seed->score += *pResults->pBestScore;
}


/***
	The function that return an aligned Rep from a cupple of seed information
	Moreover all Struct have to be init before
***/
Rep alignd(int32_t pos1, int32_t pos2, int32_t len, char *sequence,  int32_t sizeseq , 
           float Xg, SCORING *pScoring, RESULTS *pResults , int8_t opt_overlap){

		      
	Rep aligned_seed;               /* the repeat */
	char *pseq1, *pseq2;            /* pointer for the first and the second copy */
	double ScoreSeed=0;             /* score of the seed */
	int32_t max_col=0, max_row=0;   /* the dist max from the edge of the sequence max for rwo and col */
	int32_t diff_max =0;            /* the maximum difference between row and col to avoid aligning the chr with itself */



	/*
		Init by the seed
	*/			

	ScoreSeed=score_seed(pos1+1, len-1, sequence, pScoring);  /* pos1+1 and len1-1, since first bases are used into the alignement */
	InitRep(pos1, pos2, len, ScoreSeed, &aligned_seed );
	diff_max = pos2-pos1;                                     /* if pos1 is bigger than pos2, then you align the chromosome with itself */
	 
	 
	 
	/*
		If the seed is already overlapping and we do not want it
	*/			

	if(opt_overlap == 0 && aligned_seed.seed_pos1+aligned_seed.seed_len >= aligned_seed.seed_pos2){

		len=(pos2-pos1);
		ScoreSeed=score_seed(pos1, len, sequence, pScoring);	
		InitRep(pos1-1, pos2-1, len+2, ScoreSeed, &aligned_seed );
		strcpy(aligned_seed.type, "Tandem");
		return aligned_seed;

	}


	/*
		Align the right side -- forward --
	*/
	pseq1 = sequence+pos1+len-1;                              /* set the sequences */
	pseq2 = sequence+pos2+len-1;
	
	
	max_col = sizeseq-(pos2+len-1)-1;                         /* you cannot align after the chromosome end */
	max_row = sizeseq-(pos1+len-1)-1;

	if(opt_overlap == 0)
		max_row = MIN2(pos2-(pos1+len-1) -1 , max_row);   /* here, you cannot go inside the next copy */
	
	
	
	align_blast2like(pseq1, pseq2, Xg, pScoring, pResults, 'f', max_row, max_col, diff_max);
	UpdateRepRight(&aligned_seed, pResults);
		
	
	/*
		Align the left side -- backward --
	*/
	pseq1 = sequence+pos1;                                     /* set the sequences */
	pseq2 = sequence+pos2;

	max_row = pos1;
	max_col = pos2;
	
	if(opt_overlap == 0)
		max_col = MIN2( aligned_seed.pos2 - (aligned_seed.pos1 + aligned_seed.len1 - 1) -2 , max_col );
	
	align_blast2like(pseq1, pseq2, Xg, pScoring, pResults, 'b', max_row, max_col, diff_max);

	UpdateRepLeft(&aligned_seed, pResults);
	
	if( aligned_seed.pos2 < aligned_seed.pos1+aligned_seed.len1 )
		strcpy(aligned_seed.type, "Overlap");
	else if( aligned_seed.pos2 == aligned_seed.pos1+aligned_seed.len1 )
		strcpy(aligned_seed.type, "Tandem");
	else if( aligned_seed.pos2 - (aligned_seed.pos1+aligned_seed.len1) < 1000 )
		strcpy(aligned_seed.type, "Close" );
	else
		strcpy(aligned_seed.type, "Distant" );

	return aligned_seed;                                         /* return the repeat */
}





/*
	Align some Inverted repeat
*/
Rep aligni(int32_t pos1, int32_t pos2, int32_t len, char *sequence, char *invsequence,\
           int32_t sizeseq ,float Xg, SCORING *pScoring, RESULTS *pResults ){

 
	Rep aligned_seed;              /* the repeat */
	char *pseq1, *pseq2;           /* pointer for the first and the second copy */
	double ScoreSeed=0;            /* score of the seed */
	int32_t max_row=0, max_col=0;  /* the dist max from the edge of the sequence max for rwo and col */
  	int32_t diff_max =0;           /* the maximum difference between row and col to 
	                                  avoid aligning a palindrome with itself. !! For Inverted it takes NEGATIVE Values !! */

	/*
		Init by th seed
	*/			
	ScoreSeed=score_seed(pos1+1, len-1, sequence, pScoring);       /* pos1+1 and len1-1, since first bases are used into the alignement */
	InitRep(pos1, pos2, len, ScoreSeed, &aligned_seed );

	if(pos1 == pos2){
		aligned_seed = Palindome(pos1, pos2, len,  sequence, pScoring);
		return aligned_seed;
	}
	/*
		Align the right side -- forward for the first copy --
	*/
	pseq1 = sequence+pos1+len-1;                                   /* set the sequences */
	pseq2 = invsequence + sizeseq - (pos2+1);
	diff_max = -(pos2-pos1-len);                                   /* for inverted repeats, we do not want palindromes
	                                                                 so put it in negative is a code to say 'inverted' */
	
	max_row = sizeseq-(pos1+len-1)-1;
	max_col = pos2;
	
	align_blast2like(pseq1, pseq2, Xg, pScoring, pResults, 'f', max_row, max_col, diff_max);
	UpdateRepRightInv(&aligned_seed, pResults);
				
	/*
		Align the left side -- backward for the first copy --
	*/
	pseq1 = sequence+pos1;                                         /* set the sequences */
	pseq2 = invsequence+sizeseq-(pos2+len-1) -1;

	max_row = pos1;
	max_col = sizeseq-(pos2+len-1)-1;

	diff_max = sizeseq;                                            /* for inverted repeats, on left, it does not really matter */
		

	align_blast2like(pseq1, pseq2, Xg, pScoring, pResults, 'b', max_row, max_col, diff_max);
	UpdateRepLeftInv(&aligned_seed, pResults);

	if(aligned_seed.pos2 == (aligned_seed.pos1+aligned_seed.len1))
		strcpy(aligned_seed.type, "Palindrome");
	else if(aligned_seed.pos2 - (aligned_seed.pos1+aligned_seed.len1) < 1000 )
		strcpy(aligned_seed.type, "Close");
	else
		strcpy(aligned_seed.type, "Distant" );


	return aligned_seed;                                             /* return the repeat */
}





/*
	Align repeats shared by 2 sequences
*/
Rep alignd_2seq(int32_t pos1, int32_t pos2, int32_t len, char *seq1, char *seq2, \
                int32_t size1, int32_t size2, float Xg, SCORING *pScoring, RESULTS *pResults ){

 
	Rep aligned_seed;              /* the repeat */
	char *pseq1, *pseq2;           /* pointer for the first and the second copy */
	double ScoreSeed=0;            /* score of the seed */
	int32_t max_row=0, max_col=0;  /* the dist max from the edge of the sequence max for rwo and col */
	int32_t diff_max =0;           /* the maximum difference between row and col to avoid aligning the chr with itself */


	/*
		Init by the seed
	*/			

	ScoreSeed=score_seed(pos1+1, len-1, seq1, pScoring);  /* pos1+1 and len1-1, since first bases are used into the alignement */
	InitRep(pos1, pos2, len, ScoreSeed, &aligned_seed );
	diff_max = size1+size2;                               /* there is not limit here */
	 

	/*
		Align the right side -- forward --
	*/
	pseq1 = seq1+pos1+len-1;                              /* set the sequences */
	pseq2 = seq2+pos2+len-1;
	
	
	max_row = size1-(pos1+len-1)-1;
	max_col = size2-(pos2+len-1)-1;           /* you cannot further than a chronmosome end */

	align_blast2like(pseq1, pseq2, Xg, pScoring, pResults, 'f', max_row, max_col, diff_max);

	UpdateRepRight(&aligned_seed, pResults);
	
	
	
	/*
		Align the left side -- backward --
	*/
	pseq1 = seq1+pos1;                                     /* set the sequences */
	pseq2 = seq2+pos2;

	max_row = pos1;
	max_col = pos2; 	                 /* you cannot further than a chronmosome begining */
	
	align_blast2like(pseq1, pseq2, Xg, pScoring, pResults, 'b', max_row, max_col, diff_max);
	UpdateRepLeft(&aligned_seed, pResults);
	
	strcpy(aligned_seed.type, "Inter");

	return aligned_seed;                                         /* return the repeat */
}



/*
	Align repeats shared by 2 sequences
*/
Rep aligni_2seq(int32_t pos1, int32_t pos2, int32_t len, char *seq1, char *invseq2, \
                int32_t size1, int32_t size2, float Xg, SCORING *pScoring, RESULTS *pResults ){


	Rep aligned_seed;              /* the repeat */
	char *pseq1, *pseq2;           /* pointer for the first and the second copy */
	double ScoreSeed=0;            /* score of the seed */
	int32_t max_row=0, max_col=0;  /* the dist max from the edge of the sequence max for rwo and col */
  	int32_t diff_max =0;           /* the maximum difference between row and col to 
	                                  avoid aligning a palindrome with itself. !! For Inverted it takes NEGATIVE Values !! */

	/*
		Init by th seed
	*/			
	ScoreSeed=score_seed(pos1+1, len-1, seq1, pScoring);       /* pos1+1 and len1-1, since first bases are used into the alignement */
	InitRep(pos1, pos2, len, ScoreSeed, &aligned_seed );



	/*
		Align the right side -- forward for the first copy --
	*/
	pseq1 = seq1+pos1+len-1;                                   /* set the sequences */
	pseq2 = invseq2 + size2 - (pos2+1);
	diff_max = size1+size2;                                    /* NO diff max for inter-repeats */
	
	max_row = size1-(pos1+len-1)-1;
	max_col = pos2;
	
	align_blast2like(pseq1, pseq2, Xg, pScoring, pResults, 'f', max_row, max_col, diff_max);
	UpdateRepRightInv(&aligned_seed, pResults);
				
	/*
		Align the left side -- backward for the first copy --
	*/
	pseq1 = seq1+pos1;                                         /* set the sequences */
	pseq2 = invseq2+size2-(pos2+len-1) -1;

	max_row = pos1;
	max_col = size2-(pos2+len-1)-1;
	
	diff_max = size1+size2;                                     /* it does not really matter */
		

	align_blast2like(pseq1, pseq2, Xg, pScoring, pResults, 'b', max_row, max_col, diff_max);
	UpdateRepLeftInv(&aligned_seed, pResults);

	if(aligned_seed.pos2 == (aligned_seed.pos1+aligned_seed.len1))
		strcpy(aligned_seed.type, "Inter");

	return aligned_seed;                                             /* return the repeat */
}                             /* return the repeat */
