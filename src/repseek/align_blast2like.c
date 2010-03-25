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
        file     : align_blast2like.c
        function : align 2 sequence by blast2 like algo
		             first 0,0  is a MATCH
                                                 
        created  : 07 jul 2003
        modified : July 2003, July 2004: turn int_32 scores into doubles;
		   Aug 2004: Update the BestScore when it is ">=" (instead of ">")
		   Sept 2005: store the number of steps
		   Jan 2006: correct a bug in memory reading
		   Mar 2006: correct a bug whenever only 1 cell was allowed + change align_blast2like (use max_row and maw_col instead of dist_max)
		 
	author   : amikezor
*****/

#include "repseek_types.h"
#include "memory_align.h"
#include "align.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
	Determine three scores (if any available)
*/
static int8_t IsScoreLeft(int32_t row, int32_t col, RESULTS pResults)
{
	if( col>pResults.Segment_begin[row] )
		return 1;
	else
		return 0;
}
static int8_t IsScoreAbove(int32_t row, int32_t col, RESULTS pResults)
{
	if( row>0 && col >= pResults.Segment_begin[row-1] && col <= pResults.Segment_end[row-1] )
		return 1;
	else
		return 0;
}
static int8_t IsScoreDiag(int32_t row, int32_t col, RESULTS pResults)
{
	if( row>0 && col > pResults.Segment_begin[row-1] && col <= pResults.Segment_end[row-1]+1 )
		return 1;
	else
		return 0;
}



/*
	In case of update bestscore infos
*/
static void update_Bestscore(int32_t row, int32_t col, RESULTS *pResults){
	pResults->pBestScore = pResults->pscore;
	pResults->BestScore_row = row;
	pResults->BestScore_col = col;
}




/*
	The fill 'matrix' function - even if there is no more matrix here
*/
#define BEGIN pResults->Segment_begin
#define END pResults->Segment_end
#define CURRENT pResults->pscore

#define nABOVE pResults->Top
#define nLEFT pResults->Left

#define SCORE_FROM_LEFT( X ) (*(CURRENT - X))
#define SCORE_FROM_ABOVE( row, col ) (*(CURRENT - (col-BEGIN[row]) - (END[row-1]-col+1) ) )
#define SCORE_FROM_DIAG( row, col ) (*(CURRENT - (col-BEGIN[row]) - (END[row-1]-col+1) -1 ) )

#define MATCH_F(A,B) ( pScoring->matrix[ CHAR2SYMB(seq1[ (A) ]) ][ CHAR2SYMB(seq2[ (B) ]) ])
#define MATCH_B(A,B) ( pScoring->matrix[ CHAR2SYMB( *(seq1-(A)) ) ][ CHAR2SYMB( *(seq2-(B)) ) ] )
#define MATCH(A,B) ( (direction=='f')? MATCH_F(A,B):MATCH_B(A,B)  )


/*
	In this version,
	Sequence 2 is on TOP and Sequence 1 on the LEFT
	In the trace matrix, a 0 means diagonal,
	                       #>0 means # steps from above and
	                       #<0 translates into # steps from left
*/

#define SCORE_TOP pResults->F

static void fill_matrices(char *seq1, char *seq2, SCORING *pScoring, RESULTS *pResults,
                          float Xg, int32_t max_col, int32_t max_row, int32_t diff_max, char direction)
{

	int32_t col, row;          /* current row and col */
	
	double Score_left=0.0;             /* score coming from left */
	double Score_diag=0.0;             /* score coming from diagonal */
	
	
	int8_t IsDiag,             /* markers to check if the cell in Diag, Above or Left does exist ? */
	      IsAbove,
	       IsLeft;
		
	
	CURRENT=pResults->scores;  /* set the current score at the begining of the scores */
	row=col=0;	           /* set the col and row */
	BEGIN[0]=END[0]=0;         /* initialize the arrays */

	
	/*
		Set the first Score - necessarily a match - the last/first bp of the seed
	*/
	*(CURRENT)=MATCH(0,0);
	pResults->traces[ 0 ] = 0;
	update_Bestscore(row, col, pResults);
	nABOVE[0]=0;
	SCORE_TOP[0]=*(CURRENT);
	                                                             /* printf("%c --> [%d,%d]:%f ",direction, row,col,*CURRENT); */
	/*
		First line
	*/
	col=1;
	CURRENT=pResults->scores+1;
	
	while( *(CURRENT-1) > *pResults->pBestScore-Xg && col<max_col ){
		
		if( CURRENT - pResults->scores >= pResults->max_scores)  /* Check Mem and if needed realloc */
			ReallocScores(pResults, 0);
		if( col >= pResults->max_col )
			ReallocAbove(pResults, 0);


		if(col==1) 	*CURRENT = *(CURRENT-1) + pScoring->gap_open;
		else  		*CURRENT = *(CURRENT-1) + pScoring->gap_ext;

		pResults->traces[ (CURRENT-pResults->scores) ] = -1;

		nABOVE[col]=0;                                              /* set the nABOVE[col] to 0 */
		SCORE_TOP[col]=*(CURRENT);
		                             /*fprintf(stderr,"[%d,%d]:%f vs. %f (col-1: %f) / Trace=%d\n", row,col,*CURRENT, *pResults->pBestScore-Xg, *(CURRENT-1),pResults->traces[ (CURRENT-pResults->scores) ]); */
		col++;                                                     /* increase col and current score */
		CURRENT++;

		if(direction == 'b' && col > diff_max )break;

	}


	END[row]=col-1;                                                /* ignore the last score (inferior to the limit) */
					                                               /* the last CURRENT value will be erased on next value */
		                                                           /* printf("Row %d was from %d to %d\n",row, BEGIN[row], END[row]); */

	/*
		Then fill the whole 'matrix' - the matrix being an array of stripes
	*/
	while( BEGIN[row] < END[row]  && row < max_row ){
	
		row++;
		
		if(row >= pResults->max_row)
			ReallocBegEnd(pResults, 0);
		
		END[row]=BEGIN[row]=BEGIN[row-1];
		col=BEGIN[row];
		
		
		/*
			lets go for the row/segment
		*/
		while( col <= max_col ){
		
				if( diff_max>0 && ((direction=='f' && row-col >= diff_max) ||        /* do not want chromosome aligned with itself */
				                   (direction=='b' && col-row >= diff_max))     )
					break;
		
				if(diff_max<0 && col+row > -(diff_max) )                             /* do not want plaindromes aligned with itself */
					break;

				
				if( CURRENT - pResults->scores >= pResults->max_scores)              /* Check Mem and if needed realloc */
					ReallocScores(pResults, 0); 
					
				if( col >= pResults->max_col)                                      
					ReallocAbove(pResults, 0);
					
				IsDiag=IsScoreDiag(row,  col, *pResults);                            /* Look for each potential score */
				IsAbove=IsScoreAbove(row,  col, *pResults);
				IsLeft=IsScoreLeft(row,  col, *pResults);
				
				
				/*
			   		Get the three scores (if they exist)
					as well as the three 
				*/
				if( IsLeft ){
					
					if( nLEFT>0 && Score_left+pScoring->gap_ext > SCORE_FROM_LEFT(1)+pScoring->gap_open ){
						Score_left += pScoring->gap_ext;
						nLEFT++;
					}
					else{
						Score_left = SCORE_FROM_LEFT(1) + pScoring->gap_open;
						nLEFT=1;
					}
				}
				else
					nLEFT = 0;
				
				
				if( IsAbove ){
					
					if( nABOVE[col]>0 && SCORE_TOP[col]+pScoring->gap_ext > SCORE_FROM_ABOVE(row,col)+ pScoring->gap_open){
						SCORE_TOP[col] += pScoring->gap_ext;
						nABOVE[col]++;
					}
					else{
						SCORE_TOP[col] = SCORE_FROM_ABOVE(row,col) + pScoring->gap_open;
						nABOVE[col]=1;
					}

				}
				else
					nABOVE[col] = 0;

				if( IsDiag )
					Score_diag=SCORE_FROM_DIAG(row,col)+MATCH(row,col);
					
				/*
					Compare the existing Scores
				*/
				if( IsDiag &&                                                                            /* The SCORE_DIAG is defined ..., then */
					 (
					    (!IsLeft && !IsAbove ) ||                                                        /* no other scores are defined */
						 (IsLeft  && !IsAbove && Score_diag>=Score_left ) ||                             /* bigger than Score_left */
						 (!IsLeft  && IsAbove && Score_diag>=SCORE_TOP[col] ) ||                         /* bigger than SCORE_TOP */ 
						 (IsLeft && IsAbove && Score_diag>=Score_left  && Score_diag>=SCORE_TOP[col] )   /* bigger than both */
					 )
				 ){
					*CURRENT = Score_diag;                                         /* diag is eq or sup to others */
					pResults->traces[ CURRENT-pResults->scores ] = 0;
				}
				else if(  IsAbove &&                                               /*  The SCORE_ABOVE is defined ..., then */
						   (
								(!IsLeft) ||                                       /* no other value */
								(IsLeft && SCORE_TOP[col]>Score_left )             /* bigger than Score_left */
							)
							
				 ){
				    *CURRENT = SCORE_TOP[col];                                     /* Score_top is better than Score_left */
	 				 pResults->traces[ CURRENT-pResults->scores ] = nABOVE[col];

				}
				else{
					*CURRENT = Score_left;                                         /* unknown : we don't know: both way are equal */
	 				 pResults->traces[ CURRENT-pResults->scores ] = -nLEFT;        /* arbitrarily, choose the left one */
				}
				
				                                                         /* printf("[%d:%d]: Curr=%f vs. Best=%f ; Trace= %d\n\n",row,col,*CURRENT,  *pResults->pBestScore-Xg, pResults->traces[ CURRENT-pResults->scores ]); */
				
				/*
					Compare the current score to the limit
				*/													    
				if( *CURRENT > *pResults->pBestScore - (double)Xg ){    /* If it is bigger than the limit */
					                                                    /* fprintf(stderr,"[%d,%d]:%f;%c ",row,col,*CURRENT,  pResults->traces[ CURRENT-pResults->scores ]);  */
					if( *CURRENT >= *pResults->pBestScore )             /* keep it, and eventually set it as the Best_Score */
						update_Bestscore(row, col, pResults);
						
					END[row]=col;                                       /* if match it ends at least here */
					CURRENT++;                                          /* get the score */
				}
				else{                                                   /* see where we are */
					                                                    /* printf("%ld:%ld -> last_row= %ld\n", row, col, END[row-1]); */
					if( col == BEGIN[row] )
						BEGIN[row]=col+1;                               /* set begin at next and ... */
						
					if( col > BEGIN[row] && col <= END[row-1]){
						CURRENT++;                                      /* then look at the next score */
					}
					
					if(col > END[row-1])
						break;                                          /* then just stop that loop */
				}
			col++;
			
		}
		                                                                /* printf("Row %d was from %d to %d\n\n",row, BEGIN[row], END[row]);*/
		CURRENT-=(col-END[row]-1);
	                                                                     /* fprintf(stderr,"/// BestScore= %f\n",*pResults->pBestScore); */
	}

	CURRENT--;
	pResults->nSegment=row+1;                                            /* Just remember how many segments (rows) */

}

/*
	The traceback function
	it always read the trace matrix in this new form
*/
static void traceback( char *seq1, char *seq2, SCORING * pScoring,  RESULTS *pResults , char direction){
	
	int32_t trace=0;       /* the current trace */

	int32_t row,           /* current row and cols */
	        col;

	int32_t count=0;       /* the total number of step of traceback */

	int32_t *ptraces= pResults->traces + (pResults->pBestScore - pResults->scores);    /* where we start, a pointer to the trace arrays */

	int n;


	if(direction == 'f')
		pResults->traceback=pResults->traceback_f;
	else if(direction == 'b')
		pResults->traceback=pResults->traceback_b;
	else
		fprintf(stderr,"traceback: direction is either 'f' or 'b', bye\n"),exit(4);


	row = pResults->BestScore_row;         /* set cursor at the extension end -- ! 0,0 if no similarity -- */
	col = pResults->BestScore_col;	
 	
	pResults->matches=0;
	pResults->alignment_len=0;
	
	while( row>0 || col>0 ){


		/*
			Set trace
		*/
		trace = *ptraces;                   /* read the current trace */


		/*
		fprintf(stderr,"[%d,%d] = %ld\n", row, col, trace);
		if(row<-2 || col<-2)
			exit(1);
		*/
						
						
		if(trace == 0){                                               /* if it is diagonal */
			
			if( ( direction=='f' && *(seq1+row)==*(seq2+col)) || \
			    ( direction=='b' && *(seq1-row)==*(seq2-col)) )
				(pResults->matches)++;
						
			ptraces -= (col-BEGIN[row]) + (END[row-1]-col+1) + 1;

			row--;
			col--;

			pResults->traceback[count++]= 'd';

		}
		else if(trace < 0 ){                                         /* we are going to left */
		
			
			ptraces += trace;
			col += trace;             /* think that trace is negative */
			
			for(n=0;n<ABS(trace);n++)
				pResults->traceback[count++]= 'l';
		
		}
		else {                                                       /*  choose above by default */
			
			for(n=0 ; n<trace ; n++){
				ptraces -= (col-BEGIN[row]) + (END[row-1]-col+1);
				row--;
				pResults->traceback[count++]= 'a';
			}

		}
			

		if(count >= pResults->max_alignlen-1)
			ReallocTraceback(pResults, 0, direction);		

	
	}
	
	pResults->traceback[count++] = 'd';   /* add the last bit from the first cell to the previous */
	pResults->matches++;                  /* it has to be a match --it is the seed end-- */
	
	pResults->alignment_len = count;
	pResults->traceback[count] = 0;

}

#undef CURRENT

#undef nABOVE
#undef nLEFT

#undef BEGIN
#undef END
#undef SCORE_FROM_LEFT
#undef SCORE_FROM_ABOVE
#undef SCORE_DIAG

#undef MATCH_F
#undef MATCH_B
#undef MATCH

/*
	If diff_max is negatif, then we are in inverted and do not want palindromes
	otherwise diff_max is positive and is defined to avoid aligning the chromosome with itself
	row is seq1 and col is seq2
*/
void align_blast2like(char *seq1, char *seq2, float Xg, SCORING *pScoring, RESULTS *pResults, char direction, int32_t max_row, int32_t max_col, int32_t diff_max){


	if(direction != 'f' && direction != 'b')
		fprintf(stderr,"align_blast2like: align 'f'orward or 'b'ackward, bye\n"),exit(4);

	fill_matrices(seq1, seq2, pScoring, pResults, Xg, max_col, max_row, diff_max, direction);
	traceback( seq1, seq2, pScoring, pResults, direction);
	
}


