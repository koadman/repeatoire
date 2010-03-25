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
 * @file   align_memory.c
 * @author amikezor <gachaz@oeb.harvard.edu>
 * @date   Jul 07 2003
 * @modif  July 2003, Oct 2003, Nov 2003 -small bug in realloc traceback, Sep 05
 * 
 * @brief  all malloc/realloc mem of the results structure
 * 
 * 
 */
 

#include "repseek_types.h"
#include "memory.h"
#include <stdio.h>
#include <stdlib.h>

#define MIN_COLS   1000
#define MIN_ROWS   1000
#define MIN_ALIGN  10000
#define MIN_SCORES 10000

/*
	All the memory you need for the alignment itself
*/
void mem_Blast2Align(RESULTS *pResults){

	pResults->scores = (double *)  MyCalloc( MIN_SCORES, sizeof(double) , "mem_Blast2Align: scores calloc error");
	pResults->traces = (int32_t *) MyCalloc( MIN_SCORES, sizeof(int32_t) , "mem_Blast2Align: traces calloc error" );

	pResults->traceback_f = (char *) MyCalloc( MIN_ALIGN, sizeof(char) , "mem_Blast2Align: traceback_f calloc error");
	pResults->traceback_b = (char *) MyCalloc( MIN_ALIGN, sizeof(char)  , "mem_Blast2Align: traceback_b calloc error");
	
	pResults->Segment_begin = (int32_t *) MyCalloc( MIN_ROWS, sizeof(int32_t)  , "mem_Blast2Align: Segment_begin calloc error");
	pResults->Segment_end   = (int32_t *) MyCalloc( MIN_ROWS, sizeof(int32_t)  , "mem_Blast2Align: Segment_end calloc error");
	
	pResults->Top = (int32_t *) MyCalloc( MIN_COLS, sizeof(int32_t)  , "mem_Blast2Align: Top calloc error");
	pResults->F   =  (double *) MyCalloc( MIN_COLS, sizeof(double)  , "mem_Blast2Align: F calloc error");
	
			
	pResults->BestScore_row=0;
	pResults->BestScore_col=0;
	
	pResults->max_scores=MIN_SCORES;
	pResults->max_col= MIN_COLS;
	pResults->max_row= MIN_ROWS;
	pResults->max_alignlen=MIN_ALIGN;

}
/*
	The alignemnt results structure
*/
void free_Blast2Align(RESULTS *pResults){

	MyFree(pResults->traceback_f, pResults->max_alignlen*sizeof(char));
	MyFree(pResults->traceback_b, pResults->max_alignlen*sizeof(char));
	MyFree(pResults->scores, pResults->max_scores*sizeof(double));
	MyFree(pResults->traces, pResults->max_scores*sizeof(int32_t));

	MyFree(pResults->Top, pResults->max_col*sizeof(int32_t));
	MyFree(pResults->F, pResults->max_col*sizeof(double));

	MyFree(pResults->Segment_begin, pResults->max_row*sizeof(int32_t));
	MyFree(pResults->Segment_end, pResults->max_row*sizeof(int32_t));
	
}

/*
	if sig == 0, then double mem,
	otherwise readjust to sig
*/
void ReallocScores(RESULTS *pResults, int32_t sig){

	int32_t old=pResults->max_scores;

	int32_t nscores = pResults->pscore - pResults->scores;
	int32_t Bscores = pResults->pBestScore - pResults->scores;
	
	
	if( ! sig )
		(pResults->max_scores)*=2;
	else
		pResults->max_scores=sig;

				
	pResults->traces = (int32_t *)MyRealloc( pResults->traces , sizeof(int32_t) * pResults->max_scores, old*sizeof(int32_t),
	                                    "ReallocScores: realloc traces error, bye" );
													
	pResults->scores = (double *)MyRealloc( pResults->scores ,  sizeof(double)  * pResults->max_scores , old*sizeof(double),
	                                    "ReallocScores: realloc scores error, bye" );
													
	pResults->pscore=(pResults->scores + nscores);         /* set correctly the current score */
	pResults->pBestScore=(pResults->scores + Bscores);     /* set correctly the current best score */

}

/*
	if sig == 0, then double mem,
	otherwise readjust to sig
*/
void ReallocTraceback(RESULTS *pResults, int32_t sig, char direction){

	int32_t old=pResults->max_alignlen;

	if( ! sig )
		(pResults->max_alignlen)*=2;
	else
		pResults->max_alignlen=sig;
		
	pResults->traceback_f = (char *)MyRealloc( pResults->traceback_f ,sizeof(char)*pResults->max_alignlen, old*sizeof(char),
	                                           "ReallocTraceback: traceback_f realloc error, bye" );

	pResults->traceback_b = (char *)MyRealloc( pResults->traceback_b ,sizeof(char)*pResults->max_alignlen, old*sizeof(char),
	                                           "ReallocTraceback: traceback_b realloc error, bye" );

	if(direction == 'f')
		pResults->traceback=pResults->traceback_f;
	else if(direction == 'b')
		pResults->traceback=pResults->traceback_b;
	else
		fprintf(stderr,"traceback: direction is either 'f' or 'b', bye\n"),exit(4);
		
}


/*
	if sig == 0, then double mem,
	otherwise readjust to sig
*/
void ReallocAbove(RESULTS *pResults, int32_t sig){
	
	int32_t old=pResults->max_col;

	if( ! sig )
		(pResults->max_col)*=2;
	else
		pResults->max_col=sig;
		
	pResults->Top = (int32_t *)MyRealloc( (void *)pResults->Top ,sizeof(int32_t)*pResults->max_col, sizeof(int32_t)*old, "ReallocAbove: realloc error, bye" );
	pResults->F = (double *)MyRealloc( (void *)pResults->F ,sizeof(double)*pResults->max_col, sizeof(double)*old, "ReallocAbove: realloc error, bye" );
	
}

/*
	if sig == 0, then double mem,
	otherwise readjust to sig
*/
void ReallocBegEnd(RESULTS *pResults, int32_t sig){

	int32_t old=pResults->max_row;
	
	if( ! sig )
		(pResults->max_row)*=2;
	else
		pResults->max_row=sig;
		
	pResults->Segment_begin = (int32_t *)MyRealloc( (void *)pResults->Segment_begin , sizeof(int32_t)*pResults->max_row ,
	                                             sizeof(int32_t)*old, "ReallocBegEnd: Segment_begin realloc error, bye");
	
	pResults->Segment_end = (int32_t *)MyRealloc( (void *)pResults->Segment_end ,sizeof(int32_t)*pResults->max_row ,
	                                             sizeof(int32_t)*old, "ReallocBegEnd: Segment_end realloc error, bye");
}


#undef MIN_SCORES
#undef MIN_COLS
#undef MIN_ALIGN
