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
 * @file   families_2seq.c
 * @author amikezor
 * @date   14 Nov 2002
 * @modif  Jul 2003 - change to take into account Mean,  Main and fraction of R-level
 * @modif  April 2004 - Recode totally with 2 TableR and only 1 meanR
 * @modif  Sept 2005
 * 
 * @brief  set the number of copy in the 2 input seq
 * 
 * 
 */

#include "repseek_types.h"
#include "memory.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/***
	construction of the table "air" (of repeats)
***/

static void Fill_Tables_2seq(int32_t size1, int32_t size2, int32_t *TableR1, int32_t * TableR2, Rep *AllReps, int32_t nRep) 
{

	int32_t u;   /* a counter */

	for(u=0; u<nRep; u++){
      		
		if(AllReps[u].pos1 == -1){continue;}  /* skip "bad repeat" */
		
		TableR1[ AllReps[u].pos1 ]+=1;
		
		if(AllReps[u].pos1 + AllReps[u].len1 < size1)
			TableR1[ AllReps[u].pos1 + AllReps[u].len1 ]-=1;
      
		TableR2[ AllReps[u].pos2 ]+=1;
		
		if(AllReps[u].pos2 + AllReps[u].len2 < size2)
			TableR2[ AllReps[u].pos2 + AllReps[u].len2 ]-=1;
	}

}



/***
	Check all repeat and set their copy number by 
	looking at the table R of each copy
***/
static void Tag_Rep_2seqs(Rep *AllReps, int32_t nReps, int32_t *TableR1, int32_t *TableR2, int32_t *HistoR, int32_t Max_table_r)
{

	Rep *pRep;
	
	int32_t i;

	int32_t mainR=0;
	int32_t max_mainR=0;
	float meanR=0, meanR2=0;
	float fraction_mainR=0.0;


	for(pRep = AllReps; pRep < AllReps+nReps; pRep++){         /* determine the mean repeat level of each repeat */
	
		memset(HistoR, 0, (Max_table_r+1)*sizeof(int32_t) );

		for(i=0; i < pRep->len1 ; i++){
			HistoR[ TableR1[pRep->pos1+i] ]++;
			meanR += TableR1[pRep->pos1+i];
		}
		
		for(i=0; i < pRep->len2 ;i++ ){
			HistoR[ TableR2[pRep->pos2+i] ]++;
			meanR2+=TableR2[pRep->pos2+i];
		}

		meanR /= (float)(pRep->len1);
		meanR2 /= (float)(pRep->len2);

		for(i=1, max_mainR = 0; i <= Max_table_r ; i++){
						
			if( HistoR[i] > max_mainR ){
				max_mainR=HistoR[i];
				mainR=i;
			}
						
		}

		fraction_mainR = (float)max_mainR / ( (float)(pRep->len1) + (pRep->len2) );


		pRep->meanR= meanR+meanR2 -2;
		pRep->mainR=0;
		pRep->fraction_mainR=0;

	}


}

/***
	Set the number of copies of each sequence in the second seq (2 means that this rep has 2 copies in seq2)
	Moreover return the table_r of the chromosome.
	In that table, each chr pos correspond to the level of repetition
	(1(unique), 2(duplicated), 3(triplicated), ...)
***/
void SetNumberRep_2seqs(Repeats AllRepeats, int8_t opt_dir, int8_t opt_inv,
                        int32_t size1, int32_t size2,
								int32_t **pT1, int32_t **pT2)
{

	int32_t max_TableR=0;    /* the maximum redundancy in the TableR1 or TableR2 */
	int32_t *HistoTableR;    /* Histogram of the R-level of a repeat --> here R-level of copy 1 into sequence 2*/
	int32_t i;               /* dummy counter */
	int32_t count;           /* dummy counter 2 */
	
	int32_t *TableR1;
	int32_t *TableR2;
	
	
	*pT1 = NULL;
	*pT2 = NULL;
	
	if( !(AllRepeats.nDirRep-AllRepeats.nDirBadRep) && 
	    !(AllRepeats.nInvRep-AllRepeats.nInvBadRep) )
		 	return;

	/* 
		Alloc Tables
	*/
	TableR1 = (int32_t *)MyCalloc(size1, sizeof(int32_t), "SetRepFamily: calloc failed, bye");
	TableR2 = (int32_t *)MyCalloc(size2, sizeof(int32_t), "SetRepFamily: calloc failed, bye");
		

	/*
		Fill those Tables
	*/
	if(opt_dir) Fill_Tables_2seq(size1, size2, TableR1, TableR2, AllRepeats.DirRep, AllRepeats.nDirRep);
	if(opt_inv)	Fill_Tables_2seq(size1, size2, TableR1, TableR2, AllRepeats.InvRep, AllRepeats.nInvRep);
	

	for(count=1, i=0; i<size1 ; i++) TableR1[i] = (count += TableR1[i]);   
	for(count=1, i=0; i<size2 ; i++) TableR2[i] = (count += TableR2[i]);   

	/*
		Get Their max values
	*/
	max_TableR=0;
	for(i=0; i<size1; i++) max_TableR=( TableR1[i]>max_TableR)? TableR1[i]:max_TableR;  	 /* set the max value */
	for(i=0; i<size2; i++) max_TableR=( TableR2[i]>max_TableR)? TableR2[i]:max_TableR;  	 /* set the max value */

	/*
		Alloc only once the histogram
	*/
	HistoTableR  = (int32_t *)MyMalloc( (max_TableR+1) * sizeof(int32_t), "SetRepFamily: malloc error" );     /* get some memory */

	/*
		Tag the repeats
	*/
	if(opt_dir)
		Tag_Rep_2seqs(AllRepeats.DirRep, AllRepeats.nDirRep,  TableR1, TableR2, HistoTableR , max_TableR);

	if(opt_inv)
		Tag_Rep_2seqs(AllRepeats.InvRep, AllRepeats.nInvRep,  TableR1, TableR2, HistoTableR , max_TableR);


	/*
		Free the histogram
	*/
	MyFree(HistoTableR, (max_TableR+1) * sizeof(int32_t));

	*pT1 = TableR1;
	*pT2 = TableR2;
	
}
