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
        file     : families.c
        function : set the R-level of the repeats (Mean, Main, and fraction of Main)
                                                 
        created  : 15 Oct 2002
        modif    : July 2003, Oct 2003, Spe 05

		  note    : 	from Air.c - Mar 28 2001 / Jul 10 2001 / Jun 2 2000	

        author   : amikezor
*****/

#include "repseek_types.h"
#include "memory.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/***
	construction of the table "air" (of repeats)
***/

static void table_R(int32_t *table_r, Rep *AllReps, int32_t nRep, int32_t size_chr) 
{

	int32_t u;   /* a counter */

	for(u=0; u<nRep; u++){
      
		if(AllReps[u].pos1 == -1){continue;}  /* skip "bad repeat" */
		
		table_r[ AllReps[u].pos1 ]+=1;
		if(AllReps[u].pos1 + AllReps[u].len1 < size_chr)
			table_r[ AllReps[u].pos1 + AllReps[u].len1 ]-=1;
      
		table_r[ AllReps[u].pos2 ]+=1;
		if(AllReps[u].pos2 + AllReps[u].len2 < size_chr)
			table_r[ AllReps[u].pos2 + AllReps[u].len2 ]-=1;
	}

}

/***
	Built up and returned a memory adress where is table R
***/
static int32_t *mem_RepTableR(Repeats AllRepeats, int8_t opt_dir, int8_t opt_inv, int32_t size_chr)
{

	int32_t *table_r=NULL;  /* the table R */
	int32_t count=0;       /* allow to count increase and decrease */
	int32_t i=0;           /* just a counter */

	table_r = MyCalloc(size_chr, sizeof(int32_t), "SetRepFamily: calloc failed, bye");

	if(opt_dir)
		table_R(table_r, AllRepeats.DirRep, AllRepeats.nDirRep, size_chr);

	if(opt_inv)
		table_R(table_r, AllRepeats.InvRep, AllRepeats.nInvRep, size_chr);

	for(count=1, i=0; i<size_chr ; i++)
      table_r[i] = (count+=table_r[i]);   

	return table_r;

}

/***
	Check all repeat and set their family size by 
	looking at the table R of each copy
***/
static void Tag_Rep(Rep *AllReps, int32_t nReps, int32_t *table_r, int32_t *HistoR, int32_t Max_table_r)
{

	Rep *pRep;
	
	int32_t i;
	
	int32_t mainR=0;
	int32_t max_mainR=0;
	float meanR=0;
	float fraction_mainR=0.0;
			
	for(pRep = AllReps; pRep < AllReps+nReps; pRep++){                              /* determine the mean repeat level of each repeat */

		memset(HistoR, 0, (Max_table_r+1)*sizeof(int32_t) );
	
		for(i=0; i < pRep->len1 ; i++)
			HistoR[ table_r[pRep->pos1+i] ]++;

		for(i=0; i < pRep->len2 ;i++ )
			HistoR[ table_r[pRep->pos2+i] ]++;
		

		for(i=2, meanR=0, max_mainR=0; i <= Max_table_r ; i++){
				
			if( HistoR[i] > max_mainR ){
				max_mainR=HistoR[i];
				mainR=i;
			}
			
			meanR += (float) i * HistoR[i];
		}
			

		fraction_mainR = (float)max_mainR / (float)(pRep->len1+pRep->len2);
	
		pRep->meanR= meanR / (float)(pRep->len1+pRep->len2);
		pRep->mainR=mainR;
		pRep->fraction_mainR=fraction_mainR;
	}


}

/***
	Set the size of each family (duplicate, triplicate, etc...)
	Moreover return the table_r of the chromosome.
	In that table, each chr pos correspond to the level of repetition
	(1(unique), 2(duplicated), 3(triplicated), ...)
***/
int32_t *SetRepFamily(Repeats AllRepeats, int8_t opt_dir, int8_t opt_inv, int32_t size_chr)
{

	int32_t *table_r;         /* the tableR itself */
	int32_t max_TableR=0;    /* max value in TableR */
	int32_t *HistoTableR;    /* histo of R-level for a repeat */
	int32_t i;               /* a dummy token */

	if( !(AllRepeats.nDirRep-AllRepeats.nDirBadRep) && 
	    !(AllRepeats.nInvRep-AllRepeats.nInvBadRep) )
		 	return NULL;

	table_r = mem_RepTableR( AllRepeats, opt_dir,opt_inv,size_chr);
	
	for(i=0; i<size_chr; i++)
		max_TableR=(table_r[i]>max_TableR)?table_r[i]:max_TableR;
	
	HistoTableR = (int32_t *)MyMalloc( (max_TableR+1) * sizeof(int32_t), "SetRepFamily: malloc error" );
	
	if(opt_dir)
		Tag_Rep(AllRepeats.DirRep, AllRepeats.nDirRep, table_r, HistoTableR, max_TableR);

	if(opt_inv)
		Tag_Rep(AllRepeats.InvRep, AllRepeats.nInvRep, table_r, HistoTableR, max_TableR);
	
	MyFree(HistoTableR, (max_TableR+1)*sizeof(int32_t) );
	
	return table_r;
}
