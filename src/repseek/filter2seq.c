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
 * @file   KMRK.filter.c
 * @author Eric Coissac <coissac@inrialpes.fr>
 * @date   Tue Mar  2 18:33:41 2004
 * @modif  Sept 2005
 * 
 * @brief  seed filtrage function
 * 
 * 
 */
#include <stdlib.h>
#include <stdio.h>
#include "filter.h"
#include "memory.h"
#include "KMRK.h"


static void table_R2(int32_t *table_r1,
		     int32_t *table_r2,
		     Seed_type *seeds,
		     int32_t nSeed,
		     int32_t size1,
		     int32_t size2);



static void table_R2(int32_t *table_r1,
		     int32_t *table_r2,
		     Seed_type *seeds,
		     int32_t nSeed,
		     int32_t size1,
		     int32_t size2)
{

	int32_t u;
	Seed_type* curSeed; 

	for(u=0, curSeed=seeds; u<nSeed;  u++, curSeed++)
	{
		if(curSeed->pos1 == -1)
			continue;
		
		table_r1[ curSeed->pos1 ]+=1;
		
		if( (curSeed->pos1 + curSeed->length) < size1)
			table_r1[ curSeed->pos1 + curSeed->length ]-=1;
      
		table_r2[ curSeed->pos2 ]+=1;
		
		if( (curSeed->pos2 + curSeed->length) < size2)
			table_r2[ curSeed->pos2 + curSeed->length ]-=1;
	}
}


void KMRK_SeedTableR2seq(AllSeeds_type *AllSeeds, 
			     int8_t opt_dir, 
			     int8_t opt_inv, 
			     int32_t size1,
			     int32_t size2,
			     int32_t **r1,
			     int32_t **r2)
{

	int32_t *table_r1=NULL;
	int32_t *table_r2=NULL;
	int32_t count=0;
	int32_t i=0;

	table_r1= MyCalloc( size1,sizeof(int32_t),  "mem_tableR: calloc failed, bye" );
	table_r2= MyCalloc( size2,sizeof(int32_t),  "mem_tableR: calloc failed, bye" );

	if(opt_dir &&  AllSeeds->nDirSeeds)
		table_R2(table_r1, table_r2, AllSeeds->dirSeeds, AllSeeds->nDirSeeds, size1, size2);

  if(opt_inv && AllSeeds->nInvSeeds)
    table_R2(table_r1,table_r2, AllSeeds->invSeeds, AllSeeds->nInvSeeds, size1, size2);

  for(count=1, i=0; i<size1 ; i++)
    table_r1[i] = (count+=table_r1[i]);   

  for(count=1, i=0; i<size2 ; i++)
    table_r2[i] = (count+=table_r2[i]);   

  *r1 = table_r1;
  *r2 = table_r2;

}



void BuiltFilterFamily_seeds2seq(AllSeeds_type *AllSeeds, 
			     int32_t *table_r1, 
			     int32_t *table_r2,
			     double min, 
			     double Max, 
			     int8_t opt_dir, 
			     int8_t opt_inv)
{		

  float MeanR1, MeanR2;
  float TotalR1, TotalR2;

  Seed_type *curSeed;
	
  int32_t u, i;
  
  int32_t N=0;
  int8_t  opt_filter=(!min && !Max)?0:1;
	
	if(opt_dir){
	
		for(u=0, N=0, curSeed=AllSeeds->dirSeeds;  u<AllSeeds->nDirSeeds; u++,curSeed++)
		{
		
			TotalR1 = 0.0;
			TotalR2 = 0.0;
			for(i=0; i < curSeed->length ; i++)
			{
				TotalR1 += table_r1[curSeed->pos1+i]; 
				TotalR2 += table_r2[curSeed->pos2+i];    /* add both repeats together */
			}

			MeanR1 = (float)TotalR1 / ( (float) curSeed->length);       
			MeanR2 = (float)TotalR2 / ( (float) curSeed->length);       
			curSeed->rmean = MeanR1 + MeanR2 - 2;
			
			if( opt_filter && 
		    ( ( (Max && MeanR1<=Max) || !Max)   && 
	   	   ( (min && MeanR1>=min) || !min) ) &&
		    ( ( (Max && MeanR2<=Max) || !Max)   && 
		      ( (min && MeanR2>=min) || !min) ) )
 
			{     /* && has a higher priority than || */
	  
				AllSeeds->dirSeeds[N].pos1    = curSeed->pos1;                       /* thus keep it */
				AllSeeds->dirSeeds[N].pos2    = curSeed->pos2;
				AllSeeds->dirSeeds[N].length  = curSeed->length;
				AllSeeds->dirSeeds[N].rmean   = curSeed->rmean;
				
				N++;
			}
		}
	
		if(opt_filter)
		{
			AllSeeds->nFilteredDirSeeds = AllSeeds->nDirSeeds - N;
			AllSeeds->nDirSeeds=N;
	 	}

  }
	
	
	if(opt_inv){
	
		for(u=0, N=0, curSeed=AllSeeds->invSeeds; u < AllSeeds->nInvSeeds; u++,curSeed++)
		{		
			TotalR1 = 0.0;
			TotalR2 = 0.0;
	
			for(i=0; i < curSeed->length ; i++)
			{
				TotalR1 += table_r1[curSeed->pos1+i]; 
				TotalR2 += table_r2[curSeed->pos2+i];                      /* add both repeats together */
			}
			
			MeanR1 = (float)TotalR1 / ( (float) curSeed->length);       
			MeanR2 = (float)TotalR2 / ( (float) curSeed->length);       
			curSeed->rmean = MeanR1 + MeanR2 -2;

			if( opt_filter && 
					( ( (Max && MeanR1<=Max) || !Max)   && 
					( (min && MeanR1>=min) || !min) ) &&
					( ( (Max && MeanR2<=Max) || !Max)   && 
					( (min && MeanR2>=min) || !min) ) )
		  {
				
				AllSeeds->invSeeds[N].pos1   = curSeed->pos1;             /* thus keep it */
				AllSeeds->invSeeds[N].pos2   = curSeed->pos2;
				AllSeeds->invSeeds[N].length = curSeed->length;
				AllSeeds->invSeeds[N].rmean  = curSeed->rmean;	
				N++;
			}
		}

		if(opt_filter)
		{
			AllSeeds->nFilteredInvSeeds = AllSeeds->nInvSeeds-N;	      /* reset the memory */
			AllSeeds->nInvSeeds=N;
		}
	}
	
	if(opt_filter)
		KMRK_compactSeeds(AllSeeds);

  return;
}

void KMRK_FamilySeeds2seq(AllSeeds_type *AllSeeds, 
		      double min, 
		      double Max, 
		      int8_t opt_dir, 
		      int8_t opt_inv, 
		      int32_t size1,
		      int32_t size2)
{

	int32_t *table_r1;  /* the seed table R */
	int32_t *table_r2;  /* the seed table R */

	if( !AllSeeds->nDirSeeds && !AllSeeds->nInvSeeds )
	  return;
			
	if(!opt_dir && !opt_inv)
		fprintf(stderr, "FilterSeeds: opt_inv or opt_dir must be set to 1, bye\n"),exit(1);

	KMRK_SeedTableR2seq( AllSeeds, 
			     opt_dir, opt_inv, 
			     size1,size2,
			     &table_r1,&table_r2);  /* calculate the seed table R */

	BuiltFilterFamily_seeds2seq(AllSeeds, 
			 table_r1, table_r2, 
			 min, Max, 
			 opt_dir, opt_inv);   /* filter seeds if they does not fit between min and Max */
	
	MyFree(table_r1, size1*sizeof(int32_t));
	MyFree(table_r2, size2*sizeof(int32_t));
	
	return;
}


