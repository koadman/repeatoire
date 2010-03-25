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
 * @author Guillaume Achaz <gachaz@oeb.harvard.edu>, Eric Coissac <coissac@inrialpes.fr>
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


static void table_R(int32_t *table_r, 
		    Seed_type *seeds,
		    int32_t nSeed,
		    int32_t size_chr);

static void table_R( int32_t *table_r,  Seed_type *seeds,   int32_t nSeed, int32_t size_chr) 
{

  int32_t u;             /* a dummy counter */
  Seed_type* curSeed;    /* the current seed */

  for(u=0, curSeed=seeds;  u<nSeed;  u++, curSeed++){
      	
		if(curSeed->pos1 == -1)
			continue;
		
		table_r[ curSeed->pos1 ]+=1;
		
		if( (curSeed->pos1 + curSeed->length) < size_chr)
			table_r[ curSeed->pos1 + curSeed->length ]-=1;
      
      table_r[ curSeed->pos2 ]+=1;
		
      if( (curSeed->pos2 + curSeed->length) < size_chr)
			table_r[ curSeed->pos2 + curSeed->length ]-=1;
      
    }
}


int32_t *KMRK_SeedTableR(AllSeeds_type *AllSeeds, 
			int8_t opt_dir, 
			int8_t opt_inv, 
			int32_t size_chr)
{

  int32_t *table_r=NULL;    /* init the so called R-table to NULL */
  int32_t count=0;          /* a token to mark each seed */
  int32_t i=0;              /* a dummy counter */

  table_r= MyCalloc( size_chr, sizeof(int32_t), "mem_tableR: calloc failed, bye" );

  if(opt_dir &&  AllSeeds->nDirSeeds)
    table_R(table_r,
	    AllSeeds->dirSeeds,
	    AllSeeds->nDirSeeds,
	    size_chr);

  if(opt_inv && AllSeeds->nInvSeeds)
    table_R(table_r,  
	    AllSeeds->invSeeds,
	    AllSeeds->nInvSeeds,
	    size_chr);

  for(count=1, i=0; i<size_chr ; i++)
    table_r[i] = (count+=table_r[i]);   

  return table_r;

}



void BuiltFilterFamily_seeds(AllSeeds_type *AllSeeds, 
			 int32_t *table_r, 
			 double min, 
			 double Max, 
			 int8_t opt_dir, 
			 int8_t opt_inv)
{		

  double MeanR;            /* the R-mean of the seed - average redundancy */
  double TotalR;           /* the sum of all r of each site of the repeat */
  Seed_type *curSeed;      /* the current seed */
	
  int32_t u, i;
  
  int32_t N=0;
  int8_t  opt_filter=(!min && !Max)?0:1;           /* filter is on if one of the two is set by the user */
	
	if(opt_dir){
	
		for(u=0, N=0, curSeed=AllSeeds->dirSeeds;  u<AllSeeds->nDirSeeds;    u++,curSeed++)
		{
		
			TotalR = 0.0;
		
			for(i=0; i < curSeed->length ; i++)
				TotalR += (double)table_r[ curSeed->pos1+i] + table_r[curSeed->pos2+i];    /* add both repeats together */
			
			MeanR = TotalR / (2.0 * (double) curSeed->length);
			
			curSeed->rmean = MeanR;
			
			if( opt_filter && 
			   ( (Max && MeanR<=Max) || !Max)     && 
			   ( (min && MeanR>=min) || !min) )
			{    
	  
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
	
		for(u=0, N=0, curSeed=AllSeeds->invSeeds; u < AllSeeds->nInvSeeds;  u++,curSeed++)
		{		
			TotalR = 0.0;
	
			for(i=0; i < curSeed->length ; i++)
				TotalR += (double)table_r[ curSeed->pos1+i] + table_r[curSeed->pos2+i];    /* add both repeats together */
			
			MeanR = TotalR / (2.0 * (double) curSeed->length);       
			curSeed->rmean = MeanR;

			if( opt_filter && 
				( (Max && MeanR<=Max) || !Max) && 
				( (min && MeanR>=min) || !min) )
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
			AllSeeds->nFilteredInvSeeds = AllSeeds->nInvSeeds-N;	                                /* reset the memory */
			AllSeeds->nInvSeeds=N;
      	}
  }
  
	
  if(opt_filter)
    KMRK_compactSeeds(AllSeeds);

	return;
}

void KMRK_FamilySeeds(AllSeeds_type *AllSeeds, 
		      double min, 
		      double Max, 
		      int8_t opt_dir, 
		      int8_t opt_inv, 
		      int32_t size_chr)
{

	int32_t *table_r;  /* the seed table R */

	if( !AllSeeds->nDirSeeds && !AllSeeds->nInvSeeds )return;
			
	if(!opt_dir && !opt_inv)
		fprintf(stderr, "FilterSeeds: opt_inv or opt_dir must be set to 1, bye\n"),exit(1);

	table_r = KMRK_SeedTableR( AllSeeds, opt_dir, opt_inv, size_chr);  /* calculate the seed table R */

	BuiltFilterFamily_seeds(AllSeeds, table_r, min, Max, opt_dir, opt_inv);   /* filter seeds if they does not fit between min and Max */
	
	MyFree(table_r, size_chr*sizeof(int32_t));
	
	return ;
}

