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
 * @file   KMRK_merge_seeds.c
 * @author Eric Coissac <coissac@inrialpes.fr>
 * @date   Wed Mar  3 11:15:57 2004
 * 
 * @brief  Merge function of overlapping seeds
 * 
 * 
 */

#include "KMRK_merge_seeds.h"

void KMRK_MergeSeeds(AllSeeds_type *AllSeeds, 
		     int8_t opt_dir, 
		     int8_t opt_inv)
{

  int32_t i;   /* the current seed */
  int32_t j;   /* the checked seed */
		 
  int32_t N;   /* the kept ones */
  Seed_type* seeds;
	
	
  if(opt_dir){
    seeds = AllSeeds->dirSeeds;
    for(i=0, N=0 ;i<AllSeeds->nDirSeeds; i++){
			
      if(seeds[i].pos1==-1)                      /* any seed at -1 is removed */
	continue;

			
      j=i+1;
			
      while( (seeds[j].pos1!=-1) && 
	     (seeds[i].pos1!=-1) && 
	     (j < AllSeeds->nDirSeeds)        &&
	     (seeds[j].pos1 < seeds[i].pos1+ seeds[i].length))
      {
			
	if( 
	   ((seeds[j].pos2 >= seeds[i].pos2)  &&
	    (seeds[j].pos2 <  seeds[i].pos2 + seeds[i].length)) ||        /* if the seeds are overlapping */
	   ((seeds[j].pos2  + seeds[j].length >= seeds[i].pos2) &&
	    (seeds[j].pos2  + seeds[j].length < seeds[i].pos2 + seeds[i].length))) 
	  {										
	    if(seeds[j].length <= seeds[i].length)
	      seeds[j].pos1=seeds[j].pos2=seeds[j].length=-1;   /* removed the smallest */
	    else
	      seeds[i].pos1=seeds[i].pos2=seeds[i].length=-1;
	  }
				  
	j++;
      }

      if(seeds[i].pos1 !=-1)
      {                             /* if this seed is not out, store it */
			
	seeds[N].pos1   = seeds[i].pos1;
	seeds[N].pos2   = seeds[i].pos2;
	seeds[N].length = seeds[i].length;
	N++;
      }
			
    }
		
    AllSeeds->nFilteredDirSeeds += AllSeeds->nDirSeeds-N;		
    AllSeeds->nDirSeeds=N;
				
			
  }

  if(opt_inv){
    seeds = AllSeeds->invSeeds;
	
    for(i=0, N=0 ;i<AllSeeds->nInvSeeds; i++){
			
      if(seeds[i].pos1==-1)
	continue;

			
      j=i+1;
      
      while( (seeds[j].pos1!=-1 ) && 
	     (seeds[i].pos1!=-1 ) && 
	     (j < AllSeeds->nInvSeeds) &&
	     (seeds[j].pos1 < seeds[i].pos1+seeds[i].length))
	{
			
	  if( 
	     ((seeds[j].pos2 >= seeds[i].pos2)  &&                             /* if the seeds are overlapping */
	      (seeds[j].pos2 <  seeds[i].pos2+seeds[i].length)) ||
	     ((seeds[j].pos2 + seeds[j].length >= seeds[i].pos2) &&
	      (seeds[j].pos2 + seeds[j].length <  seeds[i].pos2+seeds[i].length))) 
	    {
	      
	      if(seeds[j].length <= seeds[i].length)
		seeds[j].pos1=seeds[j].pos2=seeds[j].length=-1;   /* removed the smallest */
	      else
		seeds[i].pos1=seeds[i].pos2=seeds[i].length=-1;
	    }
				  
	  j++;
	}

      if(seeds[i].pos1!=-1)
	{                             /* if this seed is not out, store it */
			
	  seeds[N].pos1   = seeds[i].pos1;
	  seeds[N].pos2   = seeds[i].pos2;
	  seeds[N].length = seeds[i].length;
	  N++;
	}
			
		}
		
    AllSeeds->nFilteredInvSeeds += AllSeeds->nInvSeeds-N;		
    AllSeeds->nInvSeeds=N;
  }	

  KMRK_compactSeeds(AllSeeds);

}
