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
        file     : read_seeds.c
        function : read seeds from a file (instead of using KMRK)
                                                 
        created  : Dec 2005
        modif    : 

        author   : amikezor
*****/


#include "repseek_types.h"
#include "KMRK_mask.h"
#include "KMRK_Seeds.h"

#include <stdio.h>
#include <stdlib.h>

/*
	This function takes a file name
	the minimum length, if dir and/or inverted should be read
	the size of the firt sequence and the size of the second sequence (if any))
*/
AllSeeds_type *read_seeds(char *seed_file, int32_t lmin, int8_t opt_dir, int8_t opt_inv, int32_t size1, int32_t size2){


	AllSeeds_type *all_seeds=NULL;    /* a seed pointer */
	FILE *fcin=NULL;                  /* file descriptor */
	
	char direction;                   /* either 'd' or 'i' */
	int8_t dir;                       /* 1 if direct, 0 otherwise */ 
	
	int32_t begin,                    /* temporary begin, end and len of a seed - tmp is needed when swap is requiered */
	        end,                      
		len,
		tmp;
	
	int32_t nseed=0;                  /* the seed_number */
		
	
	all_seeds = KMRK_allocSeeds(all_seeds, 1, opt_dir, opt_inv);    /* get seeds */
	
	/*all_seeds = KMRK_allocSeeds(all_seeds, 1, 1, 1);*/    /* get seeds */
	/*all_seeds = KMRK_allocSeeds(all_seeds, 1, 0, 1);*/    /* get seeds */
			      	
	
	fcin = fopen(seed_file, "r");
	if(!fcin)fprintf(stderr, "error whilst trying to read seed_file %s, please check it; bye\n", seed_file),exit(2);
	
	while( fscanf(fcin, "%c%ld%ld%ld", &direction, (long *) &begin, (long *) &end, (long *) &len) == 4){
				
		dir=1;                 /* by default it is direct -- for the switch -- */
				
		switch(direction){
		
			case 'i': 
				if(begin>end){ tmp=begin; begin=end; end=tmp; }
				dir=0;

			case 'd':
				if( (dir && !opt_dir) || (!dir && !opt_inv) )break;
				if(len>=lmin)
				{
					if( begin<1 || 
					    ( !size2 && end+len-1>size2 ) ||
					    ( size2 && (end+len-1>size2 || begin+len-1>size1) )
					   )
					KMRK_pushSeed(all_seeds, begin-1, end-1, len, dir);
				}
				break;
				
			default:
				fprintf(stderr, "direction of seed at line %d is not valid, bye",nseed+1) , exit(4);
		
		}
		
		while( fgetc(fcin) != '\n');                                 /* remove other  infos at the end of line */
		
		nseed++;
	}


	fclose(fcin);
	KMRK_compactSeeds(all_seeds);

	return all_seeds;
}


