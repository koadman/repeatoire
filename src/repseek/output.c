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
/***
        file     : output.c
        function : write out all information
                                                 
        created  : 02 Oct 2002
        modif    : July 2003 ,  Oct 2003, July04 (score is now a double), dec 2005

        author   : amikezor
***/

#include "repseek_types.h"
#include "KMRK_Seeds.h"
#include "align.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>




/***
	Write out all repeats
***/

#define FPF fprintf(fout,

void WriteRep(FILE *fout,  Repeats AllRepeats, int8_t opt_dir, int8_t opt_inv,  char opt_shape, int32_t sizeofchr)
{
	
	int32_t i;
	int32_t spacer;
	
	if(!opt_dir && !opt_inv)
		fprintf(stderr, "WriteRep: need opt_dir or opt_inv\n"),exit(4);

	if(opt_shape!='l' && opt_shape!='c')
		fprintf(stderr, "WriteRep: need opt_shape set to 'l' or 'c'\n"),exit(4);

	if(opt_dir){
	
		#define REP(i) AllRepeats.DirRep[i]

		for(i=0; i<AllRepeats.nDirRep; i++){
			
			if(REP(i).pos1 == -1)continue;
			
			if(opt_shape == 'l')
				spacer = REP(i).pos2 -1 - (REP(i).pos1+REP(i).len1-1);
			else
				spacer=MIN( REP(i).pos2-1 - (REP(i).pos1+REP(i).len1-1) , sizeofchr-1 - (REP(i).pos2+REP(i).len2-1) + REP(i).pos1 );
			
			FPF "%s.dir\t%ld\t%ld\t%ld\t%ld\t%ld\t",REP(i).type, (long) REP(i).pos1+1, (long) REP(i).pos2+1, (long) REP(i).len1, (long) REP(i).len2, (long) spacer);
			FPF "%ld-%ld-%ld-%.2f\t", (long) REP(i).seed_pos1+1, (long) REP(i).seed_pos2+1, (long) REP(i).seed_len,  REP(i).seed_meanR);
			FPF "%.3f\t%.2f\t", 100*(float)REP(i).match/(float)REP(i).align_len, REP(i).score );
			FPF "%.2f\t%ld\t%.2f\n", REP(i).meanR, (long) REP(i).mainR, REP(i).fraction_mainR );
			
		}
		
		#undef REP
	}
	
	if(opt_inv){
	
		#define REP(i) AllRepeats.InvRep[i]

		for(i=0; i<AllRepeats.nInvRep; i++){	

			if(REP(i).pos1 == -1)continue;

			if(opt_shape == 'l')
				spacer = REP(i).pos2 -1 - (REP(i).pos1+REP(i).len1-1);
			else
				spacer=MIN( REP(i).pos2-1 - (REP(i).pos1+REP(i).len1-1) , sizeofchr-1 - (REP(i).pos2+REP(i).len2-1) + REP(i).pos1 );


			FPF "%s.inv\t%ld\t%ld\t%ld\t%ld\t%ld\t",REP(i).type,(long) REP(i).pos1+1, (long) REP(i).pos2+1, (long) REP(i).len1, (long) REP(i).len2,(long) spacer);
			FPF "%ld-%ld-%ld-%.2f\t", (long)REP(i).seed_pos1+1,  (long)REP(i).seed_pos2+1,  (long)REP(i).seed_len,  REP(i).seed_meanR);
			FPF "%.3f\t%.2f\t", 100*(float)REP(i).match/(float)REP(i).align_len, REP(i).score );
			FPF "%.2f\t%ld\t%.2f\n", REP(i).meanR, (long) REP(i).mainR, REP(i).fraction_mainR );
		}
		
		#undef REP
	}
}

/***
	Write out all repeats
**/
void WriteRep_2seqs(FILE *fout, Repeats AllRepeats, int8_t opt_dir, int8_t opt_inv, int32_t size1, int32_t size2)
{
	
	int32_t i;
	int32_t diff2edge, diff1, diff2;
	
	
	if(!opt_dir && !opt_inv)
		fprintf(stderr, "WriteRep2: need opt_dir or opt_inv\n"),exit(4);

	if(opt_dir){
	
		#define REP(i) AllRepeats.DirRep[i]

		for(i=0; i<AllRepeats.nDirRep; i++){
			
			if(REP(i).pos1 == -1)continue;
			
			diff1 = ( REP(i).pos1+1>size1/2 )? REP(i).pos1+1-size1/2 : REP(i).pos1+1;
			diff2 = ( REP(i).pos2> size2/2 )? REP(i).pos2 -size2/2  : REP(i).pos2;
			diff2edge = ABS(diff1 - diff2);
			
			FPF "%s.dir\t%ld\t%ld\t%ld\t%ld\t%ld\t","Interseq",(long) REP(i).pos1+1, (long) REP(i).pos2+1, (long) REP(i).len1, (long) REP(i).len2, (long)diff2edge);
			FPF "%ld-%ld-%ld\t%.3f\t%.2f\t",(long)REP(i).seed_pos1+1, (long) REP(i).seed_pos2+1, (long) REP(i).seed_len,
			                           100*(float)REP(i).match / (float)REP(i).align_len,  REP(i).score );
			FPF "%.2f\t%ld\t%.2f\n", REP(i).meanR, (long) REP(i).mainR, REP(i).fraction_mainR );
		}
		
		#undef REP
	}
	
	
	if(opt_inv){
	
		#define REP(i) AllRepeats.InvRep[i]

		for(i=0; i<AllRepeats.nInvRep; i++){	

			if(REP(i).pos1 == -1)continue;

			diff1 = (REP(i).pos1+1>size1/2)?REP(i).pos1+1-size1/2:REP(i).pos1+1;
			diff2 = (REP(i).pos2> size2/2)?REP(i).pos2-size2/2:REP(i).pos2;
			diff2edge = ABS(diff1 - diff2);


			FPF "%s.inv\t%ld\t%ld\t%ld\t%ld\t%ld\t","Interseq", (long)REP(i).pos1+1, (long) REP(i).pos2+1 , (long) REP(i).len1, (long) REP(i).len2, (long) diff2edge);
			FPF "%ld-%ld-%ld\t%.3f\t%.2f\t", (long) REP(i).seed_pos1+1, (long) REP(i).seed_pos2+1, (long) REP(i).seed_len,
			                           100*(float)REP(i).match / (float)REP(i).align_len,  REP(i).score);
			FPF "%.2f\t%ld\t%.2f\n", REP(i).meanR, (long) REP(i).mainR, REP(i).fraction_mainR );
		}
		
		#undef REP
	}

}

/***
	Write out the table_r (degree of repetition of each seq position)
***/

void WriteTableR(FILE *fout, int32_t *table_r, int32_t max)
{

	int32_t i;

	if( table_r == NULL )
	  return;

	for(i=0;i<max; i++)
		if(table_r[i]>=2)
			fprintf(fout,
				"TableR %ld %ld\n",
				(long) i+1, 
				(long) table_r[i]);
}

/***
	Write out the seeds that will be used for extension phase
***/

void WriteSeeds(FILE *fout,  
		     AllSeeds_type *AllSeeds, 
		     int8_t opt_dir, 
		     int8_t opt_inv)
{
  Seed_type *curSeed;
	
	int32_t i;                    /* just a dummy counter */

	if(!opt_dir && !opt_inv)
		fprintf(stderr, "WriteSeeds: need opt_dir or opt_inv\n"),exit(4);

	if(opt_dir){
		for(i=0,curSeed = AllSeeds->dirSeeds; i<AllSeeds->nDirSeeds; i++,curSeed++)
		  FPF "d\t%ld\t%ld\t%ld\t%.2f\n", (long) curSeed->pos1 + 1, (long) curSeed->pos2 + 1, (long) curSeed->length, curSeed->rmean);
	}
	
	if(opt_inv){
		for(i=0,curSeed = AllSeeds->invSeeds;i<AllSeeds->nInvSeeds;i++,curSeed++)
		 FPF "i\t%ld\t%ld\t%ld\t%.2f\n", (long) curSeed->pos1 + 1, (long) curSeed->pos2 + 1, (long) curSeed->length, curSeed->rmean);
	}
}


#undef FPF


