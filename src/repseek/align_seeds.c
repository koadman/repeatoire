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
        file     : align_seeds.c
        function : from seeds to repeats
        created  : jul 7 2003
        modif    : July 2003, Oct 2003, Nov 2003 --change the reset_Rep(),Feb 2004
		   July 2004  : turn scores into doubles
		   Nov 2004   : correct a bug with opt_overlap --now removed too small and overlapping seeds--
		   Sep 05
		   March 06: debug and make optimization with reset_seeds()
		   March 06: change reset_repet() to let the user choose how stringent should be the overlap
		   May 07: update reset_repet() to to remove a small bug
		   Dec 09: correct a bug -- thanks to M Krawczykc 

		  author   : amikezor
*****/

#include "KMRK_Seeds.h"

#include "repseek_types.h"
#include "align.h"
#include "sequence.h"
#include "memory.h"
#include "memory_align.h"
#include "sort.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>




/*
	A function to remove all seeds that fall in the 'f'orward alignemnt of the current seed
*/
static void reset_seeds(AllSeeds_type *PtrAllSeeds, char orientation, int32_t current_seed, Rep aligned_seed, RESULTS * pResults)
{
	int32_t max_seed=0,        /* max of inv or dir */
	      s=current_seed+1;    /* the seed counter */
	
	
	char *ptraceback=NULL;


	int32_t Current_pos1,    /* the pos in the current forward alignemnt */
	        Current_pos2;
	

	int32_t mem_pos1=0, mem_pos2=0;                         /* values used for optimization (avoid always going back to start) */
	char *mem_trace=NULL;
	
	int32_t len_traceback = strlen(pResults->traceback_f);  /* the end of the traces are the begining of the alignment */

	
	/*
		Treat separatly direct and inverted
	*/
	if(orientation == 'd'){
	
	
		max_seed = PtrAllSeeds->nDirSeeds;
		mem_pos1 = 0;                          /* before the first seed checked */

		while( s<max_seed && PtrAllSeeds->dirSeeds[s].pos1 < aligned_seed.pos1+aligned_seed.len1 ){
		
			if( PtrAllSeeds->dirSeeds[s].pos1 != -1 &&
			    PtrAllSeeds->dirSeeds[s].pos2 > aligned_seed.seed_pos2 &&
			    PtrAllSeeds->dirSeeds[s].pos2 <= aligned_seed.pos2+aligned_seed.len2 ){        /* then it is a potential candidate */
				
				if( mem_pos1 == 0 ){                               /* first round */

					Current_pos1 = aligned_seed.seed_pos1 + aligned_seed.seed_len ;   /* the first bit after the seed */
					Current_pos2 = aligned_seed.seed_pos2 + aligned_seed.seed_len ;
					ptraceback = pResults->traceback_f + len_traceback - 2;           /* at the end --last cell is the seed end-- */
					
				}else{
					Current_pos1=mem_pos1;                  /* subsequent rounds */
					Current_pos2=mem_pos2;
					ptraceback = mem_trace;
				
					mem_pos1 = 0;
				}

				
				while( ptraceback >= pResults->traceback_f &&
				       Current_pos1 <= PtrAllSeeds->dirSeeds[s].pos1 &&
				       Current_pos2 <= PtrAllSeeds->dirSeeds[s].pos2 ){
				
					switch( *ptraceback ){
						case 'd':
							Current_pos1++;
							Current_pos2++;
							break;
						case 'a':
							Current_pos1++;
							break;
						case 'l':
							Current_pos2++;
							break;
						default:
							fprintf(stderr,"reset_seeds: expect traceback to be 'd','a' or 'l', NOT %c !! bye\n", *ptraceback),exit(6);
					}

					if( Current_pos1 ==  PtrAllSeeds->dirSeeds[s].pos1 && mem_pos1 == 0 ){
						mem_pos1 = Current_pos1;
						mem_pos2 = Current_pos2;
						mem_trace = ptraceback-1;           /* ptraceback-1, is the next one */
					}
					
					if(Current_pos1 ==  PtrAllSeeds->dirSeeds[s].pos1 && Current_pos2 == PtrAllSeeds->dirSeeds[s].pos2 ){					
						
						/*
						fprintf(stderr, "Remove %d %d-%d-%d USING %d %d-%d-%d \n", 
						                 s+1, PtrAllSeeds->dirSeeds[s].pos1+1, PtrAllSeeds->dirSeeds[s].pos2+1, PtrAllSeeds->dirSeeds[s].length,
								 current_seed+1, PtrAllSeeds->dirSeeds[current_seed].pos1+1, PtrAllSeeds->dirSeeds[current_seed].pos2+1, PtrAllSeeds->dirSeeds[current_seed].length);
						
						*/
						PtrAllSeeds->dirSeeds[s].pos1 = -1;
					}
					ptraceback--;
				}
			}
			
			s++;
		}
	
	}
	else if (orientation == 'i'){
	
		max_seed = PtrAllSeeds->nInvSeeds;
		mem_pos1 = 0;
				
		while( s<max_seed && PtrAllSeeds->invSeeds[s].pos1 < aligned_seed.pos1+aligned_seed.len1 ){

			if(PtrAllSeeds->invSeeds[s].pos1 != -1  &&
			   PtrAllSeeds->invSeeds[s].pos2 >= aligned_seed.pos2  &&
			   PtrAllSeeds->invSeeds[s].pos2 <= aligned_seed.seed_pos2 ){

				if( ! mem_pos1 ){                               /* first round */

					Current_pos1 = aligned_seed.seed_pos1+aligned_seed.seed_len;  /* the first bit after the seed */
					Current_pos2 = aligned_seed.seed_pos2-1;                      /* the first bit before the seed 2nd copy */
					ptraceback=pResults->traceback_f+len_traceback-2;             /* the last cell is the seed end */
					
				}else{

					Current_pos1=mem_pos1;                  /* subsequent rounds */
					Current_pos2=mem_pos2;
					ptraceback = mem_trace;
				
					mem_pos1 = 0;
				}
								
				while( ptraceback >= pResults->traceback_f &&
				       Current_pos1 <= PtrAllSeeds->invSeeds[s].pos1 &&
				       Current_pos2 >= PtrAllSeeds->invSeeds[s].pos2 ){
				
				
					switch( *ptraceback ){
						case 'd':
							Current_pos1++;
							Current_pos2--;
							break;
						case 'a':
							Current_pos1++;
							break;
						case 'l':
							Current_pos2--;
							break;
						default:
							fprintf(stderr,"reset_seeds: expect traceback to be 'd','a' or 'l', NOT %c !! bye\n",*ptraceback),exit(6);
					}					
					
					
					if( Current_pos1 ==  PtrAllSeeds->invSeeds[s].pos1 && mem_pos1 == 0 ){
						mem_pos1 = Current_pos1;
						mem_pos2 = Current_pos2;
						mem_trace = ptraceback-1;
					}

					if( Current_pos1 ==  PtrAllSeeds->invSeeds[s].pos1 && Current_pos2 == PtrAllSeeds->invSeeds[s].pos2+PtrAllSeeds->invSeeds[s].length-1 )
						PtrAllSeeds->invSeeds[s].pos1 = -1;

					ptraceback--;

				} /* endof while( ptraceback >= pResults-> ... */
			   }

			s++;
		}

	}
	else
		fprintf(stderr,"reset_seeds: orientation has to be 'd' or 'i', bye\n");

}




/*
	Remove Repeats which are have the same pos1, pos2, len1 and len2
*/
static int32_t reset_Rep(Rep *RepTable, int32_t nRep, float merge_repeats ){

	int32_t u, i,     
	      removed=0;                     
	
	int32_t overlap_copy1,
	    overlap_copy2;
	
	if(merge_repeats>1 || merge_repeats<0)
		fprintf(stderr, "reset_Rep:  0.0 <= fraction_overlap <= 1.0\n"),exit(7);
	
	for(u=0; u<nRep-1; u++){
	
		i=u+1;
		
		if(RepTable[u].pos1 == -1)
			continue;
		
		while( i<nRep &&
		       RepTable[u].pos1 != -1 && 
		       RepTable[i].pos1 <= RepTable[u].pos1+RepTable[u].len1
		      ){

			if(RepTable[i].pos1 == -1){
				i++;
				continue;
			}
			
			
						overlap_copy1 = MIN2( (RepTable[i].pos1+RepTable[i].len1), (RepTable[u].pos1+RepTable[u].len1) ) - MAX2(RepTable[i].pos1, RepTable[u].pos1 );
			overlap_copy2 = MIN2( (RepTable[i].pos2+RepTable[i].len2), (RepTable[u].pos2+RepTable[u].len2) ) - MAX2(RepTable[i].pos2, RepTable[u].pos2 );

			if(  overlap_copy1 >= (int32_t) merge_repeats*RepTable[i].len1 && overlap_copy1 >= (int32_t) merge_repeats*RepTable[u].len1 && 
			     overlap_copy2 >= (int32_t) merge_repeats*RepTable[i].len2 && overlap_copy2 >= (int32_t) merge_repeats*RepTable[u].len2
			       
			){
			
				if(RepTable[i].score > RepTable[u].score){                                  /* keep the best score */
					
					/*
					fprintf(stderr,
					        "REMOVE_REP: %d %d %d %d %d-%d-%d by %d %d %d %d %d-%d-%d\n", 
					         RepTable[u].pos1+1, RepTable[u].pos2+1, RepTable[u].len1, RepTable[u].len2,  RepTable[u].seed_pos1+1, RepTable[u].seed_pos2+1, RepTable[u].seed_len,  
					         RepTable[i].pos1+1, RepTable[i].pos2+1, RepTable[i].len1, RepTable[i].len2,  RepTable[i].seed_pos1+1, RepTable[i].seed_pos2+1, RepTable[i].seed_len
						 );
					*/
					RepTable[u].pos1=RepTable[u].pos2=RepTable[u].len1=RepTable[u].len2=-1; 
				}
				else{
					/*
					fprintf(stderr,
					        "REMOVE_REP: %d %d %d %d %d-%d-%d by %d %d %d %d %d-%d-%d\n", 
					         RepTable[i].pos1+1, RepTable[i].pos2+1, RepTable[i].len1, RepTable[i].len2,  RepTable[i].seed_pos1+1, RepTable[i].seed_pos2+1, RepTable[i].seed_len,  
					         RepTable[u].pos1+1, RepTable[u].pos2+1, RepTable[u].len1, RepTable[u].len2,  RepTable[u].seed_pos1+1, RepTable[u].seed_pos2+1, RepTable[u].seed_len
						 );
					*/
					RepTable[i].pos1=RepTable[i].pos2=RepTable[i].len1=RepTable[i].len2=-1;
				}
				removed++;
				
			}
			i++;
		}
	}
	
	return removed;
}		


/***
	get Repeats from seeds.
	This function is the only one exported and is called by main()
***/
Repeats align_seeds(AllSeeds_type *AllSeeds, char *sequence,  float Xg, float gap_open, float gap_ext, char matrix_type,\
                    int16_t Lmin, int8_t opt_dir, int8_t opt_inv, int8_t opt_overlap, float merge_repeats ){

	
	int32_t sequence_size;     /* guess... */
 
	Repeats AllRepeats;        /* A structure that contain all Repeats */
		
	int32_t s, r=0;            /* to count the 's'eeds and 'r'epeats */

	SCORING *pScoring;           /* Scoring - get scoring matrix and gap-open, gap-ext */
        pScoring = (SCORING*) malloc(sizeof(SCORING));
        
        //fprintf(stderr,"%f\n",pScoring->matrix[25][25]);
	RESULTS Results;           /* Results all the results of the alignemnt */
	
	char *invsequence=NULL;
	
 	if(!opt_dir && !opt_inv)
		fprintf(stderr, "align_seeds: At least one of the direct and inverted is needed, bye\n"),exit(1);
		 
	AllRepeats = mem_Repeats(AllSeeds->nDirSeeds, AllSeeds->nInvSeeds);    /* Init and Mem */
	mem_Blast2Align(&Results);                                             /* get mem for alignment */
	

	if( matrix_type == 'i')
		identity_matrix(sequence,gap_open, gap_ext, pScoring);                         /* Load identity matrix into memory */
	else if(matrix_type == 'l' )
	{
		log_matrix(sequence, gap_open, gap_ext, pScoring);
        }
	sequence_size = (int32_t)strlen(sequence);
	
	if(opt_dir){
	
		for(s=0 , r=0; s < AllSeeds->nDirSeeds; s++){

			if(s%10 == 0)
				fprintf(stderr,"direct: %d seeds ->  %d repeats\r",s,r);
										
			if(!opt_overlap &&
			    AllSeeds->dirSeeds[s].pos1 + AllSeeds->dirSeeds[s].length > AllSeeds->dirSeeds[s].pos2 &&
			    AllSeeds->dirSeeds[s].pos2- AllSeeds->dirSeeds[s].pos1 < Lmin
			)
				AllSeeds->dirSeeds[s].pos1 = -1;                                         /* overlapping and too small, remove it */

			if(AllSeeds->dirSeeds[s].pos1 == -1 )continue;                                     /* thus it is a removed seed */		

			AllRepeats.DirRep[r] = alignd(AllSeeds->dirSeeds[s].pos1,
			                              AllSeeds->dirSeeds[s].pos2,
						       AllSeeds->dirSeeds[s].length,
						       sequence, sequence_size ,Xg, 
						       pScoring, &Results, opt_overlap);


			AllRepeats.DirRep[r].seed_meanR  = AllSeeds->dirSeeds[s].rmean;
			
 			reset_seeds(AllSeeds, 'd', s, AllRepeats.DirRep[r], &Results);
 
			r++;
		}

		AllRepeats.nDirRep = r;

		if(AllRepeats.nDirRep != 0){
		
			AllRepeats.DirRep = (Rep *)MyRealloc(AllRepeats.DirRep, AllRepeats.nDirRep*sizeof(Rep),
                                    	        AllSeeds->nDirSeeds*sizeof(Rep) ,"align_seeds : realloc DirRep error");  
															  			
			qsort_1then2_repeats(AllRepeats.DirRep, AllRepeats.nDirRep);

			AllRepeats.nDirBadRep = reset_Rep(AllRepeats.DirRep, AllRepeats.nDirRep, merge_repeats);
		}
		else{
			MyFree(AllRepeats.DirRep, AllSeeds->nDirSeeds*sizeof(Rep) );
			AllRepeats.DirRep=NULL;
		}
	}
	if(opt_inv){
		
		/*
			Get Inv sequence
		*/
		invsequence = (char *)MyMalloc((sequence_size+1)*sizeof(char),"align_seeds: invsequence malloc error, bye");
		invseq(sequence,invsequence);
		
		for(s=0, r=0; s < AllSeeds->nInvSeeds; s++){
		
			if( s%10 == 0)
				fprintf(stderr,"inverted: %d seeds ->  %d repeats\r",s,r);
			
			if(AllSeeds->invSeeds[s].pos1 == -1)continue;              
			                          
			AllRepeats.InvRep[r] = aligni(AllSeeds->invSeeds[s].pos1,
			                              AllSeeds->invSeeds[s].pos2,
						      AllSeeds->invSeeds[s].length,
						      sequence, invsequence,sequence_size,
						      Xg, pScoring, &Results);
			
			AllRepeats.InvRep[r].seed_meanR  = AllSeeds->invSeeds[s].rmean;

			if(AllRepeats.InvRep[r].len1 >= Lmin){                                      /* if not do not keep small palindromes */
				reset_seeds(AllSeeds, 'i', s, AllRepeats.InvRep[r],  &Results);   
				r++;	
			}
			
		}
		
		AllRepeats.nInvRep = r;
		
		if(AllRepeats.nInvRep != 0){
			AllRepeats.InvRep = (Rep *)MyRealloc(AllRepeats.InvRep, AllRepeats.nInvRep*sizeof(Rep),
			                                    AllSeeds->nInvSeeds*sizeof(Rep),"align_seeds : realloc InvRep error");
			
			qsort_1then2_repeats(AllRepeats.InvRep, AllRepeats.nInvRep);
			
			AllRepeats.nInvBadRep = reset_Rep(AllRepeats.InvRep, AllRepeats.nInvRep, merge_repeats);
		}
		else{
			MyFree( AllRepeats.InvRep, AllSeeds->nInvSeeds*sizeof(Rep) );
			AllRepeats.InvRep=NULL;
		}
		
		MyFree(invsequence, sequence_size*sizeof(char) );
		
	}
	free_Blast2Align(&Results);                       /* free the RESULTS struct */		

	return AllRepeats;
}




/***
	get Repeats from seeds.
	This function is the only one exported and is called by main()
***/
Repeats align_seeds_2seq(AllSeeds_type *AllSeeds, char *seq1, char *seq2,  float Xg, float gap_open, float gap_ext, char matrix_type,\
                    int16_t Lmin, int8_t opt_dir, int8_t opt_inv , float merge_repeats){

	
	int32_t size1, size2;      /* sequences size */
 
	Repeats AllRepeats;        /* A structure that contain all Repeats */
		
	int32_t s, r=0;            /* to count the 's'eeds and 'r'epeats */

	SCORING Scoring;           /* Scoring - get scoring matrix and gap-open, gap-ext */
	RESULTS Results;           /* Results all the results of the alignemnt */
	
	char *invseq2=NULL;        /* a pointer to the inevrted sequence 2 */
	
	
 	if(!opt_dir && !opt_inv)
		fprintf(stderr, "align_seeds: At least one of the direct and inverted, bye\n"),exit(1);


	size1=(int32_t)strlen(seq1);
	size2=(int32_t)strlen(seq2);

		 
	AllRepeats = mem_Repeats(AllSeeds->nDirSeeds, AllSeeds->nInvSeeds);  /* Init and Mem */
	mem_Blast2Align(&Results);                                           /* get mem for alignment */
	

	if( matrix_type == 'i')
		identity_matrix(seq1,gap_open, gap_ext, &Scoring);                         /* Load identity matrix into memory */
	else if(matrix_type == 'l' )
		log_matrix(seq1, gap_open, gap_ext, &Scoring);


	if(opt_dir){
		
		for(s=0 , r=0; s < AllSeeds->nDirSeeds; s++){

			if(s%10 == 0)
				fprintf(stderr,"direct: %d seeds ->  %d repeats\r",s,r);
		
			if(AllSeeds->dirSeeds[s].pos1 == -1 )continue;                        /* thus it is a removed seed */		

			AllRepeats.DirRep[r] = alignd_2seq(AllSeeds->dirSeeds[s].pos1,
			                                 AllSeeds->dirSeeds[s].pos2,
													   AllSeeds->dirSeeds[s].length,
													   seq1, seq2,
													   size1, size2,
													   Xg, &Scoring, &Results);
													

			AllRepeats.DirRep[r].seed_meanR  = AllSeeds->dirSeeds[s].rmean;
			
			reset_seeds(AllSeeds, 'd', s, AllRepeats.DirRep[r], &Results);

			r++;
		}

		AllRepeats.nDirRep = r;

		if(AllRepeats.nDirRep != 0){
		
			AllRepeats.DirRep = (Rep *)MyRealloc(AllRepeats.DirRep, AllRepeats.nDirRep*sizeof(Rep),
                                    	        AllSeeds->nDirSeeds*sizeof(Rep) ,"align_seeds : realloc DirRep error");  
															  			
			qsort_1then2_repeats(AllRepeats.DirRep, AllRepeats.nDirRep);

			AllRepeats.nDirBadRep = reset_Rep(AllRepeats.DirRep, AllRepeats.nDirRep, merge_repeats);
		}
		else{
			MyFree(AllRepeats.DirRep, AllSeeds->nDirSeeds*sizeof(Rep) );
			AllRepeats.DirRep=NULL;
		}
	
	}
	
	if(opt_inv){
		
		/*
			Get Inv sequence
		*/
		invseq2 = (char *)MyMalloc(size2*sizeof(char),"align_seeds: malloc error, bye");
		invseq(seq2,invseq2);
		
		for(s=0, r=0; s < AllSeeds->nInvSeeds; s++){
		
			if( s%10 == 0)
				fprintf(stderr,"inverted: %d seeds ->  %d repeats\r",s,r);
			
			if(AllSeeds->invSeeds[s].pos1 == -1)continue;              
			                          
			AllRepeats.InvRep[r] = aligni_2seq(AllSeeds->invSeeds[s].pos1,
			                                   AllSeeds->invSeeds[s].pos2,
							   AllSeeds->invSeeds[s].length,
							   seq1, invseq2,
							   size1, size2,
							   Xg, &Scoring, &Results);
			
			AllRepeats.InvRep[r].seed_meanR  = AllSeeds->invSeeds[s].rmean;

			reset_seeds(AllSeeds, 'i', s, AllRepeats.InvRep[r],  &Results);
			r++;	
			
			
		}
		
		AllRepeats.nInvRep = r;
		
		if(AllRepeats.nInvRep != 0){
			AllRepeats.InvRep = (Rep *)MyRealloc(AllRepeats.InvRep, AllRepeats.nInvRep*sizeof(Rep),
			                                    AllSeeds->nInvSeeds*sizeof(Rep),"align_seeds : realloc InvRep error");
			
			qsort_1then2_repeats(AllRepeats.InvRep, AllRepeats.nInvRep);
			
			AllRepeats.nInvBadRep = reset_Rep(AllRepeats.InvRep, AllRepeats.nInvRep, merge_repeats);	    
		}
		else{
			MyFree( AllRepeats.InvRep, AllSeeds->nInvSeeds*sizeof(Rep) );
			AllRepeats.InvRep=NULL;
		}
		
		MyFree(invseq2, size2*sizeof(char) );
		
	}
		
	free_Blast2Align(&Results);                       /* free the RESULTS struct */		

	return AllRepeats;
}

