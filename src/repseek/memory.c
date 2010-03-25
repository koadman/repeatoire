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
        file     : memory.c
        function : All about memory of the KMR, Seeds and Repeats
		             All MyMalloc() series is about follwoing how mauch memory has been used
		          
        created  : 19 Sep 2003
        modif    : Oct 2003, Feb 2004
        modif    : Dec 2004 <EC> ; Corrected Memory declaration
        author   : amikezor
*****/

#include <stdio.h>
#include <stdlib.h>
#include "repseek_types.h"
#include "memory.h"


MemUsage Memory;


/*
	Functions to count the memory usage all along
	dybamic allocation and free
*/
void PrintMem(char *Comment){

	extern MemUsage Memory;

	fprintf(stderr,"\n%s\nMemory Usage\n\t* Max is: %d bytes, %.2f Kb, %.2f Mb\n\t* Cur is: %d bytes, %.2f Kb, %.2f Mb\n",
	Comment,
	Memory.Max, (float)Memory.Max/1024, (float)Memory.Max/(1024*1024),
	Memory.Current, (float)Memory.Current/1024, (float)Memory.Current/(1024*1024));
}
void PrintMaxMem( void ){

	extern MemUsage Memory;

	if(Memory.Max < 1024)
		fprintf(stderr,"Max Memory Usage.. %d bytes\n", Memory.Max);
	else if(Memory.Max < 1024*1024)
		fprintf(stderr,"Max Memory Usage.. %.2f Kilobytes\n", (float)Memory.Max/1024);
	else if(Memory.Max < 1024*1024*1024)
		fprintf(stderr,"Max Memory Usage.. %.2f Megabytes\n", (float)Memory.Max/(1024*1024));
	else
		fprintf(stderr,"Max Memory Usage.. %.2f Gigabytes\n", (float)Memory.Max/(1024*1024*1024));
}
void Update_Mem(int32_t Value){

	extern MemUsage Memory;

	Memory.Current += Value;
	Memory.Max = (Memory.Current>Memory.Max)?Memory.Current:Memory.Max;
	
}
void Init_Mem(int32_t Value){

	extern MemUsage Memory;

	Memory.Current = Value;
	Memory.Max = Value;
}


/*
	Replace functions of dynamic allocation 
	to allow the tracking of memory usage
*/
void *MyMalloc( int32_t size , char *Error ){

	void *pointer;

	/* fprintf(stderr,"malloc: %s %d\n",Error, size); */

	pointer = malloc(size);
	if(!pointer)fprintf(stderr,"%s\n",Error),exit(3);

	Update_Mem(size);                                     /* fprintf(stderr,"Mem += %ld - malloc %s\n", size, Error); */
	
	return pointer;
}
void *MyCalloc( int32_t number, int32_t TypeSize , char *Error ){

	void *pointer;

	/* fprintf(stderr,"calloc: %s %d of size %d\n", Error, number, TypeSize); */
	
	pointer = calloc(number, TypeSize);
	if(!pointer)fprintf(stderr,"%s\n",Error),exit(3);

	Update_Mem(number*TypeSize );                       /* fprintf(stderr,"Mem += %ld - calloc %s\n", number*TypeSize, Error); */
	
	return pointer;
}
void MyFree( void *Pointer, int32_t size){

	if(!Pointer)return;
	
	/* fprintf(stderr,"free: %d\n", size); */

	free(Pointer);
	Pointer=NULL;

	Update_Mem(-size);                                   /* fprintf(stderr,"Mem -= %ld - free\n", size); */
}
void *MyRealloc( void *Pointer, int32_t newsize, int32_t oldsize, char *Error){


	/* fprintf(stderr,"realloc: %s %d (from %d)\n", Error, newsize, oldsize); */

	Pointer = realloc(Pointer,newsize);
	if(!Pointer)fprintf(stderr,"%s\n",Error),exit(3);

	Update_Mem( newsize-oldsize );                       /* fprintf(stderr,"Mem += %ld - realloc %s\n", (newsize-oldsize), Error); */
	return Pointer;
}


/*
	Deal with Stacks structure for KMR

void MallocStack(Stacks *s, int32_t number, int32_t *histo, int32_t AllValues){

	int32_t i;
		
	s->nStacks = number;
	s->nValues =  AllValues;
		
	s->ppStacks = (int32_t **)MyMalloc( number * sizeof(int32_t *), 
	                                  "MallocStack: ppStacks  malloc error, bye\n");
	s->lenStacks = (int32_t *)MyMalloc( number * sizeof(int32_t),
	                                "MallocStack: lenStacks malloc error, bye\n");
	s->ppStacks[0] = (int32_t *)MyMalloc( AllValues * sizeof(int32_t),
	                                  "MallocStack: ppStacks[0] malloc error, bye\n");
	s->lenStacks[0]=0;
	
	for(i=1;i < number; i++){
		s->lenStacks[i]=0;
		s->ppStacks[i] = s->ppStacks[i-1] + histo[i] ;	
	}
}

void FreeStack( Stacks *p){
	MyFree(p->ppStacks[0]  , p->nValues*sizeof(int32_t) );
	MyFree(p->ppStacks  ,    p->nStacks*sizeof(int32_t *) );
	MyFree(p->lenStacks ,    p->nStacks*sizeof(int32_t));
}
*/


/*
	Deal with the Seeds part


void free_Seeds(Seeds AllSeeds)    
{

	if( AllSeeds.nDirSeeds ){
		MyFree(AllSeeds.DirPos1, AllSeeds.nDirSeeds*sizeof(int32_t));
		MyFree(AllSeeds.DirPos2, AllSeeds.nDirSeeds*sizeof(int32_t));
		MyFree(AllSeeds.DirLen, AllSeeds.nDirSeeds*sizeof(int32_t));
		MyFree(AllSeeds.DirMeanR, AllSeeds.nDirSeeds*sizeof(float));
	}
	
	if(AllSeeds.nInvSeeds ){
		MyFree(AllSeeds.InvPos1, AllSeeds.nInvSeeds*sizeof(int32_t)); 
		MyFree(AllSeeds.InvPos2, AllSeeds.nInvSeeds*sizeof(int32_t));
		MyFree(AllSeeds.InvLen, AllSeeds.nInvSeeds*sizeof(int32_t));
		MyFree(AllSeeds.InvMeanR, AllSeeds.nInvSeeds*sizeof(float));
	}

}
*/
/*
	Malloc if it is the first time otherwise readjust

void AllocSeeds(Seeds *AllSeeds, int32_t size, int32_t old_size, int8_t opt_dir, int8_t opt_inv){

	if(opt_inv != 1 && opt_dir != 1)
		fprintf(stderr,"AllocSeeds: requiere at least one of opt_dir and opt_inv to be 1\n"),exit(4);
		
	if(opt_dir == 1){
		if( AllSeeds->DirPos1 == NULL){
			AllSeeds->DirPos1 = (int32_t *)MyCalloc(size , sizeof(int32_t),"AllocSeeds: Alloc for DirPos1 failed, bye");
			AllSeeds->DirPos2 = (int32_t *)MyCalloc(size , sizeof(int32_t),"AllocSeeds: Alloc for DirPos2 failed, bye");
		}
		else{
			AllSeeds->DirPos1 = (int32_t *)MyRealloc(AllSeeds->DirPos1, size * sizeof(int32_t), old_size* sizeof(int32_t),
			                                       "AllocSeeds: realloc for DirPos1 failed, bye");
			AllSeeds->DirPos2 = (int32_t *)MyRealloc(AllSeeds->DirPos2, size * sizeof(int32_t), old_size* sizeof(int32_t),
			                                       "AllocSeeds: realloc for DirPos2 failed, bye");
		}
	}

	if(opt_inv == 1){
		if( AllSeeds->InvPos1 == NULL){
			AllSeeds->InvPos1 = (int32_t *)MyCalloc(size , sizeof(int32_t), "AllocSeeds: Alloc for InvPos1 failed, bye");
			AllSeeds->InvPos2 = (int32_t *)MyCalloc(size , sizeof(int32_t), "AllocSeeds: Alloc for InvPos2 failed, bye");
		}
		else{
			AllSeeds->InvPos1 = (int32_t *)MyRealloc(AllSeeds->InvPos1, size * sizeof(int32_t), old_size* sizeof(int32_t),
			                                       "AllocSeeds: realloc for InvPos1 failed, bye");
			AllSeeds->InvPos2 = (int32_t *)MyRealloc(AllSeeds->InvPos2, size * sizeof(int32_t), old_size* sizeof(int32_t),
			                                       "AllocSeeds: realloc for InvPos2 failed, bye");
		}
	}
}
*/

/*
	Deal with the Repeats structure
*/
Repeats mem_Repeats(int32_t Ndir, int32_t Ninv){

	Repeats AllRepeats;             /* All Repeats structure */
	
	AllRepeats.nDirRep = Ndir;      /* set the number of repeats to the number of seeds */
	AllRepeats.nInvRep = Ninv;      
	AllRepeats.nDirBadRep = 0;      /* set the "bad" repet (included into another rep) as 0 */
	AllRepeats.nInvBadRep = 0;
	                                
	if(AllRepeats.nDirRep)
		AllRepeats.DirRep = (Rep *)MyMalloc( (AllRepeats.nDirRep)*sizeof(Rep), "init_Repeats: repdir malloc error" );
	else
		AllRepeats.DirRep = NULL;

	if(AllRepeats.nInvRep)
		AllRepeats.InvRep = (Rep *)MyMalloc( (AllRepeats.nInvRep)*sizeof(Rep), "init_Repeats: repinv malloc error" );
	else
		AllRepeats.InvRep = NULL;

	return AllRepeats;
}
void free_Repeats(Repeats AllRep)
{
	if(AllRep.nDirRep)
		MyFree(AllRep.DirRep, AllRep.nDirRep*sizeof(Rep));
	if(AllRep.nInvRep)
		MyFree(AllRep.InvRep, AllRep.nInvRep*sizeof(Rep));
}
