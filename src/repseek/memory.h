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
 * @file   memory.h
 * @author Achaz G
 * @date   April 2004
 * 
 * @brief  header for memory alloc/dealloc
 * modif    : Dec 2004 <EC> ; Corrected Memory declaration
 * 
 * 
 */



#ifndef _MEMORY_H_
#define _MEMORY_H_

#include "repseek_types.h"

typedef struct  {           /********** Memory Usage structure **************/

	int32_t Max;
	int32_t Current;

} MemUsage;



#include <stdio.h>
#include <stdlib.h>

#if defined(DMALLOC)
#include <dmalloc.h>
#endif

/********** **********

      Global Variable(s)
		
 ********** **********/

extern MemUsage Memory;  /* Instance of the global variable for memory tracking */

/*
	Follow memory usage
*/
void PrintMem(char *Comment);
void PrintMaxMem( void );
void Update_Mem(int32_t Value);
void Init_Mem(int32_t Value);
/*
	All Alloc/Free to follow of memory usage
*/
void *MyMalloc( int32_t size , char *Error );
void *MyCalloc( int32_t number, int32_t TypeSize , char *Error );
void MyFree( void *Pointer, int32_t size);
void *MyRealloc( void *Pointer, int32_t newsize, int32_t oldsize, char *Error);
/*
	For Stacks

void MallocStack(Stacks *s, int32_t number, int32_t *histo, int32_t AllValues);
void FreeStack( Stacks *p);

	For Seeds

void free_Seeds(Seeds AllSeeds);
void AllocSeeds(Seeds *AllSeeds, int32_t size, int32_t old_size, int8_t opt_dir, int8_t opt_inv);
*/

/*
	For Repeats
*/
Repeats mem_Repeats(int32_t Ndir, int32_t Ninv);
void free_Repeats(Repeats AllRep);




/*
	Not used anymore, but just in case
*/


#define KMRK_MALLOC(var,type,size,message) {       \
            var = (type*) malloc(size);            \
            if (!var)                              \
              {                                    \
		fprintf(stderr,"%s\n",message);    \
		exit(4);                           \
	      }                                    \
            }

#define KMRK_CALLOC(var,type,length,message) {             \
            var = (type*) calloc(length,sizeof(type));     \
            if (!var)                                      \
              {                                            \
		fprintf(stderr,"%s\n",message);            \
		exit(4);                                   \
	      }                                            \
            }
                
#define KMRK_REALLOC(var,type,size,message) {              \
            var = (type*) realloc(var,size);               \
            if (!var)                                      \
              {                                            \
		fprintf(stderr,"%s\n",message);            \
		exit(4);                                   \
	      }                                            \
            }


#endif  /* _MEMORY_H_*/
