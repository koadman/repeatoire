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
 * @file   sort.c
 * @author Guillaume Achaz <gachaz@oeb.harvard.edu>
 * @date   oct 03 2002
 * @modif  dec 13 2002 - add the sort_repeats() and qsort_1then2_repeats()
           Oct 2003
			  April 2003  - Int32->int32_t
			  July 2004, change Qsort into double_Qsort
 * 
 * @brief   sort two arrays by values of the first one (int or Rep) - Implemented from Numrec 92 pp.333/336
 * 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include "repseek_types.h"
#include "memory.h"

/***
	sort two array by the order of the first one
	Qsort is implemented following the algorithm explained in Numrec C, pp 332-335
	Under a certain size, a part of the array is sorted by the naive method (straight insertion)
***/

#define NSTACK 1000                           /* the number of subarray that could be memorized */
#define SWITCH_METHOD 7                       /* under that size switch to a Insertion method */
#define SWAP(a,b) {tmp=(a); (a)=(b); (b)=tmp;}

static void qsort2(int32_t *array1, int32_t *array2, int32_t n){

	
	int32_t *Stack;       /* where remaining values of subarrays are temporarily stocked */
	int32_t nStack=0;     /* How many of these values are in the Stack. Is never odd */
	
	int32_t end=n-1,      /* the last value of the array to sort */
	     beg=0,        /* the first value ---  */
		  postbeg;      /* the second value ---  */

	int32_t val_postbeg;  /* the value of the postbeg position - the one which is used for partion-exchange */

	int32_t demi_len;     /* half of the distance between beg and end */
	
	int32_t i,            /* counter from postbeg to right */
	     j;            /* counter from end to left */

	int32_t val1_i,       /* for insertion stock temporarily a value */ 
		  val2_i;
	
	int32_t tmp;          /* used for the SWAP macro */
	
	
	Stack = (int32_t *)MyMalloc(NSTACK*sizeof(int32_t) , "qsort2: not enough memory for Stack") ; 
		
	while(1){ 

		if( end-beg+1 > SWITCH_METHOD){
	
			demi_len = (end-beg) >> 1 ; 
			postbeg = beg+1;
			
			SWAP( array1[beg+demi_len], array1[postbeg] );
			SWAP( array2[beg+demi_len], array2[postbeg] );
			
			if(array1[beg] > array1[postbeg]){            /* rearrange to have  beg <= postbeg <= end */ 
				SWAP( array1[beg], array1[postbeg] );
				SWAP( array2[beg], array2[postbeg] );
			}

			if(array1[beg] > array1[end]){
				SWAP( array1[beg], array1[end] );
				SWAP( array2[beg], array2[end] );
			}

			if(array1[postbeg] > array1[end]){
				SWAP( array1[postbeg], array1[end] );
				SWAP( array2[postbeg], array2[end] );
			}
			
			
			i = postbeg;
			j = end;
						
			val_postbeg =  array1[postbeg];
			
			while(1)                                   /* enter the partition exchange process */
				{
					do i++; while( array1[i] < val_postbeg );
					do j--; while( array1[j] > val_postbeg );
					
					if(j<i) break;
					
					SWAP( array1[i], array1[j] );
					SWAP( array2[i], array2[j] );
				}
						
			SWAP( array1[postbeg] , array1[j] );   /* place the postbeg value into j */
			SWAP( array2[postbeg] , array2[j] );
			
			if(nStack+2 > NSTACK)
				fprintf(stderr, "qsort2: not enough Stack... sorry bye\n"),exit(1);
			
			if(end-i+1 >= j-beg){
				Stack[nStack++] = i;                /* stock in Stack the largest and go with the smallest */
				Stack[nStack++] = end;
				
				end = j-1;
			}
			else{
				Stack[nStack++] = beg;
				Stack[nStack++] = j-1;
				
				beg = i;
			}
			
		}
		else{                                      /* Under a certain size, switch to the straight insertion method */
	
			
			for(i=beg+1; i <= end ; i++){
			
				val1_i = array1[i];
				val2_i = array2[i];
				
				for(j=i-1;j>=beg;j--){
				
					if(array1[j] < val1_i)break;
					
					array1[j+1] = array1[j];
					array2[j+1] = array2[j];
				}
				array1[j+1]=val1_i;
				array2[j+1]=val2_i;
			}

			if(nStack==0)break;          /* this si the end  - exit the process */
			
			end = Stack[--nStack];       /* go for the next round with the stacked parameters */
			beg = Stack[--nStack];
		}

	 }

	MyFree(Stack, NSTACK*sizeof(int32_t));

}


/*
	A simple qsort as described in numrec
*/
void double_Qsort(double *array1, int32_t n){

	
	double *Stack;        /* where remaining values of subarrays are temporarily stocked */
	int32_t nStack=0;     /* How many of these values are in the Stack. Is never odd */
	
	int32_t end=n-1,      /* the last value of the array to sort */
	     beg=0,           /* the first value ---  */
		  postbeg;         /* the second value ---  */

	double val_postbeg;  /* the value of the postbeg position - the one which is used for partion-exchange */

	int32_t demi_len;     /* half of the distance between beg and end */
	
	int32_t i,            /* counter from postbeg to right */
	     j;              /* counter from end to left */

	double val1_i;       /* for insertion stock temporarily a value */ 
	
	double tmp;          /* used for the SWAP macro */
	
	
	Stack= (double *)MyMalloc(NSTACK*sizeof(double) , "Qsort: not enough memory for Stack" ); 
		
	while(1){ 

		if( end-beg+1 > SWITCH_METHOD){
	
			demi_len = (end-beg) >> 1 ; 
			postbeg = beg+1;
			
			SWAP( array1[beg+demi_len], array1[postbeg] );
			
			if(array1[beg] > array1[postbeg]){            /* rearrange to have  beg <= postbeg <= end */ 
				SWAP( array1[beg], array1[postbeg] );
			}

			if(array1[beg] > array1[end]){
				SWAP( array1[beg], array1[end] );
			}

			if(array1[postbeg] > array1[end]){
				SWAP( array1[postbeg], array1[end] );
			}
			
			
			i = postbeg;
			j = end;
						
			val_postbeg =  array1[postbeg];
			
			while(1)                                   /* enter the partition exchange process */
				{
					do i++; while( array1[i] < val_postbeg );
					do j--; while( array1[j] > val_postbeg );
					
					if(j<i) break;
					
					SWAP( array1[i], array1[j] );
				}
						
			SWAP( array1[postbeg] , array1[j] );   /* place the postbeg value into j */
			
			if(nStack+2 > NSTACK)
				fprintf(stderr, "Qsort: not enough Stack... sorry bye\n"),exit(1);
			
			if(end-i+1 >= j-beg){
				Stack[nStack++] = i;                /* stock in Stack the largest and go with the smallest */
				Stack[nStack++] = end;
				
				end = j-1;
			}
			else{
				Stack[nStack++] = beg;
				Stack[nStack++] = j-1;
				
				beg = i;
			}
			
		}
		else{                                      /* Under a certain size, switch to the straight insertion method */
	
			
			for(i=beg+1; i <= end ; i++){
			
				val1_i = array1[i];
				
				for(j=i-1;j>=beg;j--){
				
					if(array1[j] < val1_i)break;
					
					array1[j+1] = array1[j];
				}
				array1[j+1]=val1_i;
			}

			if(nStack==0)break;          /* this si the end  - exit the process */
			
			end = Stack[--nStack];       /* go for the next round with the stacked parameters */
			beg = Stack[--nStack];
		}

	 }

	MyFree(Stack, NSTACK*sizeof(double));

}

#undef NSTACK
#undef SWITCH_METHOD
#undef SWAP


/*
	Sort all repeat by an array
	Just modify from qsort2
*/

#define SWAP(a,b)  { tmprep=(a); (a)=(b); (b)=tmprep; }
#define NSTACK 1000                           /* the number of subarray that could be memorized */
#define SWITCH_METHOD 7                       /* under that size switch to a Insertion method */

static int32_t GetPos1( Rep repet ){
	return repet.pos1;
}

static int32_t GetPos2( Rep repet ){
	return repet.pos2;
}

static void qsort_repeats(Rep *AllReps, int32_t n, int32_t pos ){
	
	int32_t *Stack;       /* where remaining values of subarrays are temporarily stocked */
	int32_t nStack=0;     /* How many of these values are in the Stack. Is never odd */
	
	int32_t end=n-1,      /* the last value of the array to sort */
	     beg=0,        /* the first value ---  */
		  postbeg;      /* the second value ---  */

	int32_t val_postbeg;  /* the value of the postbeg position - the one which is used for partion-exchange */

	int32_t demi_len;     /* half of the distance between beg and end */
	
	int32_t i,            /* counter from postbeg to right */
	     j;            /* counter from end to left */

	Rep val_i;         /* for insertion stock temporarily a value */ 
	Rep tmprep;        /* used for the SWAP macro */
		
	int32_t (* PtrFuncPos)( Rep repet );
		
	if( pos == 1)
		PtrFuncPos =  GetPos1;
	else if( pos == 2)
		PtrFuncPos =  GetPos2;
	else
		fprintf(stderr,"qsort_repet: pos should be 1 or 2\n"),exit(4);
	
	Stack=(int32_t *)MyMalloc(NSTACK*sizeof(int32_t) , "qsort_repet: not enough memory for Stack" ); 

	while(1){

		if( end-beg+1 > SWITCH_METHOD){
	
			demi_len = (end-beg) >> 1 ; 
			postbeg = beg+1;
			
			SWAP( AllReps[beg+demi_len], AllReps[postbeg] );
			
			if( PtrFuncPos(AllReps[beg]) > PtrFuncPos(AllReps[postbeg]) ){            /* rearrange to have  beg <= postbeg <= end */ 
				SWAP( AllReps[beg], AllReps[postbeg] );
			}
			
			if(  PtrFuncPos(AllReps[beg]) >  PtrFuncPos(AllReps[end]) ){
				SWAP( AllReps[beg], AllReps[end] );
			}
			
			if(  PtrFuncPos(AllReps[postbeg]) >  PtrFuncPos(AllReps[end]) ){
				SWAP( AllReps[postbeg], AllReps[end] );
			}
			
			i = postbeg;
			j = end;
						
			val_postbeg =   PtrFuncPos(AllReps[postbeg]);
			
			while(1)                                   /* enter the partition exchange process */
				{
					do i++; while(   PtrFuncPos(AllReps[i]) < val_postbeg );
					do j--; while(   PtrFuncPos(AllReps[j]) > val_postbeg );
					
					if(j<i) break;
					
					SWAP( AllReps[i], AllReps[j] );
				}
						
			SWAP( AllReps[postbeg] , AllReps[j] );   /* place the postbeg value into j */
			
			if(nStack+2 > NSTACK)
				fprintf(stderr, "qsort_repet: not enough Stack... sorry bye\n"),exit(1);
			
			if(end-i+1 >= j-beg){
				Stack[nStack++] = i;                /* stock in Stack the largest and go with the smallest */
				Stack[nStack++] = end;
				
				end = j-1;
			}
			else{
				Stack[nStack++] = beg;
				Stack[nStack++] = j-1;
				
				beg = i;
			}
			
		}
		else{                                      /* Under a certain size, switch to the straight insertion method */
	
			
			for(i=beg+1; i <= end ; i++){
			
				val_i = AllReps[i];
				
				for(j=i-1;j>=beg;j--){
				
					if(  PtrFuncPos(AllReps[j]) <= PtrFuncPos(val_i) )break;
					
					AllReps[j+1] = AllReps[j];
				}
				AllReps[j+1]=val_i;
			}

			if(nStack==0)break;          /* this si the end  - exit the process */
			
			end = Stack[--nStack];       /* go for the next round with the stacked parameters */
			beg = Stack[--nStack];
		}

	 }

	MyFree(Stack, NSTACK*sizeof(int32_t));

}

#undef NSTACK
#undef SWITCH_METHOD
#undef SWAP


/***
	This function sort the first array, and for equal value, use the second one
	calling the qsort2() function
***/
void qsort2_1then2(int32_t *array1, int32_t *array2, int32_t n){

	int32_t i;
	int32_t beg=0, end=n;

	char token=0;
	
	qsort2(array1, array2, n);

	for(i=1, token=0; i < n-1 ; i++){	
	
		if(array1[i] == array1[i-1]){
			if(!token) beg = i-1;
			token = 1;
		}
		
		if(array1[i] != array1[i-1])
			if(token){
				token=0;
				end=i;
				qsort2(array2+beg, array1+beg, end-beg);
			}
	}
	
	if(token){
		end=i;
		qsort2(array2+beg, array1+beg, end-beg);
	}
	
}


/***
	This function sort by pos1 and for equal value, use pos2
	calling the qsort2_repeats() function
***/

void qsort_1then2_repeats(Rep *AllReps, int32_t n){

	int32_t i;
	int32_t beg=0, end=n;

	char token=0;
	
	qsort_repeats(AllReps, n, 1 );

	for(i=1, token=0; i < n-1 ; i++){	
	
		if( AllReps[i].pos1 == AllReps[i-1].pos1 ){
			if(!token) beg = i-1;
			token = 1;
		}
		
		if(AllReps[i].pos1 != AllReps[i-1].pos1)
			if(token){
				token=0;
				end=i;
				qsort_repeats(AllReps+beg, end-beg, 2 );
			}
	}
	
	if(token){
		end=i;
		qsort_repeats(AllReps+beg, end-beg, 2 );
	}
	
}


