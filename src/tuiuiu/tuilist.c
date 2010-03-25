/* ***********************************************************************

   Copyright (C) 2010  Maria Federico, Pierre Peterlongo, Gustavo Sacomoto
   All Rights Reserved.

   This file is part of TUIUIU.

   TUIUIU is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   TUIUIU is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with TUIUIU.  If not, see <http://www.gnu.org/licenses/>.

************************************************************************/

#include <stdlib.h>
#include <stdio.h>

typedef struct TUILIST {
   int begin;
   int end;
   struct TUILIST *prox;
} tuilist;

tuilist *NewList(void)
{
   return NULL;
}

tuilist *Add(tuilist *l, int begin, int end)
{
   tuilist *new = (tuilist *)malloc(sizeof(tuilist));    
   
   new->begin = begin;
   new->end = end;
   new->prox = l;
   
   return new;
}

void FreeList(tuilist *l)
{
   tuilist *aux;
   
   while (l != NULL)
   {
      aux = l;
      l = l->prox;
      free(aux);
   }
}

tuilist *InvertList(tuilist *l)
{
   tuilist *aux, *inv = NULL; 
  
   while (l != NULL)
   {
      aux = inv;
      inv = l;
      l = l->prox;
      inv->prox = aux; 
   }
   
   

   return inv;      
}



inline int reval (int val, int middle){return (2*middle) - val -1;}

// manage the results while the reverse complement is also stored.
// For instance on a sequence of length 5, positions stored are the following:
// i   = 0  1  2  3  4  5  6  7  8  9 
// pos = 0  1  2  3 -4 -4 -3 -2 -1 -0 (with - standing for the reverse position)
// so from 0 to (size/2)-1   pos=i
// and from (size/2) to size -1: pos = 2*size - i - 1
// so results in a list like: 
// [0-2]->[6-7]->/ means that regions 0-2 is kept and reverse of region 6-7 (that is 2-3) is also kept.

tuilist * ManageResultsForReverse(tuilist *l, int size)
{
  tuilist * first_rev=NULL;	/* list from reversed */
  tuilist * init = l;		/* keep a track of the init list l */
  tuilist * rev;
  int middle = size/2;		/* size is necessary even */

  int this_end;

  tuilist * res;
  int number_in_direct = 0;

  ///////////////////// DEBUG     //////////////////////////

/*   fprintf(stderr," before all list l\n"); */
/*   while(l != NULL){ */
/*     fprintf(stderr,"[%d %d] ", l->begin, l->end); */
/*     l = l -> prox; */
/*   } */
/*   printf("\n"); */
/*   l = init; */
  ///////////////////// DEBUG     //////////////////////////
  
  // modify the list if begin or end is bigger than middle
  while(l != NULL){
    // 3 cases:
    //   - begin and end are smaller or equal to  middle. Do nothing (not represented in the code)
    //   - begin is smaller or equal to middle and end if bigger:
    //        . create a cell begin -> middle
    //        . link it to a new cell [0, rev(end)]
    //        . link this new cell to the historical next of the cell
    //   - begin and end are bigger than middle, replace by [rev(end), rev(begin)]
    if(l->begin<=middle && l->end<=middle) number_in_direct++;
    
    if(l->begin<=middle && l->end>middle){
      number_in_direct++;
      this_end = l->end;
      l->end=middle-1;
      l->prox = Add(l->prox, reval(this_end,middle), middle-1);
      first_rev = l->prox;
    }
    if(l->begin > middle && l->end>middle){
      this_end = l->end;
      l->end = reval(l->begin,middle);
      l->begin = reval (this_end, middle);
      if(first_rev==NULL) first_rev = l;
    }
    
    l=l->prox;
  }
  l = init;
  rev = InvertList(first_rev);
  first_rev = rev;

  // create a new list to be returned;
  // in list l, number_in_direct have to be considered
  // in list rev, all cells are considered
  // this lists have to be merged in order to sort their content.
  // for instance:
  // l : (1 10) -> (35 50) -> (70 100) -> ... (with number_in_direct = 3)
  // rev: (15 17) -> (40 61) -> (65 110) -> /
  // res : (1 10) -> (15 17) -> (35 61) -> (65 110) -> /
  res = NewList();
  int current_begin;
  int current_end;
  int inside_l=0;
  

  ///////////////////// DEBUG     //////////////////////////
/*   fprintf(stderr," list l, (size limited to %d)\n", number_in_direct); */
/*   while(l != NULL){ */
/*     fprintf(stderr,"[%d %d] ", l->begin, l->end); */
/*     l = l -> prox; */
/*   } */
/*   printf("\n"); */
/*   l = init; */

/*   fprintf(stderr," list rev \n"); */
/*   first_rev = rev; */
/*   while(rev != NULL){ */
/*     fprintf(stderr,"[%d %d] ", rev->begin, rev->end); */
/*     rev = rev -> prox; */
/*   } */
/*   printf("\n"); */
/*   rev = first_rev; */
  ///////////////////// END DEBUG //////////////////////////

  // walk the two lists, mergin them
  while(inside_l<number_in_direct && l != NULL && rev != NULL){
    if(l->begin < rev->begin)  current_begin = l->begin;
    else  current_begin = rev->begin;
    
    // goto the good end cell.
    while(inside_l<number_in_direct && l != NULL && rev != NULL){

      // case l totally in rev
      // l        ---------
      // rev  --------------------
      if(l->begin >= rev->begin && l->end <= rev->end){l = l->prox; inside_l++; continue;}


      // case rev totally in l
      // l   --------------------
      // rev       ---------
      if(rev->begin > l->begin && rev->end < l->end){rev = rev->prox; continue;}

      // case starting by l but rev overlaps (or touches) l
      // l    ------------
      // rev         ------------
      if(l->begin <= rev->begin && l->end >= rev->begin){l = l->prox; inside_l++; continue;}

      // case starting by rev but l overlaps (or touches) rev
      // l           ------------
      // rev  ------------
      if(rev->begin <= l->begin && rev->end >= l->begin){rev = rev->prox; continue;}

      // case starting by l and non overlapping with rev
      // l   ---------
      // rev             ----------
      if(l->begin < rev->begin && l->end < rev->begin){
	current_end=l->end;
	l = l->prox;
	inside_l++;
	break;
      }
      
      // case starting by rev and non overlapping with l
      // l               ----------
      // rev --------
      if(rev->begin < l->begin && rev->end < l->begin){
	current_end=rev->end;
	rev = rev->prox;
	break;
      }
      fprintf(stderr,"unrecheable code\n");
      // in theory, code unreachable here!
    } // end searching the goo end cell
    // end of list cases:
    
    // case end of rev
    // l    ------------
    // rev   
    if(rev==NULL){
      current_end = l->end;
      l = l->prox;
      inside_l++;
    }

    // case end of l
    // l
    // rev  -----------
    else if(inside_l == number_in_direct){
      current_end=rev->end;
      rev = rev->prox;
    }    
    res = Add(res, current_begin, current_end);
  }

  // finishing each list
  
  // case end of rev
  // l    ------------
  // rev   
  if(rev==NULL)
    while(inside_l<number_in_direct && l != NULL){
      res = Add(res, l->begin, l->end);
      l=l->prox;
      inside_l++;
    }
    
  // case end of l
  // l
  // rev  -----------
  if(inside_l == number_in_direct)
    while(rev !=NULL){
      res = Add(res, rev->begin, rev->end);
      rev = rev->prox;
    }
  
  FreeList(init);
  FreeList(rev);
  res = InvertList(res);
  return res;
}

/* int main(void)
{
   tuilist *l = NewList(), *w; 

   l = Add(l, 0, 1);
   l = Add(l, 1, 2);
   l = Add(l, 2, 3);
   l = Add(l, 3, 4);
   
   w = l;
   while(w != NULL)
   {
      printf("[%d, %d]", w->begin, w->end);
      w = w->prox;
   }
   printf("\n\n");
   l = InvertList(l); 
   w = l;
   while(w != NULL)
   {
      printf("[%d, %d] ", w->begin, w->end);
      w = w->prox;
   }    
   printf("\n\n");
   FreeList(l);
   return 0;
} */
