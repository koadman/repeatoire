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

//#include "index_str.h"
//#include "tuilist.h"
//#include "itree.h"
//#include "util.h"
#include "mono_and_multi_commons.h"
#include <stdlib.h>
#include <stdio.h>
extern int tuiuiu_z;
//extern int bin;
//extern int bin_size;
//extern float w;
#define BEGIN(i, bin, z, N) (i + N - ((bin + 1) << z))
//#define END(i, bin, z, w, N) (i + w + N - (bin << z))
#define END(begin, w, bin_size) (begin + w + bin_size - 1)
#define MAX(a, b) (a > b ? a : b)
#define PLCS_type short

char method[]="excelent mono multiPass";

int erro;
int max;
//int **c;

int *a = NULL;				/* label of each word on A */
int *b = NULL;				/* label of each word on B */
//int minSize;			/* minimal size between a and b */

int * listLSI=NULL;			/* to avoid to re-alloc many times */
int * occurrencesList;

int limit = 14;

/*
//        A
//   ------------
//   |\  \
//   | \  \
//   |  \  \
// B |   \  \
//   |    \  \
//   |     \  \
//   |      \  \
//   |       ---
*/
///////////////////////////////////////////////////////////////////////////
/////////// PART ABOUT LCS COMPUTATIONS ///////////////////////////////////
/////////// SHOULD BE IN AN OTHER FILE  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////


			     
// return the index of the first occurences in the ordered array bigger (or equal) than 
// the value
inline int findFirstAfter (const int * tab, int start, int stop, const int value)
{
   int pos, pivot;
   
   if (tab[stop] < value) return -1; 
   
   while ((stop - start) >= limit)
   {
      pivot = (start + stop) >> 1;
      if(tab[pivot] < value)
	 start = pivot;
      else
	 stop = pivot;
   }
  
   for (pos = start; pos <= stop; pos++)
      if(tab[pos] >= value) return pos;
   return -1;
}


/* return the index of the first occurences in the inverted ordered
   array lower (or equal) than the value */
inline int findFirstLowerAfter (const int * tab, int start, int stop, const int value)
{
   int pos, pivot;
   
   if (value < tab[stop]) return -1;    

   while ((stop - start) >= limit)
   {
      pivot = (start + stop) >> 1;
      if (tab[pivot] > value)
         start = pivot;
      else
         stop = pivot;
   }
   
   for (pos = start; pos <= stop; pos++)
      if(tab[pos] <= value) return pos;
   return -1;
}



// compute the longuest increasing subsequence of the list
int simplelisparallelogram (){

  
  int nbInList = 0;		// number of lines in the list
  int posInList;
  int indexInOccurrenceList=0;
  
  // each occurrences 
  while(occurrencesList[indexInOccurrenceList]!=-1){ 
    //    fprintf(stderr,"cur = %d\n", cur->info); /* DEB */
    if(nbInList==0) posInList=-1;	/* empty list */

    else posInList = findFirstAfter (listLSI, 0, nbInList - 1, occurrencesList[indexInOccurrenceList]); /* first position bigger than biggest number in list (if possible) */
	
    if(posInList == -1) 
      listLSI[nbInList++] = occurrencesList[indexInOccurrenceList];
    //	if(nbInList>=p) return 1;
    
    else listLSI[posInList] = occurrencesList[indexInOccurrenceList];
    indexInOccurrenceList++; // go to the next position on b value
  }
	
  return nbInList; 
}




///////////////////////////////////////////////////////////////////////////
/////////// END PART ABOUT LCS COMPUTATIONS ///////////////////////////////
///////////////////////////////////////////////////////////////////////////


/**
 * Returns the number of number of elements in the tree that are good parallelograms
 */
int Nelems(int *begin, int *end, itree *l, int i, int w, int bin_size, int N, empty_block goodWindows[], int curr_pass_num, int curr_parall_min, int curr_parall_max)
{
  //  printf("Nelems w: %d\n",w);
  //  printf("Nelems tuiuiu_z: %d\n",tuiuiu_z);
  int result, n, intervalSize, beg_currBin, real_end_currBin, parall_ind, j, neb, begTemp;

   if (l == NULL)
      return 0;
   else
   {
     result = Nelems(begin, end, l->right, i, w, bin_size, N, goodWindows, curr_pass_num, curr_parall_min, curr_parall_max);

     //     printf("Nelems2 w: %d\n",w);
     //     printf("Nelems2 tuiuiu_z: %d\n",tuiuiu_z);
     //real_end_currBin = (*end);

     beg_currBin = BEGIN(i, l->bin, tuiuiu_z, N);
     //     printf("BEGIN tuiuiu_z: %d\n",tuiuiu_z);
     //real_end_currBin = END(i, l->bin, tuiuiu_z, w, N);
     real_end_currBin = END(beg_currBin, w, bin_size);
     //     printf("END tuiuiu_z: %d\n",tuiuiu_z);
     if(real_end_currBin > N)
       real_end_currBin = N;

     neb = 0;

     //if (BEGIN(i, l->bin, tuiuiu_z, N) >= (*end))
     //     printf("beg_currBin\n");
     if(beg_currBin >= (*end))
       //bin is non-overlapping
       {
	 if(real_end_currBin - beg_currBin >= (w-erro))
	   {
	     /* EMPTY BLOCK */
	     // a block has to be count only if it is not an
	     // empty block
	     //
	     // in the first pass: the empty-block info is available
	     // in goodWindows only for blocks of index less than the
	     // block containing the current windows (currBin)
	     //
	     // in successive passes: the empty-block info is significative
	     // also for blocks of index greater than currBin
	     //
	     // in either cases, if the \\ of the bin is the same
	     // of parall of the window we cannot use the fact that it is empty 
	     // because it depends also from the current window (in particular
	     // in pass successive to the first one, it was not empty at previous pass 
	     // because it contains at least the current window (kept in previous pass)
	     // [but it could be reset to empty block from some previous window in the
	     // same \\ in the current pass]
	     
	     parall_ind = beg_currBin >> tuiuiu_z;
	     //printf("parall_ind >>\n");
	     if(((parall_ind == curr_parall_min) || (parall_ind == curr_parall_max)) 
		||
		((curr_pass_num==1) && 
		 (((parall_ind < curr_parall_min) && ((goodWindows[parall_ind]).notEmpty!=0))
		  || (parall_ind > curr_parall_max)))
		|| 
		((curr_pass_num > 1) && ((goodWindows[parall_ind]).notEmpty!=0))) {
	     
	       result++;
	       
	       //*begin = BEGIN(i, l->bin, tuiuiu_z, N);
	       //(*begin) = beg_currBin;
	       //MOD
	       (*begin) = beg_currBin + (w-erro);
	       (*end) = real_end_currBin;
	     }
	   }
       }
     else // bin could be overlapping 
       {
	 if(beg_currBin > (*begin))
	   (*begin) = beg_currBin;
	 //	 printf("intervalSize\n");
	 //intervalSize = (END(i, l->bin, tuiuiu_z, w, N)) - (*begin);
	 intervalSize = real_end_currBin - (*begin);
	 
	 if (intervalSize >= (w-erro)) // bin is non-overlapping 
	   { 
	     n = intervalSize / (w-erro);         // the weight of each interval
	     
	     begTemp = (*begin);
	     
	     for(j=1; j<=n; j++){
	       
	       parall_ind = begTemp >> tuiuiu_z;
	       /* EMPTY BLOCK */
	       //               printf("parall_ind\n");
	       if(((parall_ind == curr_parall_min) || (parall_ind == curr_parall_max)) 
		  ||
		  ((curr_pass_num==1) && 
		   (((parall_ind < curr_parall_min) && ((goodWindows[parall_ind]).notEmpty!=0))
		    || (parall_ind > curr_parall_max)))
		  || 
		  ((curr_pass_num > 1) && ((goodWindows[parall_ind]).notEmpty!=0))) {
		 
		 result++;
		 
		 neb=j;
	       }
	       
	       begTemp = begTemp + (w-erro);
	     }
	     
	     //result = result + n;
	     
	     //(*begin) = (*begin) + n * (w-erro);  /* the rest */
	     (*begin) = (*begin) + neb * (w-erro);  /* the rest */ 
	     //neb is the LAST not empty block
	     (*end) = real_end_currBin;
	   }
	 //else{
	 // bin is overlapping
	 //}
       }
     
     //(*end) = real_end_currBin;
     //     printf ("return result\n");
     return result + Nelems(begin, end, l->left, i, w, bin_size, N, goodWindows, curr_pass_num, curr_parall_min, curr_parall_max);
   }
}

/**
 * Computes the Longest common sequence between seq[A..] and seq[B..]. 
 * Sequences are translated to the k-mer labels.
 */  
int  LCS(const int begin_A,
	       const int begin_B,
	       const int w,
	       const int tam,
	       const int k,
	       const char seq[],
	       const index_str * k_factor_ind)
{
  int i, A;
  int res;			/* Store the result */
  int firstIndex, lastIndex;	/* index positions in k_factor_ind->list */
  int smalestPositionInPara;	/* first possible occurrence of a word in a parrallelogram */
  int biggestPositionInPara;    /* last possible occurrence of a word in a parrallelogram */
  int indexNumber;		/* used for walking the k_factor_ind->list array */
  int indexInOccurrenceList=0;	/* position in the occurrencesList array */
  int nextNpos = 0;
  extern char ACTGnumber[];
  A = 0; 
  for (i = 0; i < k; i++)
  {
    /* N handling */ 
    if (seq[begin_A + i] == 'N')
       nextNpos = i + k + 1;
         
    /* list of positions for one word in a */
    A = A + (ACTGnumber[(int)seq[begin_A + i]] * pot4(i));
  }
   
  for (; i <= w; i++) /* all a positions */ //BUG <= instead of <
  {
     if (i >= nextNpos) /* N handling */
     {  
       firstIndex = k_factor_ind->ind[A];	/* first occurrence of A in n */
       lastIndex = k_factor_ind->ind[A+1]; /* last (+1) occurrence of A in n */
       smalestPositionInPara = begin_B + i - k; /* smalest position in the parrallelogram */
       biggestPositionInPara = begin_B + i + tam -1 -k; /* biggest reachable position in the parrallelogram  */
       
       // Finds in k_factor_ind->list the first index who
       //  - is in [firstIndex, lastIndex]
       //  - such that k_factor_ind->tuilist[firstIndex] <= biggestPositionInPara
       indexNumber = findFirstLowerAfter (
					  k_factor_ind->tuilist,
				           firstIndex,
				           lastIndex - 1,
				           biggestPositionInPara);
    
        if(indexNumber != -1) /* at least one occurrence found */
        {	
           while(indexNumber<lastIndex && /* still for the good A */
		 k_factor_ind->tuilist[indexNumber]>=smalestPositionInPara) /* still in the parrallelogram */
	     occurrencesList[indexInOccurrenceList++] = k_factor_ind->tuilist[indexNumber++]; /* add in the occurrences list */
	//addAtEnd(occurrencesList, k_factor_ind->tuilist[indexNumber++]); /* add in the occurrences list */
        } 
     } 
     /* shift the current A word */
     A = (A >> 2) + (ACTGnumber[(int)seq[begin_A + i ]] * pot4(k - 1));
     
     /* N hadling */
     if (seq[begin_A + i] == 'N')
        nextNpos = i + k + 1;     
  } // end all a positions

  // last position in occurrencesList
  occurrencesList[indexInOccurrenceList]=-1;
  
  // compute LIS on the occurrencesList
  
  res = simplelisparallelogram(occurrencesList);
  //printf ("lcs %5d %5d %2d\n", begin_A, begin_B, res); 
  return res;
}


/**
 * Construct a tree of good parallelograms.
 * WARNING ! here the term good parallelogram seems unappropriate as we also compute LCS. We should better use perfect parallelogram.
 * Would merit more detailled comments from brasilians. Gustavo and Alair...
 */
void goodPar(itree *l, int nextBegin[], PLCS_type pLCS[], int begin, int end, int p, int w, int N, 
             int k, int tam, char seq[], itree **l2, const index_str* k_factor_ind)
{
  PLCS_type plcs; 
  
   if (l != NULL)
   { 
     
     if (begin < nextBegin[l->bin])
       {
	 if (pLCS[l->bin] <= 0){
	   AddTree(l->bin, l2);
	 }   
       }
     else
       {   
	 pLCS[l->bin] = plcs = p - LCS(begin, BEGIN(begin, l->bin, tuiuiu_z, N), w, tam, k, seq, k_factor_ind);
	 if (plcs<=0) {
	   AddTree(l->bin, l2);
	   nextBegin[l->bin] = begin + 1 - plcs;
	 } 
	 else 
	   nextBegin[l->bin] = begin + plcs;
       }
   
     goodPar(l->left, nextBegin, pLCS, begin, end, p, w, N, k, tam, seq, l2, k_factor_ind);
     goodPar(l->right, nextBegin, pLCS, begin , end, p, w, N, k, tam, seq, l2, k_factor_ind); 
   }
}


/** 
 * Hear of the algorithm, computes the link between windows and
 * returns a list of positions contaiing the information about the
 * filtered sequences
*/
tuilist *Filter(int N, int k, int p, int e, int r, int bin_size, 
             int w, index_str *k_factor_ind, char seq[],
	     empty_block goodWindows[], int curr_pass_num)
{
  int  *bins, i, j, last = 0, next = 0, begin = 0, end = 0, nextNpos = 0, lastNpos = 0, curr_parall_min=0, curr_parall_max=0, curr_parall;
  PLCS_type *pLCS;
  int *nextBegin;
  extern char ACTGnumber[];
   //   int nb_tests_good_par=0;	/* DEB */

   // ADDED BY PIERRE
   //occurrencesList = (int *) malloc (sizeof(int)*(((bin_size + e + 1)*w)+1));  
   int maxOccListSize = ((bin_size + e + 1)*w)+1;
   maxOccListSize *= maxOccListSize;
   //occurrencesList = (int *) malloc (sizeof(int)*maxOccListSize);
   occurrencesList = (int *) calloc (maxOccListSize, sizeof(int));

   //listLSI = (int*) malloc (sizeof(int)*(w-k+1+1));
   listLSI = (int*) calloc ((w-k+1+1), sizeof(int));
   // END ADDED BY PIERRE
   
   // Initialisation
   tuilist *result = NewList();
   itree *l = NewTree(), *l2;   
   
   
   erro = e;
   max = ((2 * N) / bin_size) + 2;
   
   bins = (int *)calloc(max, sizeof(int));
   pLCS = (PLCS_type *)calloc(max, sizeof(PLCS_type));  /* pLCS = p - lastLCS */
   nextBegin = (int *)calloc(max, sizeof(int)); 
      
   // Compute k factors that occur in the first window (window at position 0)
   // compute first k-mfactor
   //printf("finestra 0 \n");

   for (i = 0; i < k; i++)
   {
      if (seq[i] == 'N')
	 nextNpos = i + k + 1;  
      next = next + ACTGnumber[(int)seq[i]] * pot4(i); 
   }
   last = next; lastNpos = nextNpos - k;

   // for each position of the first window.
   for (i = k; i <= w; i++)
   { 
      lastbin = -1;
      
      if (i >= nextNpos)
	// for each occurrence of the k factor called 'next'
	for (j = k_factor_ind->ind[next]; j < k_factor_ind->ind[next + 1]; j++)
	  // comment todo
	  CountBins((i - k) - k_factor_ind->tuilist[j] + N - 1, e, p, bins, &l); 

      // jump the 'N' positions
      if (seq[i] == 'N')
	 nextNpos = i + k + 1;

      // compute the next k-factor value (tricks on the bit shift for the computation)
      next = (next >> 2) + (pot4(k - 1) * ACTGnumber[(int)seq[i]]);
     
   } 
   //   printf("finestra 0 \n");
   // if the first window has enougth 'friends' (zone respecting the necessary conditions)
   if (Nelems(&begin, &end, l, 0, w, bin_size, N, goodWindows, curr_pass_num, curr_parall_min, curr_parall_max) >= r)
   { 
      begin = 0; end = 0; l2 = NewTree(); 
      // the current par. is added in the tree of good parallelograms
      //      printf("goodPar\n");
      goodPar(l, nextBegin, pLCS, 0, w, p, w, N, k, bin_size + e, seq, &l2, k_factor_ind);

      if (Nelems(&begin, &end, l2, 0, w, bin_size, N, goodWindows, curr_pass_num, curr_parall_min, curr_parall_max) >= r)
	ReportPgram(0, w, &result, goodWindows, curr_parall_min, curr_parall_max, curr_pass_num);
      /*EMPTY BLOCK + MULTI PASS*/
      //reset of goodWindows[curr_parall] from 1 to 0 during passes
      // successive to the first pass, this reset is done only if 
      // the notEmpty field was modified during a previous pass with
      // respect to current pass
      else
	for(curr_parall = curr_parall_min; curr_parall <=curr_parall_max; curr_parall++)
	  if(((goodWindows[curr_parall]).pass_num<curr_pass_num) && ((goodWindows[curr_parall]).notEmpty==1)){
	    (goodWindows[curr_parall]).notEmpty=0;
	    (goodWindows[curr_parall]).pass_num=curr_pass_num;
	  }
      /**/
      FreeTree(l2);
   }
   /*EMPTY BLOCK + MULTI PASS*/
   //reset of goodWindows[curr_parall] from 1 to 0 during passes
   // successive to the first pass, this reset is done only if
   // the notEmpty field was modified during a previous pass with
   // respect to current pass
   else
   {
     //for(curr_parall = curr_parall_min; curr_parall <=curr_parall_max; curr_parall++)
     //       printf("empty block else %d\n",(curr_parall));
       if(((goodWindows[curr_parall]).pass_num<curr_pass_num) && ((goodWindows[curr_parall]).notEmpty==1)){
	 (goodWindows[curr_parall]).notEmpty=0;
	 (goodWindows[curr_parall]).pass_num=curr_pass_num;
       }    

   }
   /**/
   
   /* Adjusting the value for the window sliding */
   nextNpos = nextNpos - w - 1;
   
   /* Andando com a janela de tamanho w*/
   for (i = 0; i < N - w; i++)
   {
     /* EMPTY BLOCK */
     //curr_parall is the index of the block containing the current window 
     curr_parall = (i+1) >> tuiuiu_z;
     curr_parall_min = curr_parall;
     curr_parall_max = curr_parall;
     if((curr_parall>0) && ((i+1) <= ((curr_parall * bin_size) + e)))
       curr_parall_min = curr_parall - 1;
     /**/
     
     //     printf("\r %d   ", i); 
      /* Contabilizando o k-factor que esta entrando na janela */
      lastbin = -1;
      if (i >= nextNpos)
         for (j = k_factor_ind->ind[next]; j < k_factor_ind->ind[next + 1]; j++)
	   CountBins(i + w - k - k_factor_ind->tuilist[j] + N, e, p, bins, &l);
      //else printf("%d ", i + w - k);
      if (seq[i + 1 + w] == 'N')
	 nextNpos = i + k + 1;
      next = (next >> 2) + (pot4(k - 1) * ACTGnumber[(int)seq[i + 1 + w]]); 
             
      /* Contabilizando o k-fator que esta saindo da janela */
      lastbin = -1;
      if (i >= lastNpos)
	//if(l!=NULL)
	  for (j = k_factor_ind->ind[last]; j < k_factor_ind->ind[last + 1]; j++){
            UncountBins(i - k_factor_ind->tuilist[j] + N - 1, e, p, bins, &l);
	}
      //else printf("%d ", i);
      if (seq[i + k] == 'N')
	 lastNpos = i + k + 1;
        
      last = (last >> 2) + (pot4(k - 1) * ACTGnumber[(int)seq[i + k]]);  
      
      /* Verificando a codicao do alinhamento multiplo */
      begin = 0; end = 0;
      //int x,y;
      int temp = Nelems(&begin, &end, l, i + 1, w, bin_size, N, goodWindows, curr_pass_num, curr_parall_min, curr_parall_max);
      //      printf("  %d   ", temp);	/* DEB */
      if (temp >= r)
      { 
	//	nb_tests_good_par++;	/* DEB */
	 begin = 0; end = 0; l2 = NewTree(); 
         goodPar(l, nextBegin, pLCS, i + 1, i + 1 + w, p, w, N, k, bin_size + e, seq, &l2, k_factor_ind);
         if (Nelems(&begin, &end, l2, i + 1, w, bin_size, N, goodWindows, curr_pass_num, curr_parall_min, curr_parall_max) >= r)
	   ReportPgram(i + 1, i + 1 + w, &result, goodWindows, curr_parall_min, curr_parall_max, curr_pass_num);
	 /*EMPTY BLOCK + MULTI PASS*/
	 //reset of goodWindows[curr_parall] from 1 to 0 during passes
	 // successive to the first pass, this reset is done only if 
	 // the notEmpty field was modified during a previous pass with
	 // respect to current pass
	 else
	   for(curr_parall = curr_parall_min; curr_parall <= curr_parall_max; curr_parall++)
	     if(((goodWindows[curr_parall]).pass_num<curr_pass_num) && ((goodWindows[curr_parall]).notEmpty==1)){
	       (goodWindows[curr_parall]).notEmpty=0;
	       (goodWindows[curr_parall]).pass_num=curr_pass_num;
	     }
	 /**/
	 
         FreeTree(l2);
	 //printf ("good %5d %5d  %5d %5d   %3d %3d\n", i+1, i+1+w, begin, end-begin+1, x, y);
      }
      /*EMPTY BLOCK + MULTI PASS*/
      //reset of goodWindows[curr_parall] from 1 to 0 during passes
      // successive to the first pass, this reset is done only if
      // the notEmpty field was modified during a previous pass with
      // respect to current pass
      else
	for(curr_parall = curr_parall_min; curr_parall <= curr_parall_max; curr_parall++)
	  if(((goodWindows[curr_parall]).pass_num<curr_pass_num) && ((goodWindows[curr_parall]).notEmpty==1)){
	    (goodWindows[curr_parall]).notEmpty=0;
	    (goodWindows[curr_parall]).pass_num=curr_pass_num;
	  }
      /**/
      //printf("??? %d\n", i);  
   }
      
   //   printf("nb_tests_good_par = %d\n", nb_tests_good_par); /* DEB */
   free(bins);
   FreeTree(l);
   free(nextBegin);
   free(pLCS);
   free(occurrencesList);
   free(listLSI);
   
   return result;
}      
