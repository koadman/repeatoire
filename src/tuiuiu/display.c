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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h> 
#include <getopt.h>
#include <time.h>
#include <string.h>
#include "tuilist.h"

////////////////////////////////////////////////////////////////////////////
///                                         multi                      ///// 
////////////////////////////////////////////////////////////////////////////    

double countResultMultiSeqs(int N, int *seqBegins, tuilist **result, int nseq)
{
  int total_kept = 0;
  int seqNumber=0;

  // fprintf(output, "\tPositions conserved by the filter:\n");

  for(seqNumber=0; seqNumber<nseq; seqNumber++){
    //fprintf(output, "> %s (Conserved after Filtration by Tuiuiu)\n", name[seqNumber]);
    tuilist * current_seq = result[seqNumber];
    while(current_seq !=NULL){
      total_kept += current_seq->end - current_seq->begin +1;
      //fprintf(output, "[%d-%d] ", current_seq->begin, current_seq->end);
      current_seq = current_seq->prox;
    }
    //fprintf(output, "\n");
  }

  return (double)total_kept / N;
}


/** 
 * output the result, respecting the multisequences format
 * @returns the percentage of conseved data
 **/
double WriteResultMultiSeqs(int N, int *seqBegins, tuilist **result, FILE *output,  int nseq, char **name)
{
 
   int total_kept = 0;   
   int seqNumber=0;

   fprintf(output, "\tPositions conserved by the filter:\n");
   
   for(seqNumber=0; seqNumber<nseq; seqNumber++){
     fprintf(output, "> %s (Conserved after Filtration by Tuiuiu)\n", name[seqNumber]); 
     tuilist * current_seq = result[seqNumber];
     while(current_seq !=NULL){
       total_kept += current_seq->end - current_seq->begin +1;
       fprintf(output, "[%d-%d] ", current_seq->begin, current_seq->end);
       current_seq = current_seq->prox;
     }
     fprintf(output, "\n");
   }
   
   return (double)total_kept / N; 
}

void saveFilteredMultiSeqWithNs(int *seqBegins, tuilist **result, char **seq, int nseq)
{
    int i = 0,seqNumber = 0;
     
   for(seqNumber=0;seqNumber<nseq;seqNumber++){
     //fprintf(output, "> %s (Conserved after Filtration by Tuiuiu (removed part are masked by 'N's)\n", name[seqNumber]); 
     tuilist * current_seq = result[seqNumber];
     int prev=0;
     while(current_seq !=NULL){
       for(i=prev;i<current_seq->begin;i++) {
	 (*seq)[i+seqBegins[seqNumber]] = 'N';
	 //fprintf(output, "N");
	 //if (++j % 100==0) fprintf(output, "\n");
       }
       
       for(i= current_seq->begin; i<= current_seq->end; i++){
	 //fprintf(output, "%c", seq[i+seqBegins[seqNumber]]);
	 //if (++j % 100==0) fprintf(output, "\n");
       }
       prev=current_seq->end+1;
       current_seq = current_seq->prox;
     }
     int size_seq = seqBegins[seqNumber+1]-seqBegins[seqNumber];
     for(i=prev;i<size_seq;i++) {
       (*seq)[i+seqBegins[seqNumber]] = 'N';
       //fprintf(output, "N");
       //if (++j % 100==0) fprintf(output, "\n");
     }
     //fprintf(output, "\n");
   }
}

/**
 * print data of the conserved regioins.
 * 
 **/
void WriteFilteredMultiSeqWithNs(int *seqBegins,  tuilist **result, char *seq, FILE *output, int nseq, char **name, int bin_size)
{
    int i = 0,j =0, seqNumber = 0;
     
   for(seqNumber=0;seqNumber<nseq;seqNumber++){
     j=0;
     fprintf(output, "> %s (Conserved after Filtration by Tuiuiu (removed part are masked by 'N's)\n", name[seqNumber]); 
     tuilist * current_seq = result[seqNumber];
     int prev=0;
     while(current_seq !=NULL){
       for(i=prev;i<current_seq->begin;i++) {
	 fprintf(output, "N");
	 if (++j % 100==0) fprintf(output, "\n");
       }
       
       for(i= current_seq->begin; i<= current_seq->end; i++){
	 fprintf(output, "%c", seq[i+seqBegins[seqNumber]]);
	 if (++j % 100==0) fprintf(output, "\n");
       }
       prev=current_seq->end+1;
       current_seq = current_seq->prox;
     }
     int size_seq = seqBegins[seqNumber+1]- bin_size -seqBegins[seqNumber];
     for(i=prev;i<size_seq;i++) {
       fprintf(output, "N");
       if (++j % 100==0) fprintf(output, "\n");
     }
     fprintf(output, "\n");
   }
}



/**
 * print data of the conserved regioins.
 * 
 **/
void WriteFilteredMultiSeqOneN(int *seqBegins, tuilist **result, char *seq, FILE *output, int nseq, char **name)
{
   int i = 0,j =0, seqNumber = 0;
     
   for(seqNumber=0;seqNumber<nseq;seqNumber++){
     j=0;
     
     fprintf(output, "> %s (Conserved after Filtration by Tuiuiu - concatenated and separated by one 'N')\n", name[seqNumber]); 
     tuilist * current_seq = result[seqNumber];
     while(current_seq !=NULL){
       fprintf(output, "N"); 
       if (++j % 100==0) fprintf(output, "\n");

       for(i= current_seq->begin; i<= current_seq->end; i++){
	 fprintf(output, "%c", seq[i+seqBegins[seqNumber]]);
	 if (++j % 100==0) fprintf(output, "\n");
       }
       current_seq = current_seq->prox;
     }
     fprintf(output, "\n");
   }

}



/**
 * print data of the conserved regioins.
 * 
 **/
void WriteFilteredMultiSeqConcatenated(int *seqBegins, tuilist **result, char *seq, FILE *output, int nseq, char **name)
{
   int i = 0,j =0, seqNumber = 0;
     
   for(seqNumber=0;seqNumber<nseq;seqNumber++){
     j=0;
     
     fprintf(output, "> %s (Conserved after Filtration by Tuiuiu - concatenated)\n", name[seqNumber]); 
     tuilist * current_seq = result[seqNumber];
     while(current_seq !=NULL){
       for(i= current_seq->begin; i<= current_seq->end; i++){
	 fprintf(output, "%c", seq[i+seqBegins[seqNumber]]);
	 if (++j % 100==0) fprintf(output, "\n");
       }
       current_seq = current_seq->prox;
     }
     fprintf(output, "\n");
   }
}



/**
 * print data of the conserved regioins.
 * 
 **/
void WriteFilteredMultiSeq(int *seqBegins, tuilist **result, char *seq, FILE *output, int nseq, char **name, int flag, int bin_size)
{
  if(flag==0) return WriteFilteredMultiSeqConcatenated(seqBegins, result, seq, output, nseq, name);
  if(flag==1) return WriteFilteredMultiSeqOneN(seqBegins, result, seq, output, nseq, name);
  if(flag==2) return WriteFilteredMultiSeqWithNs(seqBegins, result, seq, output, nseq, name, bin_size);
  
  
  fprintf(stderr, "bad option value -N %d, write  Ns instead of filtered nucleotides in result file\n", flag); 
  WriteFilteredMultiSeqWithNs(seqBegins, result, seq, output, nseq, name, bin_size);
}






////////////////////////////////////////////////////////////////////////////
///                                         mono                       ///// 
////////////////////////////////////////////////////////////////////////////    


double countResult(int N, tuilist *result)
{
  int total_kept = 0;
  //fprintf(output, "\tPositions conserved by the filter:\n\t\t");
  while (result != NULL)
    {
      //fprintf(output, "[%d-%d] ", result->begin, result->end);
      total_kept = total_kept + (result->end - result->begin) + 1;
      result = result->prox;
    }
  //   fprintf(output,"\n\nTotal Kept: %f\n", (double)total_kept / N);
  //fprintf(output,"\n");
  return (double)total_kept / N;
}



double WriteResult(int N, tuilist *result, FILE *output)
{
   int total_kept = 0;   
   fprintf(output, "\tPositions conserved by the filter:\n\t\t");
   while (result != NULL)
   {
      fprintf(output, "[%d-%d] ", result->begin, result->end);
      total_kept = total_kept + (result->end - result->begin) + 1;
      result = result->prox;
   }
   //   fprintf(output,"\n\nTotal Kept: %f\n", (double)total_kept / N);
   fprintf(output,"\n");
   return (double)total_kept / N;
}   


void saveFilteredSeqWithNs(int N, char **seq, tuilist *result){
  int i;
   
  //fprintf(output, "> %s (Filtered by tuiuiu)\n", name);  
  //j=0;
  for (i = 0; i < N; i++)
    {
      if (result != NULL && i == result->begin)
	{
	  for (; i <= result->end; i++){
	    //fprintf(output, "%c", seq[i]);
	    //if (++j%100==0) fprintf(output, "\n");
	  }
	  i--;
	  result = result->prox;
	} 
      else {
	//fprintf(output, "N");
	(*seq)[i] = 'N';
	//if (++j%100==0) fprintf(output, "\n");
      }
    }
  //fprintf(output, "\n");
}


void WriteFilteredSeqWithNs(int N, char *seq, tuilist *result, FILE *output, char *name){
  int i,j;
   
  fprintf(output, "> %s (Filtered by tuiuiu)\n", name);  
  j=0;
  for (i = 0; i < N; i++)
    {
      if (result != NULL && i == result->begin)
	{
	  for (; i <= result->end; i++){
	    fprintf(output, "%c", seq[i]);
	    if (++j%100==0) fprintf(output, "\n");
	  }
	  i--;
	  result = result->prox;
	} 
      else {
	fprintf(output, "N");
	if (++j%100==0) fprintf(output, "\n");
      }
    }
  fprintf(output, "\n");
}
  

void WriteFilteredSeqWithOneN(int N, char *seq, tuilist *result, FILE *output, char *name){
  int i,j;
   
  fprintf(output, "> %s (Filtered by tuiuiu)\n", name);  
  j=0;
  int inter=0;
  for (i = 0; i < N; i++)
    {
      if (result != NULL && i == result->begin)
	{
	  inter=0;
	  for (; i <= result->end; i++){
	    fprintf(output, "%c", seq[i]);
	    if (++j%100==0) fprintf(output, "\n");
	  }
	  i--;
	  result = result->prox;
	} 
      else {
	if(inter==0){
	  fprintf(output, "N");
	  if (++j%100==0) fprintf(output, "\n");
	  inter=1;
	}
      }
    }
  fprintf(output, "\n");
}
  

void WriteFilteredSeqConcatenated(int N, char *seq, tuilist *result, FILE *output, char *name){
  int i,j;
   
  fprintf(output, "> %s (Filtered by tuiuiu)\n", name);  
  j=0;
  for (i = 0; i < N; i++)
    {
      if (result != NULL && i == result->begin)
	{
	  for (; i <= result->end; i++){
	    fprintf(output, "%c", seq[i]);
	    if (++j%100==0) fprintf(output, "\n");
	  }
	  i--;
	  result = result->prox;
	} 
      else ;
    }
  
  fprintf(output, "\n");
}
  

void WriteFilteredSeq(int N, char *seq, tuilist *result, FILE *output, char *name, int flag)
{
  if(flag==0) return WriteFilteredSeqConcatenated(N, seq, result, output, name);
  if(flag==1) return WriteFilteredSeqWithOneN(N, seq, result, output, name);
  if(flag==2) return WriteFilteredSeqWithNs(N, seq, result, output, name);
  
  fprintf(stderr, "bad option value -N %d, write  Ns instead of filtered nucleotides in result file\n", flag); 
  WriteFilteredSeqWithNs(N, seq, result, output, name);
}

