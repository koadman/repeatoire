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
//#include <getopt.h>
#include <time.h>
#include <string.h>
#include "mono_and_multi_commons.h"

inline void CountBins(int d, int e, int p, int bins[], itree **l)
{
  int bin, binstop;
  extern int tuiuiu_z;
   for ( bin = d >> tuiuiu_z, binstop = (d + e) >> tuiuiu_z; bin <= binstop; bin ++ ){
     
     if (bin > lastbin) {
     
       bins[bin]++;
       lastbin = bin;
       if (bins[bin] == p){	 
	   AddTree(bin, l);
       }
	 
     }
   }
}


inline void UncountBins(int d, int e, int p, int bins[], itree **l)
{
   int bin, binstop;
   extern int tuiuiu_z;
   for ( bin = d >> tuiuiu_z, binstop = (d + e) >> tuiuiu_z; bin <= binstop; bin ++ ){

     if (bin > lastbin) {
       bins[bin]--;
       lastbin = bin;
       if (bins[bin] == p - 1)
	 RemoveTree(bin, l);	   
     }
   }
}


void ReportPgram(int begin, int end, tuilist **result, empty_block goodWindows[], int curr_parall_min, int curr_parall_max, int curr_pass_num)
{
  int curr_parall;
  /* EMPTY BLOCK */
  // goodWindows update when the current window is kept
  for(curr_parall = curr_parall_min; curr_parall <=curr_parall_max; curr_parall++){
    (goodWindows[curr_parall]).notEmpty=1;
    (goodWindows[curr_parall]).pass_num=curr_pass_num;
  }
  /**/

  end--;
  if (*result != NULL)
   {
      if ((*result)->end >= begin)
      {
	 (*result)->end = end;
         return;
      }
   } 
   *result = Add(*result, begin, end);
}    


/*
 * returns the time in seconds since the begining of the processuss
 */
inline double cpuTime()
{
  return clock() / (double) CLOCKS_PER_SEC; //time(NULL);
}

/*
  set position just after the first symbol read
  @returns the first non space symbol read
*/
char firstNonSpaceRead(FILE *file){
  char base;
  for(base = (char) fgetc(file); base == ' ' || base == '\n' || base=='\t'; base = (char) fgetc(file)) 
    if(base== EOF){
      return EOF;
    }
  return base;
}


/**
  Read all lines starting with > (allowing spaces before the >
  Set the file offset at the beginning of the first line that doesn't start with a >
  Fill the variable name with the value of the comment read. 
  If no comment is read, file name with the defaultmessage value
  @returns 1 if at least one line starting with > was read, else return 0
**/
int readcomments(FILE *file, char **name, char * default_message){
  long before_comments = ftell(file);
  int sizename=0;
  const int chunk=128;
  
       
  char base = firstNonSpaceRead(file);
  if(base == EOF) return 0;
  
  if(base!='>'){		/* no comment :-) */
    *name = (char *) malloc(128*sizeof(char));
    strcpy (*name, default_message);
    fseek(file, before_comments, SEEK_SET); /* set the offset at the beggining of the line */
    return 0;
  }

  // if we are here, at least one comment was read

  while(1){
    for ( base = (char) fgetc(file);
	  base != '\n';
	  base = (char) fgetc(file) ) {
      if(base == EOF) break;
      if (sizename % chunk == 0)
	*name = (char *) realloc( *name, (sizename+chunk)*sizeof(char) );
      (*name)[sizename] = base;
      sizename++;
    }

    // check next line
    before_comments = ftell(file);
    base = firstNonSpaceRead(file);
    if(base!='>') {			      /* the real sequence starts */
      fseek(file, before_comments, SEEK_SET); /* set the offset at the beggining of the line */
      break;
    }
  }
  *name = (char *) realloc( *name, (sizename+1)*sizeof(char) );
  (*name)[sizename] = '\0';
  
  return 1;
}

int readsequence(FILE *file, char**seq, char **name){
  char base;
  int sizeseq;
  const int chunk=128;
  
  sizeseq = 0;
  *seq = NULL;
  *name = NULL;

  readcomments(file, name, "Comment was missing, this one was generated by tuiuiu");
  
  
  for ( base = toupper((char) fgetc(file));
	base != '>' && base != EOF;
	base = toupper((char) fgetc(file)) )
    if ( base == 'A' || base == 'C' || base == 'T'
	 || base == 'G' || base == 'N') {  /* N means any base */
      if (sizeseq % chunk == 0)
	*seq = (char *) realloc( *seq, (sizeseq+chunk)*sizeof(char) );
      (*seq)[sizeseq] = base;
      sizeseq++;
    }
  
  //*seq = (char *) realloc( *seq, (sizeseq+1)*sizeof(char) );
  //(*seq)[sizeseq] = '\0';
  
  if (base!=EOF) ungetc(base,file);	/* to allow reading multifasta files */
  
  //if(sizeseq>0) printf("%s: %d bases\n", *name, sizeseq);

  return sizeseq;
}

inline char reverse (const char c){
  switch (c){
  case 'A':
  case 'a': return 'T';
  case 'C':
  case 'c': return 'G';
  case 'G':
  case 'g': return 'C';
  case 'T':
  case 't': return 'A';
  }
  return 'N';
}

// Read a sequence of length N. 
// returns a sequence of length 2N: 
//   - From 0 to N-1: direct sequence
//   - From N to 2N-1: reverse complemented sequence
// For instance if ATACC if given, this function returns ATACCGGTAT
int readsequenceAndReverse(FILE *file, char**seq, char **name){
  char base;
  int i, sizeseq;
  const int chunk=128;
  
  sizeseq = 0;
  *seq = NULL;
  *name = NULL;

  readcomments(file, name, "Comment was missing, this one was generated by tuiuiu");
  
  
  for ( base = toupper((char) fgetc(file));
	base != '>' && base != EOF;
	base = toupper((char) fgetc(file)) )
    if ( base == 'A' || base == 'C' || base == 'T'
	 || base == 'G' || base == 'N') {  /* N means any base */
      if (sizeseq % chunk == 0)
	*seq = (char *) realloc( *seq, (sizeseq+chunk)*sizeof(char) );
      (*seq)[sizeseq] = base;
      sizeseq++;
    }
  
  /// 
  // store the reverse complement of the sequence
  *seq = (char *) realloc( *seq, (2*sizeseq)*sizeof(char) );
  
  for (i=0;i<sizeseq;i++) (*seq)[i+sizeseq] = reverse ((*seq)[sizeseq-i-1]);

  sizeseq *=2;

  *seq = (char *) realloc( *seq, (sizeseq+1)*sizeof(char) );
  (*seq)[sizeseq] = '\0';
  
  if (base!=EOF) ungetc(base,file);	/* to allow reading multifasta files */
  
  // if(sizeseq>0) printf("%s: %d bases\n", *name, sizeseq);

  // printf("%s\n", *seq);		/* DEB */

  return sizeseq;
}

int reverseSequence(char**seq, int sizeseq){
  
  int i;
    
  /// 
  // store the reverse complement of the sequence
  *seq = (char *) realloc( *seq, (2*sizeseq)*sizeof(char) );
  
  for (i=0;i<sizeseq;i++) (*seq)[i+sizeseq] = reverse ((*seq)[sizeseq-i-1]);

  sizeseq *=2;

  *seq = (char *) realloc( *seq, (sizeseq+1)*sizeof(char) );
  (*seq)[sizeseq] = '\0';
  
  return sizeseq;
}


/**
  jump all lines starting with > (allowing spaces before the >
  @returns 1 if at least one line starting with > was read, else return 0
**/
int jumpcomments(FILE *file){
  long before_comments = ftell(file);
  
       
  char base = firstNonSpaceRead(file);
  if(base == EOF) return 0;
  
  if(base!='>'){		/* no comment :-) */
    fseek(file, before_comments, SEEK_SET); /* set the offset at the beggining of the line */
    return 0;
  }

  // if we are here, at least one comment was read

  while(1){
    for ( base = (char) fgetc(file);
	  base != '\n';
	  base = (char) fgetc(file) ) {
      if(base == EOF) break;
    }

    // check next line
    before_comments = ftell(file);
    base = firstNonSpaceRead(file);
    if(base!='>') {			      /* the real sequence starts */
      fseek(file, before_comments, SEEK_SET); /* set the offset at the beggining of the line */
      break;
    }
  }
  
  return 1;
}

/**
 * Returns 0 if no other sequence are after the current offset in file, else return the length of the sequence
 */
int isThereASequence (FILE *file){
  char base;
  int sizeseq;
  const long before_comments = ftell(file);
  
  sizeseq = 0;

  jumpcomments(file);
  

  for ( base = toupper((char) fgetc(file));
	base != '>' && base != EOF;
	base = toupper((char) fgetc(file)) )
    if ( base == 'A' || base == 'C' || base == 'T'
	 || base == 'G' || base == 'N') {  /* N means any base */
      sizeseq++;
    }
  
  if (base!=EOF) ungetc(base,file);	/* to allow reading multifasta files */
  fseek(file, before_comments, SEEK_SET);
  return sizeseq;
}

