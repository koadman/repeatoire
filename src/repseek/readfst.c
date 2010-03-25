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
        file     : readfst.c
        function : what you need to read the fst files

        created  : Oct 03 2003
        modif    : may 2008 - count Xs in the sequence

			note    : copy and paste from libiofst, and change the malloc() to MyMalloc()
        author   : amikezor
*****/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "repseek_types.h"
#include "sequence.h"
#include "memory.h"


static int32_t size_oneseq(FILE *fseq)
{
  int32_t size, count=0, bonux=0,      /* size of the seq, count number of letter/line, and final carac(\n or nothing) */
     begin_oct, end_oct;               /* begin and and end in octet */

  int temp;                            /* used for temp memory */

  while( fgetc(fseq) != '\n');         /* get out the first line [comments] */

  begin_oct = ftell(fseq);             /* where begin the second line */

  while( fgetc(fseq) != '\n')
    count++;                           /* get the size of the second line */


	fseek(fseq, -1, SEEK_END);         /* Do we have ACGT\n[EOF] or ACGT[EOF]   */
	bonux = (fgetc(fseq) == '\n')?1:0;
  
   fseek(fseq,-count-bonux,SEEK_END);

  while( fgetc(fseq) != '\n');

  end_oct = ftell(fseq);                /* where finish the ante-last line */

  size = (end_oct-begin_oct)/(count+1);
  size*=count;
   
  count = 0;
  while((temp=fgetc(fseq)) != EOF && temp != '\n' )
    count++;                             /* get the last line */

  size+=count;

  fseek(fseq,0,SEEK_SET);                /* go back to the begin of the file */

  return size;
}


static int32_t get_curfst(char *sequence, FILE *fseq)
{

	int32_t lseq;
	int letter;
  
	lseq = 0;
  
	while(fgetc(fseq) != '\n');  /* lecture de la premiere ligne */


	while ((letter=fgetc(fseq)) != EOF || letter== '>')
		if (letter !='\n')
			sequence[lseq++]=letter;

	if(letter == '>')
		fseek(fseq, -1, SEEK_CUR);

  sequence[lseq]=0;
  return lseq;
  
}

static int32_t numseq(FILE *ffile)
{

	int32_t number=0;
	int tmp;


  if(fgetc(ffile) != '>')
    fprintf(stderr,"Check your fasta sequences\n"), exit(4);
 
	rewind(ffile);                            /* go back at the begining */

	while( ( tmp=fgetc(ffile) ) != EOF )
		if(tmp == '>'){
			number++;
			while(fgetc(ffile) != '\n');
		}
	
	
	rewind(ffile);  /*return at the file's begining*/

	return number;
}




char *file2mem(char *file, int32_t *size)  /*le pointeur size permet de recup la size */
{

  FILE *fseq;
  char *seq;

  fseq = fopen(file,"r");
  if(fseq == NULL)
    fprintf(stderr,"Flux error, check the fasta seq file\n"),exit(3);


	if( numseq(fseq) == 1)
	  *size = size_oneseq(fseq);                           /* return length of seq */
	else
		fprintf(stderr,"there is not exactly ONE sequence in that file\n"),exit(2);

  seq = (char *)MyMalloc((*size+1)*sizeof(char),"get_chr: Malloc error");

  get_curfst(seq, fseq); 

  fclose(fseq);

  return seq;
}


char *fuse2files_mem(char *file1, char *file2, int32_t *size1, int32_t *size2)
{

	FILE *fseq1;    /* pointer to the first file */
	FILE *fseq2;    /* pointer to the second one */
	char *seq;      /* the sequence seq1_X_seq2 */

	int32_t size;       /* the total size = size1+size2+1 */
	
	fseq1 = fopen(file1,"r");
	fseq2 = fopen(file2,"r");
	if(!fseq1 ||  !file2)fprintf(stderr,"twofiles2mem: Flux error, check the fasta seq files\n"),exit(3);

	if( numseq(fseq1) == 1)
	  *size1 = size_oneseq(fseq1);                           /* return length of seq1 */
	else
		fprintf(stderr,"twofiles2mem: there is not only 1 sequence in file %s\n", file1),exit(2);


	if( numseq(fseq2) == 1)
	  *size2 = size_oneseq(fseq2);                           /* return length of seq2 */
	else
		fprintf(stderr,"twofiles2mem: there is not only 1 sequence in file %s\n", file2),exit(2);

	size = (*size1) + 1 + (*size2);

	seq = (char *)MyMalloc((size+1)*sizeof(char),"twofiles2mem: malloc error");

	if (get_curfst(seq, fseq1) != *size1)
	fprintf(stderr,"twofiles2mem: sequence of file %s is not clean\n",file1),exit(4);

	seq[ *size1 ]='X';
	
	if (get_curfst(seq + (*size1)+1, fseq2) != *size2)
	fprintf(stderr,"twofiles2mem: sequence of file %s is not clean\n",file2),exit(4);
	
	fclose(fseq1);
	fclose(fseq2);

	return seq;
}

char *readFastaSeq(const char *file, int32_t *size, int32_t *nX)
{
	
  char *sequence;
  
  /*fprintf(stderr,"get sequence... %s   ",file);*/
  
  sequence = file2mem(file, size); 
  (void) UpperSequence(sequence);                        /* turn into UPPER sequence */
  nonACGTXtoN(sequence);                                 /* set all non ACGTX into N */
  (void )CheckSeq(sequence, "ACGTX");                    /* check the number of bad symbols */
  
  *nX=countX(sequence);
  

  return sequence;
}

