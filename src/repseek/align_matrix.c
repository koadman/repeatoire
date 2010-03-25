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
        file     : matrix.c
        function : create a log(freq) / identity matrix for DNA
                                                 
        created  : 15 Oct 2002
        modif    : nov 12 2002: Set all score of X against anything to -1000; (usefull for mask and twoseq)
		  modif    : nov 26 2002: get my new matrix format PtrConst->pScoring (for removing seqaln and my S&W)
		  modif    : Oct 2003  
		  modif    : July 2004  : change the formula of the log matrix, mainly Match(N,i) is changed
		  modif    : July 2004  : turn scores into doubles.
		  
        author   : amikezor
*****/

#include "repseek_types.h"
#include "align.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/***
	calculate the frequency of all symbol
		in DNA sequence
***/
static int32_t get_letter(int32_t letters[5], char *seq)
{

  int32_t lseq=0;
  int8_t letter;
   
	letters[0] = letters[1] = letters[2] = letters[3] = letters[4] = 0;
	
  while( (letter = *seq) != 0)
  {

    lseq++;

    if (letter == 'A')
      letters[0]++;
    if (letter == 'C')
      letters[1]++;
    if (letter == 'G')
      letters[2]++;
    if (letter == 'T')
      letters[3]++;
    
    seq++;
  }
  
  letters[4] = lseq-letters[0]-letters[1]-letters[2]-letters[3]-letters[4];

  return lseq;
}

/***
	From a float return the closest int
***/
static int float2int(float a){

	if(a<0){
		if( ( a-(float)(int)a ) > -0.5 )
			return (int)a;
		else
			return ((int)a-1);
	}
	else{
		if( ( a-(float)(int)a ) < 0.5 )
			return (int)a;
		else
			return ((int)a+1);
	}
		
}


/*
	Print the matrix if you want
*/
void printout_matrix2( SCORING *pScoring ){

	int i,j;

	fprintf(stderr,"    A    B    C    D    E    F    G    H    I    J    K    L    M    N    O    P"
	"    Q    R    S    T    U    V    W    X    Y    Z\n");
	for (i=0; i<26; i++) {
		fprintf(stderr,"%c ", SYMB2CHAR(i) );
		for (j=0; j<26; j++)
			if(pScoring->matrix[i][j]<100 && pScoring->matrix[i][j]>-100)
				fprintf(stderr,"%5.2f",pScoring->matrix[i][j]);
			else
				fprintf(stderr,"%5.1g",pScoring->matrix[i][j]);

		fprintf(stderr,"\n");
	}
}


/*
	Print the matrix if you want
*/
void printout_matrix( SCORING *pScoring ){

	int i,j;
	char carac[6]={'A', 'C', 'G', 'T', 'N', 'X'};

	fprintf(stderr,"     ");
	for (i=0; i<6; i++)
		fprintf(stderr,"%6c", carac[i]);
	fprintf(stderr,"\n");
	
	for (i=0; i<6; i++) {
	
		fprintf(stderr,"%-6c ", carac[i] );
		
		for (j=0; j<6; j++)
			if(pScoring->matrix[CHAR2SYMB(carac[i])][CHAR2SYMB(carac[j])]<100 && 
			   pScoring->matrix[CHAR2SYMB(carac[i])][CHAR2SYMB(carac[j])]>-100)
				fprintf(stderr,"%6.2f",pScoring->matrix[CHAR2SYMB(carac[i])][CHAR2SYMB(carac[j])]);
			else
				fprintf(stderr,"%6.1g",pScoring->matrix[CHAR2SYMB(carac[i])][CHAR2SYMB(carac[j])]);

		fprintf(stderr,"\n");
	}
}



/***
	Set up the matrix MATCH(i,j) = 1/2 * log4(pi*pj) - thought to be the more corrective
	get all symbol frequency and estimate match(N,i) = pi * log(pi);
***/
void log_matrix(char *sequence, float gap_open, float gap_ext, SCORING *pScoring){

	int32_t letters[5]={0,};           /* count of every letter */
	double freq[4]= {0.0,};            /* frequency of nucleotides ACGT */
	char carac[4]={'A','C','G','T'};

	int32_t size;                      /* size of the seq including Ns ans Xs */
	double temp;                       /* make the read easier */
 	int32_t size_eff;                  /* size_eff is the real number of non N nucleotide */
	int16_t i=0, j=0, u;               /* dummy counters */

        //punt: segmentation fault here if pScoring not allocated before passing reference
	for(i=0;i<26;i++)
	{
		for(j=0;j<26;j++)
		{
                  //fprintf(stderr,"");
		  pScoring->matrix[i][j]=0;
		}
	}
        //fprintf(stderr,"get_letter\n");
	size = get_letter(letters, sequence);             /* count all letters */

	size_eff = size-letters[4];                       /* count the size of ACGT symbols */
	for(i=0;i<4;i++)                                  /* get the frequencies */
		freq[i]=(double)letters[i]/(double)size_eff;

	//fprintf(stderr,"set matrix score\n");
	/*
		Set the matrix score
	*/
	for(i=0;i<4;i++){
		for(j=0;j<4;j++){ 
			temp = 0.5 * log10( freq[i]*freq[j] )/log10( 4.0 );
			pScoring->matrix[ CHAR2SYMB(carac[i]) ][ CHAR2SYMB(carac[j]) ] = (i==j)? -temp: temp;
		}
		temp = -  freq[i] * log10( freq[i] )/log10( 4.00 );
		pScoring->matrix[CHAR2SYMB('N')][CHAR2SYMB(carac[i])] = temp;
		pScoring->matrix[CHAR2SYMB(carac[i])][CHAR2SYMB('N')] = temp;

	}

	pScoring->matrix[CHAR2SYMB('N')][CHAR2SYMB('N')] = 0.25;
	
	 /*
	 	set all X against anything to -1e6
	 */
	for(i=0;i<26;i++)             
		pScoring->matrix[i][CHAR2SYMB('X')] = pScoring->matrix[CHAR2SYMB('X')][i] = -1e6; 


	/*
		Estimate the mean score of the matrix
	*/
 	pScoring->expect=0.0;
	for(i=0;i<4;i++)
	{
		for(u=0;u<4;u++)             /* print out log(freq*freq) into the pScoring->matrix */
		{
			pScoring->expect +=  freq[i]*freq[u]*pScoring->matrix[CHAR2SYMB(carac[i])][CHAR2SYMB(carac[u])];
		}
	}

	pScoring->gap_open = (gap_open<0)? gap_open : -gap_open;
	pScoring->gap_ext  = (gap_ext<0)?  gap_ext  : -gap_ext;


	/*
		printout_matrix( pScoring );

	*/
}

/*
	This matrix is
		anything is +/- 1
		N is 1
		X is -1e6
*/
void identity_matrix(char *sequence, float gap_open, float gap_ext,SCORING *pScoring){

	int16_t i,j,u;
	int32_t letters[5]={0,};
	double freq[4]= {0.0,};          /* frequency of nucleotides ACGT */

	char carac[4] = {'A','C','G','T'};

	int32_t size;
	
	int32_t size_eff;                /* size_eff is the real number of non N nucleotide */
	

	for(i=0;i<26;i++)                                 /* set the matrix to 0 */
		for(j=0;j<26;j++)
			pScoring->matrix[i][j]=0;


	size = get_letter(letters, sequence);             /* count all letters */
	size_eff = size-letters[4];                       /* count the size of ACGT symbols */

	for(i=0;i<4;i++)                                  /* get the frequencies */
		freq[i]=(double)letters[i]/(double)size_eff;
		
		
	for(i=0;i<4;i++)
		for(u=0;u<4;u++)
			pScoring->matrix[CHAR2SYMB(carac[i])][CHAR2SYMB(carac[u])]=(u==i)? 1 : -1;

	for(i=0;i<26;i++)              /* set all N against anything to +1 */
			pScoring->matrix[i][CHAR2SYMB('N')] = pScoring->matrix[CHAR2SYMB('N')][i] = 0.25; 

	for(i=0;i<26;i++)              /* set all X against anything to -1e6 */
			pScoring->matrix[i][CHAR2SYMB('X')] = pScoring->matrix[CHAR2SYMB('X')][i] = -1e6; 


 	pScoring->expect=0.0;
	for(i=0;i<4;i++)
		for(u=0;u<4;u++){             /* print out log(freq*freq) into the pScoring->matrix */
			pScoring->expect +=  freq[i]*freq[u]*pScoring->matrix[CHAR2SYMB(carac[i])][CHAR2SYMB(carac[u])];
	}
	

	pScoring->gap_open = (gap_open>0)? -gap_open : gap_open;
	pScoring->gap_ext  = (gap_ext>0 )? -gap_ext  : gap_ext;

	/*
	printout_matrix(  pScoring );
	*/	

}
