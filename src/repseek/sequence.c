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
 * @file   sequence.c
 * @author amikezor
 * @date   20 Dec 2002
 * @modif  Oct 2003, April 2004, may 2008 - count Xs in the sequence
 * 
 * @brief  sequence manipulation
 * 
 * 
 */


#include <stdio.h>
#include <string.h>
#include "repseek_types.h"

/* -------------------------------------------- */
/*             check sequence symbols           */
/* -------------------------------------------- */
int8_t CheckSeq(char *seq, char *alpha)
{
	int32_t ntot, nbad;
	
	
	for ( ntot = nbad = 0; *seq ; seq++, ntot++)
		if (! strchr(alpha, *seq))
			nbad++;

	if (nbad){
		fprintf(stderr, "*WARNING* %d bad symbols on %d total ... ", nbad, ntot);
		return 0;
	}
	     
	return 1;
}

/* -------------------------------------------- */
/*   change all non A,C,G,T,X char into N       */
/* -------------------------------------------- */
void nonACGTXtoN(char *seq)
{
	char *pseq;
	
	for( pseq=seq; *pseq; pseq++ )
		if (! strchr("ACGTX", *pseq) )
			*pseq='N';

		 
	return;
}



/* -------------------------------------------- */
/*   count all X in a sequence                  */
/* -------------------------------------------- */
int32_t countX(char *seq)
{
	char *pseq;
	int32_t nX=0;
	
	for( pseq=seq; *pseq; pseq++ )
		if ( *pseq == 'X' )
			 nX++;
		 
	return nX;
}




/* -------------------------------------------- */
/*          uppercase sequence				      */
/* -------------------------------------------- */

#define IS_LOWER(c) (((c) >= 'a') && ((c) <= 'z'))
#define TO_UPPER(c) ((c) - 'a' + 'A')

void UpperSequence(char *seq){

    char *cseq;

	for (cseq = seq ; *cseq ; cseq++) 
		if (IS_LOWER(*cseq))
			*cseq = TO_UPPER(*cseq);
		 
}
 
#undef IS_LOWER
#undef TO_UPPER



void invseq(char *seqsrc, char *seqdest)
{

	char *pseqsrc,
		  *pseqinv;

	pseqsrc=seqsrc;
	pseqinv=seqdest+strlen(seqsrc)-1;

  while(*pseqsrc){
     	  
		if(*pseqsrc == 'A' || *pseqsrc == 'a'){*pseqinv = 'T';}
		else if(*pseqsrc == 'C' || *pseqsrc == 'c'){*pseqinv = 'G';}
		else if(*pseqsrc == 'G' || *pseqsrc == 'g'){*pseqinv = 'C';}
		else if(*pseqsrc == 'T' || *pseqsrc == 't'){*pseqinv = 'A';}
		else if(*pseqsrc == 'X' || *pseqsrc == 'x'){*pseqinv = 'X';}
	   else{*pseqinv = 'N';}
		
		pseqinv--;
		pseqsrc++;
		
    }
	

	seqdest[strlen(seqsrc)]=0;

}


