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
 * @file   smin.c
 * @author Guillaume Achaz <gachaz@gmail.com>
 * @date   Jan 2006
 * @modif  Mar 06 - a bug while plugin the regressions numbers-
 * @modif May 08 - do not count Xs from the sequence
 *
 * @brief  from a sequence and a pvalue, get an Smin
 * 
 * 
 */

#include "repseek_types.h"
#include <string.h>
#include <math.h>

static float get_dGC(char * sequence, int32_t len )
{
	int32_t letters[5] = {0,};
	int32_t i=0;

	for(i=0;i<len; i++){
	
		if (sequence[i] == 'A')
			letters[0]++;
		else if (sequence[i] == 'C')
			letters[1]++;
		else if (sequence[i] == 'G')
			letters[2]++;
		else if (sequence[i] == 'T')
			letters[3]++;
		else
			letters[4]++;
	}

	return ABS( 50.0 - 100.0*(float)( letters[1] + letters[2] ) / (float)( len - letters[4] ) );
}

static float get_dGC_2seq(char * sequence1, int32_t len1, char *sequence2, int32_t len2 )
{

	int32_t letters[5] = {0,};
	int32_t i=0;

	for(i=0;i<len1; i++){
	
		if (sequence1[i] == 'A')
			letters[0]++;
		else if (sequence1[i] == 'C')
			letters[1]++;
		else if (sequence1[i] == 'G')
			letters[2]++;
		else if (sequence1[i] == 'T')
			letters[3]++;
		else
			letters[4]++;
	}
	
	for(i=0;i<len2; i++){
	
		if (sequence2[i] == 'A')
			letters[0]++;
		else if (sequence2[i] == 'C')
			letters[1]++;
		else if (sequence2[i] == 'G')
			letters[2]++;
		else if (sequence2[i] == 'T')
			letters[3]++;
		else
			letters[4]++;
	}

	
	return ABS( 50.0 - 100.0*(float)( letters[1] + letters[2] ) / (float)( len1+len2 - letters[4] ) );
}




static double get_log_gamma(float dGC, int32_t len){

	double log_gamma=0.0;
	log_gamma = - 0.0029688*dGC*dGC + 0.0727745*dGC -  0.0309583*log(len)*log(len)  ;
	return log_gamma;
}
static double get_log_p(float dGC, int32_t len){

	double log_p=0.0;
	log_p = 0.0001789 *dGC*dGC - 0.002409*dGC + 0.006838*log(len) - 1.060;
	return log_p;
}

/*
	Using the formula of Watrerman & Vingron
	log(-log(F)) = log(gamma*m*n) + S*log(p)
	for unique sequence it becomes
	log(-log(F)) = log(gamma*n*(n-1)/2) + S*log(p)
*/

static double get_Smin(float dGC, int32_t len, double pvalue){

	double Smin=0.0;
	
	double log_log_minusF = 0.0,
	       log_gamma      = 0.0,
	       log_p          = 0.0;


	log_log_minusF = log( -log( 1.00 - pvalue ) );
	log_gamma = get_log_gamma( dGC, len );
	log_p = get_log_p( dGC, len );	
	
	Smin = log_log_minusF - log_gamma - log( (double)len * (len-1.0) / 2.0 );
	Smin /= log_p;
	
	return Smin;
}


static double get_Smin_2seq(float dGC, int32_t len1, int32_t len2, double pvalue){

	double Smin=0.0;
	
	double log_log_minusF = 0.0,
	       log_gamma      = 0.0,
	       log_p          = 0.0;


	log_log_minusF = log( -log( 1.00 - pvalue ) );
	log_gamma = get_log_gamma( dGC, (len1+len2)/2 );    /* average dGC and average len */
	log_p = get_log_p( dGC, (len1+len2)/2 );	
	
	
	Smin = log_log_minusF - log_gamma - log( (double)len1 * len2 );  /* when two sequences are considered,
	                                                                    the whole matrix [l1 l2] can potentially have repeats
									 */
	Smin /= log_p;
	
	return Smin;
}







double compute_smin(char *sequence, double pvalue, int32_t seq_len, int32_t nX){

	double Smin=0.0;
	float dGC=0.0;


	dGC = get_dGC(sequence, seq_len);
	Smin = get_Smin( dGC, seq_len-nX, pvalue);
	
	return Smin;
}



double compute_smin_2seq(char *sequence1, char *sequence2, int32_t size1,  int32_t size2,   double pvalue, int32_t nX, int32_t nX2){

	double Smin=0.0;
	float dGC=0.0;
	
	dGC = get_dGC_2seq(sequence1, size1, sequence2, size2 );
	
	Smin = get_Smin_2seq( dGC, size1-nX, size2-nX2, pvalue);
	
	return Smin;
}




void UpdateRepeats(Repeats *AllRepeats, double Score_min){

	int32_t i;

	if( AllRepeats->nDirRep )
		for(i=0;i<AllRepeats->nDirRep;i++)
			if(AllRepeats->DirRep[i].pos1 != -1 && AllRepeats->DirRep[i].score<Score_min){
				AllRepeats->DirRep[i].pos1 = -1;
				AllRepeats->DirRep[i].pos2 = -1;
				AllRepeats->DirRep[i].len1 = -1;
				AllRepeats->DirRep[i].len2 = -1;
				AllRepeats->DirRep[i].score = -1;
				AllRepeats->nDirBadRep++;
			}

	if( AllRepeats->nInvRep )
		for(i=0;i<AllRepeats->nInvRep;i++)
			if(AllRepeats->InvRep[i].pos1 != -1 && AllRepeats->InvRep[i].score<Score_min){
				AllRepeats->InvRep[i].pos1 = -1;
				AllRepeats->InvRep[i].pos2 = -1;
				AllRepeats->InvRep[i].len1 = -1;
				AllRepeats->InvRep[i].len2 = -1;
				AllRepeats->InvRep[i].score = -1;
				AllRepeats->nInvBadRep++;
			}
}


