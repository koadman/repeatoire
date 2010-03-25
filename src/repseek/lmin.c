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
	file     : lmin.c
	function : get lmin from a pvalue
	created  : Oct 11 2002
				
	note   :   almost all functions comes from lari_stat.c
              Apr. 97 <AV> first draft		
              April 98 <ed> K formula added	
				  Oct 2002 <amikezor> - subset of lari_stat - add LminRiskSeq()
*****/


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "repseek_types.h"
#include "memory.h"

#define NALPHA 4

typedef struct {
    int32_t nalpha;
    int32_t total_symb, val_symb, bad_symb;
    int32_t count[NALPHA+1];
    float freq[NALPHA+1];
    char  alpha[NALPHA+2];
} FreqTable;

#undef NALPHA



/***
	Few function coming from numrec as it
***/

static void nrerror(char *error_text)

{
	fprintf(stderr,"> Numerical Recipes run-time error...\n");
	fprintf(stderr,"> %s\n", error_text);
	exit(1);
}

static float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

#define ITMAX 100
#define EPS 3.0e-7

static void gser(float *gamser, float a, float x, float *gln)
{
	int n;
	float sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}

#define FPMIN 1.0e-30

static void gcf(float *gammcf, float a, float x, float *gln)
{
	int i;
	float an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

#undef ITMAX
#undef EPS
#undef FPMIN

static float gammq(float a, float x)
{
	float gamser=0.0,
	      gammcf=0.0,
	      gln=0.0;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}




/* ----------------------------------------------- */
/* compute symbol frequencies			   */
/* from seq and alphabet			   */
/* ----------------------------------------------- */

int InitFreqTable(char *seq, char *alpha, FreqTable *table)
{
    int   i;
    float norm;
    char  *s, *found;

    strcpy(table->alpha, alpha);
    
    table->nalpha = strlen(alpha);
    
    for (i = 0 ; i <= table->nalpha ; i++)
		table->count[i] = 0;
	
    table->total_symb = 0;       	
    for (s = seq; *s; s++) {
		table->total_symb++;      	
		if ( (found=strchr(alpha, *s) ))
	    	table->count[(int)(found - alpha)]++;
		else
	    	table->count[table->nalpha]++; 
    }
    
    table->bad_symb = table->count[table->nalpha];
    table->val_symb = table->total_symb - table->bad_symb;
    
    norm = MAX(1., table->val_symb);

    for (i = 0 ; i < table->nalpha ; i++)
		table->freq[i] = (float) table->count[i] / norm;

    return table->val_symb;	
}

/* ----------------------------------------------- */
/* compute x^n					   */
/* ----------------------------------------------- */

static double sPower(double x, int n)
{
    double pow = 1.;
    
    while (n-- > 0)
	pow *= x;
    
    return pow;    
}

/* ----------------------------------------------- */
/* compute lambdaR				   */
/* ----------------------------------------------- */

static double LambdaR(FreqTable *table, int r)
{
    int	    i;
    double  lamb = 0.;
    
    for (i = 0 ; i < table->nalpha ; i++)
		lamb += sPower((double)table->freq[i], r);    

    return lamb;
}

/* ----------------------------------------------- */
/* log of binomial coeff 			   */
/* ----------------------------------------------- */

static double sLogBinoCof(int n, int p)
{
    double lb = 0.;
    
    while (p >= 1)
		lb += log((double)n-- / (double)p--);
    
    return lb;    
}

/* ----------------------------------------------- */
/* compute a(x,y)				   */
/* ----------------------------------------------- */

static double sAlpha(double x, double y)
{
    double x1, x2;
    
    x1 = (log(1. - x) + log(y)) / log(x);
    x2 = 0.57772 / log(x);
    
    return -(x1 + x2);
}


/* ----------------------------------------------- */
/* compute E(L(N,r))                               */
/* ----------------------------------------------- */

static double ExpectL(FreqTable *table, int N, int r)
{
    double x1, x2, lambr;
    
    lambr = LambdaR(table, r);
    x1    = sLogBinoCof(N, r) / -log(lambr);
    x2    = sAlpha(lambr, lambr);
    
    return x1 + x2 + 0.5;
}

/* ----------------------------------------------- */
/* compute Var(N,r)				   */
/* ----------------------------------------------- */

double VarianceL(FreqTable *table, int r)
{
    double x;
    
    x = 1. / log(LambdaR(table, r));
    
    return 1.645 * x * x;
}
/* ----------------------------------------------- */
/* compute recursive rlogN			   */
/* ----------------------------------------------- */

static void finN(indx, r, data, res)
    int *indx,  r, *data;
    double *res;
{
    int i;
    double aux;

    for (i = 0, aux=1. ; i < r ; i++)
		aux *=  data[indx[i]];

	*res += aux;
	
		
}

static void recursN(indx, s, r, rcurr, ipos, data, res) 
	int	*indx, *data;
	double *res;
	int	s, r, rcurr, ipos;
{
	int   i;
	
	if (rcurr == r) 	/* end of recursion 	*/
	    finN(indx, r, data, res);
	else {			/* body of recursion	*/
	
	   for (i = ipos ; i < s ; i++) {
	       indx[rcurr] = i;
	       recursN(indx, s, r, rcurr+1, i+1, data, res);
	   }
	}
}

/* ----------------------------------------------- */
/* compute recursive Rlambda			   */
/* ----------------------------------------------- */

static void finR(indx, r, data, res)
    int 		*indx,  r;
    double 		*res;
    FreqTable 	*data;
{
    int i, j;
    double aux, aux2;

	for (j=0, aux2=0; j<4; j++){
		for (i = 0, aux=1. ; i < r ; i++)
			aux *= (double) data[indx[i]].freq[j];
		aux2 += aux;
		}
		
	if (aux2>*res)
		*res = aux2;

}

static void recursR(indx, s, r, rcurr, ipos, data, res) 
	int	*indx, s, r, rcurr, ipos;
	double *res;
	FreqTable *data;
{
	int   i;
	
	if (rcurr == r) 	/* end of recursion 	*/
	    finR(indx, r, data, res);
	else {			/* body of recursion	*/
	
	   for (i = ipos ; i < s ; i++) {
	       indx[rcurr] = i;
	       recursR(indx, s, r, rcurr+1, i+1, data, res);
	   }
	}
}

/* ----------------------------------------------- */
/* binomial coeff 				    */
/* ----------------------------------------------- */

static double sBinoCof(int n, int p)
{
    double b = 1.;
    
    while (p >= 1)
	b *= ((double)n-- / (double)p--);
    
    return b;    
}

/* ----------------------------------------------- */
/* compute E(K(N,r))                               */
/* ----------------------------------------------- */



static double ExpectK(FreqTable *table, int r, int s)
{
    double  x2, x3, lambr;
    int	    *indx, *data, i, nelt;
    double x1;
    
    if (s!=r)
    	nelt = ceil(sBinoCof(s, r));
    else
    	nelt = 1;
    
    indx = (int32_t *)MyCalloc( nelt, sizeof(int32_t), "ExpectK out of memory" ) ;
    data = (int32_t *)MyCalloc( s, sizeof(int32_t), "ExpectK out of memory" ) ;
	 
    x1=0.;
    
    for (i=0; i<s; i++)
    	data[i] = table[i].total_symb;
    
    recursN(indx, s, r, 0, 0, data, &x1);

	lambr=0;
	
    recursR(indx, s, r, 0, 0, table, &lambr);
    
    x2 = log(x1) / -log(lambr);
    x3 = sAlpha(lambr, lambr);
    
    MyFree(indx, nelt* sizeof(int32_t));
    MyFree(data, s* sizeof(int32_t));
    
    return x2 + x3 + 0.5;
}

/* ----------------------------------------------- */
/* compute VarK (N,r)				   */
/* ----------------------------------------------- */

static double VarianceK(FreqTable *table, int r)
{
    double x;
    
    x = 1. / log(LambdaR(table, r));
    
    return 1.645 * x * x;
}

/* ----------------------------------------------- */
/* compute Erf(x)				   */
/* ----------------------------------------------- */
static double sErf(double x)
{
    double gamp;

    gamp = 1.0 - gammq(0.5, x * x);
    return ((x < 0.) ? -gamp : gamp);
}

/* ----------------------------------------------- */
/* compute cumulative of Normal			   */
/* ----------------------------------------------- */
#define SQRT2 1.4142135624
static double sCumNormal(double x, double mu, double sigma)
{
    double xx;
    
    xx = (x - mu) / sigma / SQRT2;

    return (1.0 + sErf(xx)) / 2.;    
}

#undef SQRT2

/* ----------------------------------------------- */
/* compute Lmin(N,r) at error level p		   */
/* ----------------------------------------------- */

#define TOLERANCE 1e-6

static double LminRisk(FreqTable *table, int N, int r, float p)
{
    double pp, mu, sigma, xmin, xmax, xcurr, pcurr;

    if (p <= 0.)
		return 1.;
    if (p >= 1.)
		return (double) N;

    pp    = 1.0 - (double) p;
    mu    = ExpectL(table, N, r);
    sigma = sqrt(VarianceL(table, r));

    xmin  = 1;
    xmax  = N;
    xcurr = (xmin + xmax) / 2. ;
    pcurr = sCumNormal(xcurr, mu, sigma);
    
    /* dichotomic search */
    
    while (fabs(pp - pcurr) > TOLERANCE) {
	if (pp > pcurr)
	   xmin = xcurr;
	else
	   xmax = xcurr;
	xcurr = (xmin + xmax) / 2. ;
	pcurr = sCumNormal(xcurr, mu, sigma);
    }

    return xcurr;        
}

/* ----------------------------------------------- */
/* compute K Lmin(N,r) at error level p		   */
/* ----------------------------------------------- */

static double KLminRisk(FreqTable *table, double rlogn, int r, int s, 
	float p, int len)
{
    double pp, mu, sigma, xmin, xmax, xcurr, pcurr;

    if (p <= 0.)
		return 1.;
    if (p >= 1.)
		return rlogn;

    pp    = 1.0 - (double) p;
    mu    = ExpectK(table, r, s);
    sigma = sqrt(VarianceK(table, r));

    xmin  = 1;
    xmax  = len;
    xcurr = (xmin + xmax) / 2. ;
    pcurr = sCumNormal(xcurr, mu, sigma);
    
    /* dichotomic search */
    
    while (fabs(pp - pcurr) > TOLERANCE) {
		if (pp > pcurr)
	   		xmin = xcurr;
		else
	   		xmax = xcurr;
		xcurr = (xmin + xmax) / 2. ;
		pcurr = sCumNormal(xcurr, mu, sigma);
    }

    return xcurr;        
}

#undef TOLERANCE 

/***
	Take a sequencem the number of occurence and the prob and
	return an associated lmin
***/

#define DNA_ALPHABET "ACGT"

static double LminRiskSeq(char *sequence, int r, float p)
{

	int nSymbols=0;
	char alphabet[32]={0,};
	FreqTable ftable;

	double Lmin=0.0;

	(void) strcpy(alphabet, DNA_ALPHABET);
	
	(void) InitFreqTable(sequence, alphabet, &ftable);
	
	nSymbols = ftable.val_symb;
	
	Lmin = LminRisk(& ftable, nSymbols, r, p);
	
	return Lmin;
}



int16_t set_lmin(float pval, char *sequence){

	float flmin;
	int16_t lmin;

	flmin = LminRiskSeq(sequence, 2, pval);
	
	if(flmin - (float)(int16_t)flmin < 0.5)
		lmin = (int16_t)flmin;
	else
		lmin = (int16_t)flmin + 1;

	return lmin;
}

/***
	Take a sequence the number of occurence and the prob and
	return an associated lmin
***/

#define DNA_ALPHABET "ACGT"

static double LminRisk2Seq(char *Sequence1, 
			   char *Sequence2, 
			   int32_t SeqSize1, 
			   int32_t SeqSize2, 
			   int r, float p)
{

	char alphabet[32]={0,};  /* an theoretical bigger alphabet than ACGT */
	FreqTable ftable[2];     /* a freq for each sequence */

	int nseq=2;              /* number of sequences */
	
	double Lmin=0.0;         /* the value of interest */

	(void) strcpy(alphabet, DNA_ALPHABET);
	
	(void) InitFreqTable(Sequence1, alphabet, &ftable[0]);           /* init freq of the first seq */
	(void) InitFreqTable(Sequence2, alphabet, &ftable[1]);           /* init freq of the second seq */

	Lmin = KLminRisk(ftable, (double) r*log(SeqSize1+SeqSize2/nseq), r, nseq, p, SeqSize1+SeqSize2);
	
	return Lmin;
}

#undef DNA_ALPHABET


int16_t set_lmin_2seqs(float pval, 
		       char *Sequence1, 
		       char *Sequence2, 
		       int32_t SeqSize1, 
		       int32_t SeqSize2)
{

	double flmin;
	int16_t lmin;
	int8_t  r=2;

	flmin = LminRisk2Seq(Sequence1, Sequence2, SeqSize1, SeqSize2, r, pval);
		
	if(flmin - (double)(int16_t)flmin < 0.5)
		lmin = (int16_t)flmin;
	else
		lmin = (int16_t)flmin + 1;

	return lmin;
}





#undef DNA_ALPHABET

