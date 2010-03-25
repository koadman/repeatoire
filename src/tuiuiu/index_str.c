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
#include "index_str.h"

void constructor(int N, int k, char seq[], index_str *k_factor_ind)
{
  //k_factor_ind->ind = (int *)malloc(sizeof(int) * (pot4(k) + 1));
  //k_factor_ind->tuilist = (int *)malloc(sizeof(int) * N);
  extern char ACTGnumber[];   
  k_factor_ind->ind = (int *)calloc((pot4(k) + 1), sizeof(int)); //MULTIPASS: clean structures
  k_factor_ind->tuilist = (int *)calloc(N, sizeof(int));

   ACTGnumber['N'] = ACTGnumber['n'] = 0;  /* it could be any value between 0 and 3 */   
   ACTGnumber['A'] = ACTGnumber['a'] = 0;
   ACTGnumber['C'] = ACTGnumber['c'] = 1;
   ACTGnumber['T'] = ACTGnumber['t'] = 2;
   ACTGnumber['G'] = ACTGnumber['g'] = 3;
   
   int i, pos = 0, nextNpos = 0, Ncounter = 0;

   //for (i = 0; i <= pot4(k); i++) //BUG: <= instead of < and remove last line
       //k_factor_ind->ind[i] = 0; 
    
   /* Contando o numero de repeticoes de cada k-fator*/
   for (i = 0; i < k; i++)
   {
      if (seq[i] == 'N')
      { 
	 Ncounter++;
         nextNpos = i + k;
      }
      pos = pos + ACTGnumber[(int)seq[i]] * pot4(i);
   }

   if (i - 1 >= nextNpos)
      k_factor_ind->ind[pos]++; 
       
   for (i = k; i < N; i++)
   {  
     if (seq[i] == 'N')
      {
         nextNpos = i + k;
         Ncounter++;
      }
      pos = (pos >> 2) + (pot4(k - 1) * ACTGnumber[(int)seq[i]]);
      if (i >= nextNpos) 
         k_factor_ind->ind[pos]++;
   }
      
   /* Acumulando e copiando*/
   for (i = 1; i <= pot4(k); i++) //BUG: <= instead of < and remove last line
      k_factor_ind->ind[i] = k_factor_ind->ind[i] + k_factor_ind->ind[i - 1];
   
   /* Guardando as posicoes iniciais de cada k-fator em tuilist */
   pos = 0; nextNpos = 0;
   for (i = 0; i < k; i++)
   {
      if (seq[i] == 'N')
         nextNpos =  i + k; 
      pos = pos + ACTGnumber[(int)seq[i]] * pot4(i);
   }

   if (i-1 >= nextNpos)
   {
      k_factor_ind->tuilist[k_factor_ind->ind[pos] - 1] = 0; 
      k_factor_ind->ind[pos]--;
   }
   
   for (i = k; i < N; i++)
   {
      if (seq[i] == 'N')
	 nextNpos = i+k;
      pos = (pos >> 2) + (pot4(k - 1) * ACTGnumber[(int)seq[i]]);
      if (i >= nextNpos)
      {
         k_factor_ind->tuilist[k_factor_ind->ind[pos] - 1] = i - k + 1;
         k_factor_ind->ind[pos]--;
      }
   }
   //printf("%d %d %d!!!\n", N, Ncounter, k);
   //k_factor_ind->ind[pot4(k)] = N - k + 1 - Ncounter;  //BUG: removed
   //printf("\n%d \n", N-k+1-Ncounter - k_factor_ind->ind[pot4(k) - 1]);
}			   

/* int main(void)
{
   char Seq[101], alpha[] = "ACTG";
   int i, Len = 100;
   index_str *k_factor_ind = (index_str *)malloc(sizeof(index_str));  
 
   for (i=0; i<Len ;i++)
     Seq[i] = alpha[ (int)(4.0*random()/(RAND_MAX+1.0)) ];
   Seq[i] = '\0';
   constructor(Len, 3, Seq, k_factor_ind);
   printf("%s\n\n", Seq);
   for (i = 0; i < Len - 3 + 2; i++)
      printf("%d %d\n", i, k_factor_ind->tuilist[i]);
   printf("\n");
   for (i = 0; i <= pot4(3); i++)
     printf("%d %d\n", i, k_factor_ind->ind[i]);     
   printf("\n");
   
   free(k_factor_ind->tuilist);
    
   free(k_factor_ind->ind); 
   return 0; 
} */
            
 
