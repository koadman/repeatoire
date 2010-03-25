/* ***********************************************************************                                                                                                                                        
   Copyright (C) 2010  Maria Federico, Pierre Peterlongo, Gustavo Sacomoto                                                                                                                All Rights Reserved.                                                                                                                                                                                           
   This file is part of TUIUIU.                                                                                                                                                                                   
   TUIUIU is free software: you can redistribute it and/or modify                                                                                                                         it under the terms of the GNU General Public License as published by                                                                                                                   the Free Software Foundation, either version 3 of the License, or                                                                                                                      (at your option) any later version.                                                                                                                                                                            
   TUIUIU is distributed in the hope that it will be useful,                                                                                                                              but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                                                                         MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                                                                                          GNU General Public License for more details.                                                                                                                                                                                                                                                                                                                                  You should have received a copy of the GNU General Public License                                                                                                                      along with TUIUIU.  If not, see <http://www.gnu.org/licenses/>.                                                                                                                                                
   *********************************************************************** */                                                                                                                                                                                                                              

#ifndef mono_and_multi_commons_h
#define mono_and_multi_commons_h

#include "libMems/Tuiuiu/list.h"
#include "libMems/Tuiuiu/itree.h"
#include "libMems/Tuiuiu/emptyBlock.h"
#include <stdlib.h>
#include <stdio.h>

#define pot4(x) (1 << (2 * (x))) /* Returns 4^x */

#define MINP  8

int * tui_possiblePositions;	/* TODD value, storing the possible positions of excelent parralelograms */
int tui_posinpossiblepositions;	/* TODD optio: position in the  possiblePositions array*/
//int * tui_possiblePositions;	/* TODD value, storing the possible positions of excelent parralelograms */
//int tui_posinpossiblepositions;	/* TODD optio: position in the  possiblePositions array*/

int tui_lastbin; /* in order to count only projections inside a parallelogram */
int tui_erro;

typedef struct INDEX {
  int *ind;
  int *tlist;
} index_str;

void constructor(int N, int k, char seq[], index_str *k_factor_ind);
inline double cpuTime();
int readsequence(FILE *file, char**seq, char **name);
int readsequenceAndReverse(FILE *file, char**seq, char **name);
int isThereASequence (FILE *file);
inline void CountBins(int d, int e, int p, int bins[], itree **l);
inline void UncountBins(int d, int e, int p, int bins[], itree **l);
void ReportPgram(int begin, int end, tlist **result, empty_block goodWindows[], int curr_parall_min, int curr_parall_max, int curr_pass_num);
int Nelems(int *begin, int *end, itree *l, int i, int w, int bin_size, int N, empty_block goodWindows[], int curr_pass_num, int curr_parall_min, int curr_parall_max);

#endif // mono_and_multi_commons_h
