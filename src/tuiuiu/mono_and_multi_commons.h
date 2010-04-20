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
#ifndef MONO_AND_MULTI_COMMONS_H
#define MONO_AND_MULTI_COMMONS_H
#include "tuilist.h"
#include "emptyBlock.h"
#include "util.h"
#include "tuiglobals.h"
//#include "tuiglobals.c"
#include <stdlib.h>
#include <stdio.h>
#define MINP  8

//int lastbin; /* in order to count only projections inside a parallelogram */


inline double cpuTime();
int readsequence(FILE *file, char**seq, char **name);
int readsequenceAndReverse(FILE *file, char**seq, char **name);
int reverseSequence(char**seq, int sizeseq);
int isThereASequence (FILE *file);
inline void CountBins(int d, int e, int p, int bins[], itree **l);
inline void UncountBins(int d, int e, int p, int bins[], itree **l);
void ReportPgram(int begin, int end, tuilist **result, 
		 empty_block goodWindows[], int curr_parall_min, int curr_parall_max, int curr_pass_num);

#endif
