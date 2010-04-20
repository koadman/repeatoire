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
#ifndef TUI_UTIL_H
#define TUI_UTIL_H

#include "itree.h"
#include "index_str.h"
#include "emptyBlock.h"
#include "tuilist.h"
#include "tuiglobals.h"

//int z;

inline void CountBins(int d, int e, int p, int bins[], itree **l);

inline void UncountBins(int d, int e, int p, int bins[], itree **l);

//void ReportPgram(int begin, int end, tuilist **result, int goodWindows[]);

void ReportPgram(int begin, int end, tuilist **result, 
		 empty_block goodWindows[], int curr_parall_min, int curr_parall_max, int curr_pass_num);

tuilist *Filter(int N, int k, int p, int e, int r, int bin_size, 
             const int w, index_str *k_factor_ind, char seq[], 
	     empty_block goodWindows[], int curr_pass_num);

#endif
