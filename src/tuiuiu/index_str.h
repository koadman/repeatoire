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

#ifndef __index_str
#define __index_str
#include "tuiglobals.h"
#define pot4(x) (1 << (2 * (x))) /* Returns 4^x */

typedef struct INDEX {
   int *ind; 
   int *tuilist;
} index_str;  

//char ACTGnumber[256];

void constructor(int N, int k, char seq[], index_str *k_factor_ind);

#endif // __index_str
