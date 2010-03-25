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

#ifndef __tuilist
#define __tuilist

typedef struct TUILIST {
   int begin;
   int end;
   struct TUILIST *prox;
} tuilist;

tuilist *NewList(void);

tuilist *Add(tuilist *l, int begin, int end);

void FreeList(tuilist *l);

tuilist *InvertList(tuilist *l);
tuilist * ManageResultsForReverse(tuilist *l, int size);


#endif
