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

#ifndef itree_h
#define itree_h
#include "tuiglobals.h"
typedef struct _itree {
  struct _itree* parent;
  struct _itree* left;
  struct _itree* right;
  int bin;    /* key */
} itree;

//float w_bin_size;

void itree_delete(itree** root, itree* z);
itree *itree_sucessor(itree *x);
void itree_insert(itree** root, itree* z);
itree *NewTree(void);
void AddTree(int bin, itree **l);
void RemoveTree (int bin, itree **l);
void FreeTree (itree *l);

#endif // __itree
