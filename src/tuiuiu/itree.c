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
#include "itree.h"
#define TRUE 1
#define FALSE 0
//extern int z;
//extern float w_bin_size;
#define BEGIN(i, bin, z, N) (i + N - 1 - ((bin + 1) << z))
#define END(i, bin, z, w, N) (i + w + N - 1 - (bin << z))


itree *NewTree(void)
{
   return NULL;
}

void AddTree(int bin, itree **l)
{
   itree *z = (itree *)malloc(sizeof(itree));
   
   z->bin = bin;
   z->parent = NULL; 
   z->left = NULL; z->right = NULL;
   /* printf("addtree\n"); */  
   itree_insert(l, z);
   /* printf("addtree!!!\n"); */  
} 

void RemoveTree (int bin, itree **l)
{
   itree *x = *l;
   
   // printf("removetree\n");
   while (x->bin != bin)
     {
       if (bin < x->bin) 
	 x = x->left;
       else
	 x = x->right;
     }

   itree_delete(l, x);
   //printf("removetree!!!\n");
}     
   
void FreeTree (itree *l)
{
   if (l != NULL)
   { 
      FreeTree(l->left);
      FreeTree(l->right);
      free(l);
   }
} 

void itree_insert(itree** root, itree* z)
{
   itree* x;
   itree* y;
   
   y = NULL;
   x = *root;
   while (x) 
   {
      y = x;
      if (z->bin < x->bin) x = x->left;
      else  x = x->right;
    
   }
   z->parent = y;
   if (!y) *root = z;
   else 
   {
     if (z->bin < y->bin) y->left = z;
     else y->right = z;
   }
}

itree *itree_sucessor(itree *x)
{
   x =  x->right; 
   while (x->left != NULL)
      x = x->left;
   return x;
} 
      

void itree_delete(itree** root, itree* z)
{
   itree *x, *y; 
  
   if (z->left == NULL || z->right == NULL) y = z;
   else y = itree_sucessor(z);
   
   if (y->left != NULL) x = y->left;
   else x = y->right;
   
   if (x != NULL) x->parent = y->parent;
   
   if (y->parent == NULL) *root = x;
   else
   {
      if (y == y->parent->left) y->parent->left = x;
      else y->parent->right = x;
   }
   if (y != z)
      z->bin = y->bin;
   free(y);
}




