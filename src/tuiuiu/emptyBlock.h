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

#ifndef __empty_block
#define __empty_block

/*EMPTY BLOCK*/
// type of elements stored in the array goodWindows
// notEmpty=1 if a block contains at least one window kept by the filter
// pass_num is the number of the pass in which the field notEmpty is updated
typedef struct EMPTYB {
   int notEmpty; 
   int pass_num;
} empty_block;  
/**/

#endif
