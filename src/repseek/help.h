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

#ifndef _help_h
#define _help_h


/******
        file     : help.c
        function : Help for all repseek series
                                                 
        created  : 02 Oct 2002
        modif    : Oct 2003

        author   : amikezor
*****/


#include <stdio.h>
#include <stdlib.h>

/********* *********

	All you need for RepSeek
	
********* *********/

void printversion_repseek();
void printusage_repseek(char **argv);
void printhelp_repseek(char **argv);
void printmorehelp_repseek(char **argv);

#endif

