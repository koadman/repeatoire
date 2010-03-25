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
/*
 *  KMRK_mask.h
 *  repseek
 *
 *  Created by Eric Coissac on 04/12/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef KMRK_MASK_H
#define KMRK_MASK_H

#include "repseek_types.h"

typedef struct {
    int32_t begin;
    int32_t end;
} masked_area_t;

typedef struct {
    int32_t reserved;
    int32_t count;
    
    masked_area_t area[1];
} masked_area_list_t;

typedef struct {
    int32_t seqcount;
    int32_t total;
    
    masked_area_list_t *sequence[1];
} masked_area_table_t;

masked_area_table_t *KMRK_ReadMaskArea(char *areafile,int32_t seqcount,int32_t complement);
char KMRK_isMasked(masked_area_table_t *mask,int32_t seq, int32_t position);

#endif
