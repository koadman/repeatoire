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
 *  KMRK_mask.c
 *  repseek
 *
 *  Created by Eric Coissac on 04/12/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "KMRK_mask.h"
#include <stdio.h>
#include <stdlib.h>
#include "memory.h"

#define MASKED_AREA_TABLE_SIZE(seqcount) (sizeof(masked_area_table_t) + (sizeof(masked_area_list_t*) * ((seqcount)-1)))
#define MASKED_AREA_LIST_SIZE(areacount) (sizeof(masked_area_list_t) + (sizeof(masked_area_t) * ((areacount)-1)))

#define AREA_COUNT_INIT (1000)

static masked_area_table_t *new_masked_area_table(int32_t seqcount, int32_t areacount);
static masked_area_list_t *new_masked_area_list(int32_t areacount);
static masked_area_list_t *realloc_masked_area_list(masked_area_list_t *list,int32_t areacount);
static int32_t push_area(masked_area_table_t* table,int32_t sequence,int32_t begin,int32_t end);
static void sort_area_table(masked_area_table_t* table);
static int32_t compare_area(const masked_area_t* v1,const masked_area_t* v2);
static int32_t search_area(const masked_area_t* v1,const masked_area_t* v2);

static masked_area_list_t *strip_area_list(masked_area_list_t* list);
static void strip_area_table(masked_area_table_t* table);


static masked_area_list_t *new_masked_area_list(int32_t areacount)
{
    masked_area_list_t *list;
    
    list = MyCalloc(1, MASKED_AREA_LIST_SIZE(areacount), "Not enougth memory for mask table");
    list->reserved=areacount;
    
    return list;
}

static masked_area_list_t *realloc_masked_area_list(masked_area_list_t *list,int32_t areacount)
{    
    list = MyRealloc(list,
                     MASKED_AREA_LIST_SIZE(areacount),
                     MASKED_AREA_LIST_SIZE(list->reserved),
                     "Not enougth memory for mask table");

    list->reserved=areacount;
    
    return list;
}

static masked_area_table_t *new_masked_area_table(int32_t seqcount, int32_t areacount)
{
    masked_area_table_t *table;
    int32_t i;
    
    table = MyCalloc(1, MASKED_AREA_TABLE_SIZE(seqcount),"Not enougth memory for mask table");
    table->seqcount=seqcount;
    
    for (i=0;i<seqcount;i++)
        table->sequence[i]=new_masked_area_list(areacount);
    
    return table;
}

static int32_t push_area(masked_area_table_t* table,int32_t sequence,int32_t begin,int32_t end)
{
    masked_area_list_t * list;
    
    if (sequence >= table->seqcount)
        return -1;
    
    list = table->sequence[sequence];

    if (list->reserved == list->count)
    {
        list = realloc_masked_area_list(list,list->reserved*2);
        table->sequence[sequence]=list;
    }
    
    list->area[list->count].begin=begin;
    list->area[list->count].end=end;
    
    list->count++;
    table->total++;

    return table->total;
}

static int32_t compare_area(const masked_area_t* v1,const masked_area_t* v2)
{
    return v1->begin - v2->begin;
}

static void sort_area_table(masked_area_table_t* table)
{
    int32_t i;
    
    for (i=0; i<table->seqcount;i++)
    {
        qsort(table->sequence[i]->area,
              table->sequence[i]->count,
              sizeof(masked_area_t),
              (int (*)(const void *, const void *))compare_area);
    }
}

static masked_area_list_t *strip_area_list(masked_area_list_t* list)
{
    int32_t i;
    int32_t j;
    int32_t count;
    int32_t newcount;
    
    count = list->count;
    newcount=count;
    
    for (i=1;i<count;i++)
    {
       /* fprintf(stderr,"\n%d->%d    %d->%d ==>",list->area[i-1].begin,list->area[i-1].end,list->area[i].begin,list->area[i].end); */
        if ((list->area[i].begin-1) <= list->area[i-1].end)
        {
            /* fprintf(stderr," joined"); */
            list->area[i].begin=list->area[i-1].begin;
            list->area[i-1].begin=-1;
            newcount--;
        }
    }
    list->count=newcount;

    for (i=0,j=0;i<count;i++)
    {
        if (list->area[i].begin>=0)
        {
            if (i!=j)
                list->area[j]=list->area[i];
            j++;
        }
    }
    
    return realloc_masked_area_list(list,newcount);

}

static void strip_area_table(masked_area_table_t* table)
{
    int32_t seq;
    int32_t oldcount;
    masked_area_list_t* list;
    
    sort_area_table(table);
    
    for (seq=0; seq < table->seqcount;seq++)
    {
        list = table->sequence[seq];
        oldcount = list->count;
        table->sequence[seq]=strip_area_list(list);
        table->total-=oldcount - table->sequence[seq]->count;
    }
}


static int32_t search_area(const masked_area_t* v1,const masked_area_t* v2)
{
    int32_t pos;
    
    pos = v1->begin;
    
    if (pos < v2->begin)
            return -1;
    
    if (pos > v2->end)
            return 1;
    
    return 0;
}


masked_area_table_t *KMRK_ReadMaskArea(char *areafile, int32_t seqcount, int32_t complement)
{
    FILE* area;
    char buffer[1000];
    char *ok;
    int32_t begin;
    int32_t end;
    int32_t sequence;
    int32_t column;
    int32_t linecount;
    masked_area_table_t *table;
    
    if (complement > 0)
        seqcount++;
    else
        complement=0;
    
    area = fopen(areafile,"r");
    linecount=0;
    table = new_masked_area_table(seqcount,AREA_COUNT_INIT);
    
    do {
        linecount++;
        ok = fgets(buffer,999,area);
        if (ok)
        {
            column = sscanf(buffer,"%d %d %d",&begin,&end,&sequence);
            if (column > 1 && begin <= end)
            {
                begin--;
                end--;
                if (column==3)
                    sequence--;
                else
                    sequence=0;

              if (sequence>1 || sequence<0)
                fprintf(stderr,"ERROR while reading mask file, line %d (see help)\n", linecount), exit(0);
              
                if (sequence && complement)
                    sequence++;
                
                push_area(table,sequence,begin,end);
                
                if (!sequence && complement)
                    push_area(table,1,complement -1 - end,complement -1 -begin);
                
            }
            
            if (column==1)
                fprintf(stderr,"EROOR in mask file reading line %d\n",linecount), exit(0);
            
        }
    } while (ok);
    
    fprintf(stderr,"\nread %d masked areas from file\n",table->total);
    strip_area_table(table);
    fprintf(stderr,"strip to %d non overlaping areas\n",table->total);
    
    return table;   
}

char KMRK_isMasked(masked_area_table_t *mask,int32_t seq, int32_t position)
{
    masked_area_t input;
    int32_t result;
    masked_area_list_t *list;
    
    if (! mask || (seq >= mask->seqcount))
        return 0;
    
    list = mask->sequence[seq];
    input.begin = position;
    
    result = bsearch(&input,
                     list->area,
                     list->count,
                     sizeof(masked_area_t),
                     (int (*)(const void *, const void *))search_area) != NULL;

        return result;
}
