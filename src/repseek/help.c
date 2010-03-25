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
/******
        file     : help.c
        function : Help for all repseek series
                                                 
        created  : 02 Oct 2002
        modif    : Oct 2003
                   Dec 2004 Adds for option -D : managing of mask for seed detection <EC>
                   Dec 2005 switch to a one mode version

        author   : amikezor
*****/

#include "repseek_types.h"

#include <stdio.h>
#include <stdlib.h>

static void help_contact(){
	fprintf(stderr,"   Contact\n"
	               "     if you need help, want to report a bug or suggest an option,\n"
		       "     please contact G. Achaz by email: achaz(at)abi.snv.jussieu.fr\n"
	               "\n"
		       );
		       
}



static void help_information(){
	fprintf(stderr,"   Information\n"
	               "\t-v           : Print 'v'ersion and exit\n"
	               "\t-h           : Print a short 'h'elp and exit\n"
	               "\t-H           : Print a larger 'H'elp on a specific mode and exit\n"
	               "\n");
}

static void help_corevalues(){
	fprintf(stderr, "   Core values (select either at least one of the four -l/-p  -L/-P option)\n"
			"      seed minimum length: lmin\n"
			"\t-l lmin      : set a minimum length for the seeds detection (no pvalue).\n"
			"\t-p prob      : set a p-value that gives a minimum seed length (Karlin & Ost equation)\n"
			"      repeat minimum score: smin\n"
			"\t-L smin      : set a minimum score for extended repeat (def is 0: no filter).\n"
			"\t-P prob      : set a p-value that gives a minimum repeat score (Waterman & Vingron regression)\n"
			"\n"
	);
}

static void help_options( ){
	
	fprintf(stderr,"   Optionnal values\n"

	               "\t-r file      : 'r'edirect output into file (otherwise stdout is used).\n"
	               "\t-T           : print out R-'T'able\n"
		       "\t               the R-Table (or R-Table2) shows all non-unique position and its degree of redundancy\n"

	               "\t-S           : output 'S'eeds and exit (do not perform extention).\n"
	               "\t-s seed_file : use seeds given in seed_file instead of detecting them\n"
		       "\t               File format is 'd|i begin end len [optionnal other fields]'.\n"
                       "\t-D mask_file : mask sequence regions during seed detection only (cannot be used with -s seed_file).\n"
                       "\t               File format is 'begin end [seq#]' (seq# is 1 or 2; nothing is treated as 1).\n"
                       "\t-R 0.##      : merge 'R'epeats when they share 0.## of their both copies (default = 0.90).\n"
		       "\t               When set to 1.0, merge repeats having exactly the same positions\n"

	               "\t-d           : detect only 'd'irect repeats.\n"
	               "\t-i           : detect only 'i'nverted repeats.\n"
	               "\t               (-d and -i are mutually exclusive).\n"
	               "\t-B           : if some seeds overlaps (i.e. in low-complexity sequence), just keep the 'B'iggest one.\n"
	               "\t-m #.#       : keep only seeds occurring at least the specified minimum number of times (min is a float).\n"
	               "\t-M #.#       : keep only seeds occurring less than the specified maximum number of times (Max is a float).\n"
	               "\t-O 0|1       : for 1 sequence, direct repeat can have their 2 copies 'O'verlapping (0: no, 1: yes -default-)\n"
	               "\t-c           : set the chromosome as 'c'ircular -default is linear- (unused for 2 sequences)\n"
	               "\t-o #.#       : set gap_open penalty -default 4.0- (cannot be change when -P is used)\n"
	               "\t-e #.#       : set gap_ext penalty -default 1.0- (cannot be change when -P is used)\n"
	               "\t-X #.#       : set 'X'g, the 'exploration value' -default 30.0-\n"
	               "\n"
		       );


}

static void help_input( void ){	

	fprintf(stderr,"   Input\n"
	               "     Input file(s) should be in a fasta format.\n"
	               "     All non-ACGT characters are considered as N, except for X in which no repeats are detected (mask)\n"
	               "\n"
		       );
}

static void help_output( void ){	

	fprintf(stderr,"   Output\n"
	               "\tRepeats are displayed in 12 COLUMNS\n"
                       "\t 01 - type of the repeat (Tandem|Close|Distant|Overlap|Palindrome|Interseq).(dir|inv)\n"
	               "\t        'Tandem'     : repeat is a perfect tandem (spacer=0)\n"
	               "\t        'Close'      : spacer is smaller than 1 kb\n"
	               "\t        'Distant'    : spacer is larger  than 1 kb\n"
	               "\t        'Overlap'    : repeat is composed of two overlapping copies (spacer<0)\n"
	               "\t        'Palindrome' : repeat is a perfect palindrome (spacer=0)\n"
	               "\t        'Interseq'   : repeat is shared by two sequences\n"
	               "\t 02 - position of the first copy\n"
	               "\t 03 - position of the second copy\n"
	               "\t 04 - length of the first copy\n"
	               "\t 05 - length of the second copy\n"
	               "\t 06 - spacer of the repeats (corrected when the sequence is circular)\n"
	               "\t 07 - characteristics of the seed that gave that repeat (pos1-pos2-len-meanR)\n"
	               "\t 08 - %%identity between the two copies (matches / alignment_length)\n"
	               "\t 09 - score of the alignement\n"
	               "\t 10 - mean-R of the repeat\n"
	               "\t 11 - mode-R of the repeat\n"
	               "\t 12 - fraction of its mode-R\n"
	               "\n"
				);

	fprintf(stderr,"\tSeeds are displayed in 5 COLUMNS\n"
	               "\t 01 - Type of the repeat (Seed).(dir|inv)\n"
	               "\t 02 - position of the first copie\n"
	               "\t 03 - position of the second copie\n"
	               "\t 04 - length of the seed\n"
	               "\t 05 - mean-R\n\n"
			);

}





void printversion_repseek(){
	fprintf(stderr,"repseek version %s (%s)\n", REPSEEK_VERSION,  REPSEEK_DATE);
}


void printusage_repseek(char **argv){

	fprintf(stderr,"General usage is '%s [-v|-h|-H] [ opts ] { core_value } file.fst [file2.fst]'\n\n",argv[0]);
	help_corevalues();
	help_information();
}

void printhelp_repseek(char **argv){

	fprintf(stderr,"General usage is '%s [-v|-h|-H] [ opts ] { core_value } file.fst [file2.fst]'\n\n",argv[0]);
	
	help_corevalues();
	help_options( );
	help_information();

	help_input();
	help_contact();
}

void printmorehelp_repseek(char **argv){

	fprintf(stderr,"/*\n"
	               "\tRepseek looks for all repeats couples within (one fasta file) or between sequences (two files)\n"
		       "\tThe method is based on a two steps method (seed detection followed by their extension into repeats)\n"
		       "\tThe statistical significance of repeats can be (i) not evaluated (-l opt),\n"
		       "\t(ii) evaluated for seeds (-p opt) or (iii) evaluated for repeats (-P opt). For further details, \n"
		       "\tplease refer to the manuscript (Achaz et al., in prep)\n"
		       "*/\n"
		       );
	
	fprintf(stderr,"General usage is '%s [-v|-h|-H] [ opts ] { core_value } file.fst [file2.fst]'\n\n",argv[0]);
	
	help_corevalues();
	help_options( );
	help_information();

	help_input();
	help_output();

	help_contact();
}

