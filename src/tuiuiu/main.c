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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h> 
#include <getopt.h>
#include <time.h>
#include <string.h>

#include "mono_and_multi_commons.h"
#include "list.h"
#include "index_str.h"
#include "itree.h" 
#include "util.h" 
#include "display.h"

int Nbases = 0;

extern char method[];
extern char version[];
extern char publication[];

 
void print_usage(){
  fprintf(stderr, "Usage : tuiuiu [OPTIONS] input_file ouput_information_file output_data_file\n"
	  " OPTIONS:\n"
	  " \t -w value : length of the window (default 100)\n"
	  " \t -e value : number of errors (substitutions, insertions and deletions) between two windows (default 12)\n"
	  " \t -r value : number of window searched (default 8)\n"
	  " \t -k value : length of the k-mers (default 6 decreasing to a good threshold)\n"
	  " \t -c       : ask to check direct AND reverse complement strands\n"
	  " \t -M       : ask to perform a filtering multi pass\n"
          " \t -N {0 1 2} : \n"
          "            2: Write Ns instead of filtered nucleotides, (default)\n"
          "            1: Concatenate concerved fragments, separated with one 'N' \n"
	  "            0: Concatenate concerved fragments\n");
  exit(1);
}

int main(int argc, char **argv)
{
  //const double multi = 0.98;  //constant to control the iteration of multi passes

   list *result;
   index_str *k_factor_ind = (index_str *)malloc(sizeof(index_str));
   FILE *input, *output, *output2;  
   char *seq, *name;
   int N, bin_size = 0, p;
   int w=100, e=12, r=8, k=6;
   int temoin, pot2, flag = 0;
   double kept = 1.0;
   //double oldKept = 1.0;
   int writeN=2;
   int rev = 0;
   //multi pass variables
   int multipass = 0, passNum=1;
   // empty block variables
   int numParal=0, g=0;
   empty_block *goodWindows;
   empty_block eb;

   //   fprintf(stderr, "method = %s\n", method); /* DEB */

   printf("Method = %s\n", method);

   long double timeStart = cpuTime(), timeStop;
   while (1)
     {
       temoin = getopt (argc, argv, "+r:e:k:w:b:N:c:M");
       if (temoin == -1){
	 break;
       }
       switch (temoin)
	 {
	 case 'w':
	   w=atoi(optarg);
	   break;
	 case 'e':
	   e=atoi(optarg);
	   break;
	 case 'r':
	   r=atoi(optarg);
	   break;
	 case 'k':
	   k=atoi(optarg);
	   if(k<1) {
	     fprintf(stderr, "Error, k must be bigger than 0"); exit(1);
	   }
	   break;
	 case 'c': 
	   rev=1;
	   break;
	 case 'M':
	   multipass=1;
	   break;
	 case 'N':
	   writeN=atoi(optarg);
	   break;
	 case 'b':
	   bin_size=atoi(optarg);
	   z = 0; pot2=1;
	   while (bin_size > 1) {
	     bin_size >>= 1;
	     z++;
	     pot2 <<= 1;
	   }
	   bin_size=pot2;
	   break;
	   /*   case 'n': */
	   /*            flag = 1; */
	   /*            break;   */
	 case '?':
	   break;
	 default:
	   printf ("Unknown option '%c'\n", temoin); print_usage();
	 }
     }
   if (argc  - optind !=3)
     {
       print_usage();
       exit(1);
     }
   
   if (!bin_size) {
     bin_size=1; z = 0;
     while (bin_size <= e)
       {
	 bin_size <<= 1;
	 z++;
       }
     /* z--;
	bin_size >>= 1; */
   }

   p = (w + 1) - k * (e + 1);
   printf( "p (filtration condition) = %3d >= %3d\r\n", p, MINP*w/100);
   if(p<MINP*w/100){
     fprintf(stderr,"Error: with such parameters, Tuiuiu is useless. Please Restart Tuiuiu with more stringent parameters.\n");
     exit(0);
   }
   input = fopen(argv[optind], "rt");
   if(input == NULL) {fprintf(stderr, "Cannot read input sequence \"%s\", exit\n",argv[optind]); exit(1);}


   if(!rev){ 
     N = readsequence(input, &seq, &name);
     seq = (char *) realloc(seq, (N+1)*sizeof(char) );
     seq[N] = '\0';
   }
   else N = readsequenceAndReverse(input, &seq, &name);

   if(isThereASequence (input)) {
     fprintf(stderr,"Error: Tuiuiu mono cannot be applied to more than exactly one sequence. Please Restart Tuiuiu Mono with exactly one sequence\n");
     exit(0);
   }
   fclose(input);

   /*EMPTY BLOCK*/
   // goodWindows is an array of elements of type empty_block [mono_and_multi_commons.h]
   // of size the number of blocks
   //goodWindows[i] is 1 iff i-th block contains
   //at least one window kept by the filter (it is
   // updated in ReportPgram function [mono_and_multi_commons.c])
   // the creation and initialization of this array
   // is made here because the info in it has to be
   // saved through multiple passes
   numParal = (N/bin_size)+2;
   goodWindows = (empty_block *)calloc(numParal, sizeof(empty_block));
   eb.notEmpty=0;
   eb.pass_num=1;
   //initialization to 0
   for(g=0; g<numParal; g++)
     goodWindows[g]=eb;
   /**/

   timeStart = cpuTime();
   constructor(N, k, seq, k_factor_ind);
   result = Filter(N, k, p, e, r, bin_size, w, k_factor_ind, seq, goodWindows, passNum); // the last parameter is to say that it is the first pass (it serves for empty block strategy)
   timeStop = cpuTime();

   
   //   if(result->end == N) result->end = N-1; // small bug correction

   result = InvertList(result);

   if(rev) {
     result = ManageResultsForReverse(result,N);
     N = N/2;			/* going back to a sequence with its concatenated rev comp. */
   }
   
     
   if(!multipass){
     // if filtering multi pass is NOT asked for, the info output file is written
     // and filtered sequences are written on output file
     /* OUTPUT INFORMATIONS */
     output = fopen(argv[optind+1], "wt");
     if(output==NULL){fprintf(stderr,"ouput_information_file \"%s\" is not writable, exit\n", argv[optind + 1]); exit(1);}
     fprintf(output,"GIVEN INFORMATIONS:\n"
	     "\tinput file: %s\n"
	     "\tinput size: %d\n"
	     "\twindow length: L=%3d\n"
	     "\tmaximal edit distance between any pair of window: d=%2d\n"
	     "\tminimal number of window: r=%2d\n"
	     "\tk-factors length: %2d\n",  argv[optind], N, w, e, r, k);
     if(rev) fprintf(output,"\tFiltration on direct AND reverse complement strands\n");
     else fprintf(output,"\tFiltration ONLY on direct strand\n");
     
     fprintf(output,"COMPUTED INFORMATIONS:\n"
	     "\tminimal number of k-factors shared by two windows: %3d\n"
	     "\tparallelograms thickness: %2d\n", p, bin_size);
     fprintf(output,"RESULTS:\n");
     kept = WriteResult(N, result, output);
     fprintf(output,"\tpercentage kept: %7.3f%%\n"
	     "\tsize kept: %7d\n"
	     "\ttime (sec): %6.2Lf\n", kept*100, (int)(kept*N+0.5), timeStop-timeStart);
     if (flag) fprintf(output, "\tN kept = %d\n", Nbases);     
     fprintf(output,"\tmethod used: %s\n"
	     "\tversion: %s\n\n",method, version);
     fprintf(output,"Thanks for using Tuiuiu. Please cite\n %s\n", publication);
     fclose(output);
     
     
     /* OUTPUT DATA */
     output2 = fopen(argv[optind + 2],"wt");
     if(output2==NULL){fprintf(stderr,"output_data_file \"%s\" is not writable, exit\n", argv[optind + 2]); exit(1);}
     WriteFilteredSeq(N, seq, result, output2, name, writeN);
     fclose(output2);
     
     
     free(k_factor_ind->ind);
     free(k_factor_ind->list);
     free(k_factor_ind); 
     FreeList(result);
     free(seq);
     free(name);
     free(goodWindows);
   }else{
     // if filtering multi pass is asked for, filter has to be re-run
     // on filtered sequences, so the filtered sequences of the first pass 
     // are not written on output file, but they are used to perform
     // a second (or even more) run of the filter on. The sequences filtered by the
     // second filtering run are finally written in the output file

     
     /* OUTPUT INFORMATIONS first pass*/
     output = fopen(argv[optind+1], "wt");
     if(output==NULL){fprintf(stderr,"ouput_information_file \"%s\" is not writable, exit\n", argv[optind + 1]); exit(1);}
     fprintf(output,"GIVEN INFORMATIONS:\n"
	     "\tinput file: %s\n"
	     "\tinput size: %d\n"
	     "\twindow length: L=%3d\n"
	     "\tmaximal edit distance between any pair of window: d=%2d\n"
	     "\tminimal number of window: r=%2d\n"
	     "\tk-factors length: %2d\n",  argv[optind], N, w, e, r, k);
     if(rev) fprintf(output,"\tFiltration on direct AND reverse complement strands\n");
     else fprintf(output,"\tFiltration ONLY on direct strand\n");
     
     fprintf(output,"\nCOMPUTED INFORMATIONS FOR PASS %d:\n"
	     "\tminimal number of k-factors shared by two windows: %3d\n"
	     "\tparallelograms thickness: %2d\n", passNum, p, bin_size);
     fprintf(output,"RESULTS:\n");
     //kept = WriteResult(N, result, output);
     kept = countResult(N, result);
     fprintf(output,"\tpercentage kept: %7.3f%%\n"
	     "\tsize kept: %7d\n"
	     "\ttime (sec): %6.2Lf\n", kept*100, (int)(kept*N+0.5), timeStop-timeStart);
     if (flag) fprintf(output, "\tN kept = %d\n", Nbases);     

     fprintf(output,"\n");

     long double totalTime = timeStop-timeStart;

     printf("Pass %d\n", passNum);

     if((kept*100)!=0){ // if at least one window is kept after the first pass

       //do{
	 passNum++;
	 
	 //oldKept = kept;
	 
	 // store filtered sequences on which we run filter for more times
	 saveFilteredSeqWithNs(N, &seq, result);
	 
	 free(k_factor_ind->ind);
	 free(k_factor_ind->list);
	 free(k_factor_ind);
	 FreeList(result);
	 
	 k_factor_ind = (index_str *)malloc(sizeof(index_str));
	 
	 // further run of filter
	 timeStart = cpuTime();
	 constructor(N, k, seq, k_factor_ind);
	 result = Filter(N, k, p, e, r, bin_size, w, k_factor_ind, seq, goodWindows, passNum); // the last parameter is to say that it is not the first pass (it serves for empty block strategy)
	 timeStop = cpuTime();
	 
	 //   if(result->end == N) result->end = N-1; // small bug correction
	 
	 result = InvertList(result);
	 
	 if(rev) {
	   result = ManageResultsForReverse(result,N);
	   N = N/2;			/* going back to a sequence with its concatenated rev comp. */
	 } 
	 
	 totalTime = totalTime + (timeStop-timeStart);
	 
	 /* OUTPUT INFORMATIONS multiple passes*/
	 fprintf(output,"\nCOMPUTED INFORMATIONS FOR PASS %d:\n"
		 "\tminimal number of k-factors shared by two windows: %3d\n"
		 "\tparallelograms thickness: %2d\n", passNum, p, bin_size);
	 fprintf(output,"RESULTS:\n");
	 kept = countResult(N, result);
	 fprintf(output,"\tpercentage kept: %7.3f%%\n"
		 "\tsize kept: %7d\n"
		 "\ttime (sec): %6.2Lf\n"
		 "\ttotal time (%d passes) (sec): %6.2Lf\n", 
		 kept*100, (int)(kept*N+0.5), timeStop-timeStart, passNum, totalTime);
	 if (flag) fprintf(output, "\tN kept = %d\n", Nbases);
	 
	 fprintf(output,"\n");
	 
	 printf("Pass %d\n", passNum);
	 //}while(/*(kept < oldKept*multi) &&*/ (kept < oldKept) && ((kept*100)!=0) && (passNum<r));
	 // first condition commented for first tests in order to evaluate param multi
	 
	 // For Todd: you find a commented do..while construct,
	 // because the code was written to perform more passes, as you can see 
	 // from the previous commented line, but we don't evaluate the parameter multi
	 // until now and in performed tests we observed that passes successive to the 
	 // second one do not improve very much the selectiveness of the filter, so for
         // now we decided to stop after only two passes (and in this case, the do--while
         // construct is unuseful) 
     }

     /* TOTAL OUTPUT INFORMATIONS */
     kept = WriteResult(N, result, output);

     fprintf(output,"\nRESULT SUMMARY \n");
     fprintf(output,"\tpercentage kept: %7.3f%%\n"
	     "\tsize kept: %7d\n", kept*100, (int)(kept*N+0.5));
     fprintf(output,
	     "\ttotal time (num passes: %d) (sec): %6.2Lf\n", passNum, totalTime);
     fprintf(output,"\tmethod used: %s\n"
	     "\tversion: %s\n\n",method, version);
       
     // close info output file 
     fprintf(output,"Thanks for using Tuiuiu. Please cite\n %s\n", publication);
     fclose(output);
  
     /* OUTPUT DATA */
     output2 = fopen(argv[optind + 2],"wt");
     if(output2==NULL){fprintf(stderr,"output_data_file \"%s\" is not writable, exit\n", argv[optind + 2]); exit(1);}
     WriteFilteredSeq(N, seq, result, output2, name, writeN);
     fclose(output2);
     
     
     free(k_factor_ind->ind);
     free(k_factor_ind->list);
     free(k_factor_ind); 
     FreeList(result);
     free(seq);
     free(name);
     free(goodWindows);
   }
     
   printf("Computation ended, output information file is : \"%s\" and output data file is \"%s\"\n", argv[optind+1], argv[optind + 2]);
   printf("Thanks for using Tuiuiu. Please cite\n %s\n", publication);
      
   return 0;
}
