/*
	<-- RepSeek, a repeat search engine -->
	Some help to compile the repseek program
*/

1. Untar, and go the newly created directory repseek

2. Edit the Makefile, and fill the fields MACHINE, MALLOC, PROTO, INSTALLDIR
   and change any other field you want to.

3. Type 'make repseek' to compile the repseek binary, 'make install' to install the program
   in your ~/bin and 'make clean' to clean all object and the binary.

4. You can create an archive of RepSeek by typing 'make archive'

5. General usage is 'repseek [-v|-h|-H] [ opts ] { core_value } file.fst [file2.fst]'

   Core values (select either -l, -p or -P option)
        -l lmin      : set a minimum length for the seeds detection (no pvalue). Also usuable along with -P.
        -p prob      : set a p-value that gives a minimum seed length (Karlin & Ost equation)
        -P prob      : set a p-value that gives a minimum repeat score (Waterman & Vingron regression)

   Information
        -v           : Print 'v'ersion and exit
        -h           : Print a short 'h'elp and exit
        -H           : Print a larger 'H'elp on a specific mode and exit


6. For more details on repseek insights, please refer to the repseek documentation

7. For License, Repseek is under Lesser GPL. For more details, see the License.txt file

