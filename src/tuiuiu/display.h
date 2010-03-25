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

double countResultMultiSeqs(int N, int *seqBegins, tuilist **result, int nseq);
double WriteResultMultiSeqs(int N, int *seqBegins, tuilist **result, FILE *output,  int nseq, char **name);
void WriteFilteredMultiSeq(int *seqBegins, tuilist **result, char *seq, FILE *output, int nseq, char **name, int flag, int bin_size);
void saveFilteredMultiSeqWithNs(int *seqBegins, tuilist **result, char **seq, int nseq);

double countResult(int N, tuilist *result);
double WriteResult(int N, tuilist *result, FILE *output);
void WriteFilteredSeq(int N, char *seq, tuilist *result, FILE *output, char *name, int flag);
void saveFilteredSeqWithNs(int N, char **seq, tuilist *result);
