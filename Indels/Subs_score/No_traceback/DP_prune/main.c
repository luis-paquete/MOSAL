/****************************************************************************

 It computes the set of non dominated scores for the "Substitution score"-#Indels
 problem without traceback (only final scores). It is a Dynamic Programming
 (DP) with pruning approach, by filling a DP matrix - P
 (with dimensions (M+1)*(N+1), where M and N are the sizes of the 1st and
 2nd sequence respectively).
 Each entry P[i,j] will store the set of states corresponding to Pareto
 optimal alignments of subsequences (a_1,...,a_i) and (b_1,...,b_j).

 ---------------------------------------------------------------------

    Copyright (c) 2013

		Maryam Abassi    (maryam@dei.uc.pt)
		Luís Paquete     (paquete@dei.uc.pt)
		Arnaud Liefooghe (arnaud.liefooghe@univ-lille1.fr)
		Miguel Pinheiro  (monsanto@biocant.pt)
		Pedro Matias     (pamatias@student.dei.uc.pt)

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA

****************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "pareto_set.h"
#include "bounds.h"

#define MAX_LENGTH		(4000)				/* maximum length of the sequences */

void read_sequence(char seq[] , char *filename);
void init_dynamic_table( );
void init_subs_table( char * filename );
void align( );
int get_score( int i, int j );
void remove_dynamic_table();


/* read sequences */
char seq1[MAX_LENGTH];
char seq2[MAX_LENGTH];


/* DP table */
Pareto_set ** P;

/* substitution score table */
int S[LEN_ALPHABET][LEN_ALPHABET];

/* sizes of the DP table (M-#lines, N-#colunms) */
int M, N;

/* maximum #states for each DP table cell (not using linked lists) */
int MAX_STATES;


/**
 * For the sake of memory, the Dynamic Programming (DP) table used
 * has only two lines, because there is no traceback step (the traceback
 * pointers do not need to be stored). Also, in each iteration, the
 * only cells that are needed for the algorithm are the ones that
 * are immediately left, up and diagonal(up and left) of the current cell.
 */
int main( int argc, char ** argv ) {

	if(argc != 5){
		printf("Usage: %s <seq1_file> <seq2_file> <subs_file> <mid_bounds_number>\n", argv[0]);
		return 0;
	}

	read_sequence(seq1, argv[1]);
	read_sequence(seq2, argv[2]);

	M = strlen(seq1);		N = strlen(seq2);

	init_subs_table( argv[3] );
	init_bounds( atoi(argv[4]) );
	init_dynamic_table();	
	align( );

	/* print the result scores */
	int i, line;
	line = M % 2;	/* determine from which line to start printing */ 
	for( i = 0; i < P[line][N].num; i++ )
		printf("%d %d\n", P[line][N].scores[i].matches , P[line][N].scores[i].indels);

	remove_dynamic_table();
	remove_bounds_tables();

	return 0;
}

void read_sequence(char seq[], char *filename){
	/**
	 * Reads the sequence contained in 'filename'.
	 * This file must follow the FASTA format.
	 */
	char fasta_header[300];
	char line[300];

	FILE *f = fopen(filename, "r");
	if(f){
		if( fscanf(f, ">%[^\n]\n", fasta_header)>0 ){	/* determine if file is valid */
			while( fscanf(f, "%s", line)!=EOF )
				strcat(seq, line);			
		}
		else{
			printf("Sequence file format is not correct!\n");
			exit(-1);
		}
		fclose(f);
	}
	else{
		fprintf(stderr, "Failed to open sequence file: %s\n", strerror(errno));
		exit(-1);
	}
}

void init_dynamic_table( ){
	/**
	 * Allocates memory for the DP table and
	 * initialize the base cases.
	 */	
	int i, j;

	P = (Pareto_set **) malloc( 2 * sizeof(Pareto_set *) );

	for(i=0; i < 2 ; i++){
		P[i] = (Pareto_set *) malloc( (N+1) * sizeof(Pareto_set) );

		/* allocate memory */
		for( j = 0; j < N+1; ++j ){
			P[i][j].scores = (VMD *) malloc( MAX_STATES * sizeof(VMD) );
			P[i][j].num = 0;
		}
	}

	/* init top line (base cases) */
	for( i = 0; i <= N ; i++ ){
		append_node( &P[0][i], 0, i);
		prune( &P[0][i], 0, i );
	}
}

void init_subs_table( char * filename ){
	/**
	 * Builds substitution table from files that have a table with the format:
	 * 			A1  A2 ... An  *
	 *   	A1   -  -   -   -  -
	 *    	A2   -  -   -   -  -
	 *     	...  -  -   -   -  -
	 *      An   -  -   -   -  -
	 *      *    -  -   -   -  -

	 * -> the Ai is the ith letter of the alphabet and each in each cell is an
	 * integer with the correspnding substitution score of the letters crossing
	 * -> The * column is optional and corresponds to the minimum score
	 */

	FILE *f;
	char *f_contents, *it, *A, c, *sep = " \r\n\t", alphabet[LEN_ALPHABET];
	int i, j, f_size, len_alphabet = 0, result;

	f = fopen(filename, "r");
	if(f){
		/* obtain file size */
		fseek(f, 0 , SEEK_END);
		f_size = ftell(f);
		rewind(f);

		/* obtain memory */
		f_contents = (char *)malloc(f_size * sizeof(char));
		if(!f_contents){ fprintf(stderr, "Failed to allocate memory: %s\n", strerror(errno)); exit(-1);}

		result = fread(f_contents, sizeof(char), f_size, f);
		if(result != f_size){fprintf(stderr, "Failed to load file: %s\n", strerror(errno)); exit(-1);}

		/* skip comments */
		it = f_contents;
		while(*it=='#')	it = strchr(it, '\n')+1;

		/* read alphabet */
		while(*it != '\n'){
			c = *it++;
			if((c>='A' && c<='Z') || (c=='*')){
				alphabet[ len_alphabet++ ] = c;
			}
		}
		it++;

		/* read scores */
		for (i = 0; i < len_alphabet; ++i){
			A = strtok(i==0 ? it:NULL, sep);
			for (j = 0; j < len_alphabet; ++j)
				S[ alphabet[i]-'*' ][ alphabet[j]-'*' ] = atoi(strtok(NULL , sep));
		}
		free(f_contents);
		fclose(f);
	}
	else{
		fprintf(stderr, "Failed to open subsitution score file: %s\n", strerror(errno));
		exit(-1);		
	}
}

void align(  ){
	/**
	 * Filling of the DP table.
	 * For simplicity, two pointers were created:
	 * 'curr' and 'prev' that point to the 
	 * current and previous line of the DP table,
	 * because there are only two lines.
	 */	
	int i, j, k, subs_score;
	int prev, curr;		/* curr - current line in the matrix */
						/* prev - previous line */

	for( i = 1; i <= M; ++i ){				

		/* exchange lines */
		prev = !(i % 2);
		curr = (i % 2);

		/* first column is a base case */
		append_node( &P[curr][0], 0, i);
		prune( &P[curr][0], i, 0 );

		for( j = 1; j <= N; ++j ){
			subs_score = get_score(i-1,j-1);
			pareto_merge( &P[curr][j], P[prev][j], P[prev][j-1], P[curr][j-1], subs_score, i, j );
		}
		
		/* reset previous line to store new scores */
		for( k = 0; k <= N; ++k )
			P[prev][k].num = 0;
		
	}
}

int get_score( int i, int j ){
	/**
	 * Returns the substitution score associated with the letters
	 * at position i in seq1 and j in seq2.
	 */
	return S[ seq1[i]-'*' ][ seq2[j]-'*' ];
}

void remove_dynamic_table(){
	/**
	 * Free DP table.
	 */
	int i, j;

	for( i = 0; i < 2; ++i ){
		for( j = 0; j < N+1; ++j )
			free( P[i][j].scores );
	}

	for( i = 0; i < 2; ++i ){
		free(P[i]);
	}
	free(P);
}