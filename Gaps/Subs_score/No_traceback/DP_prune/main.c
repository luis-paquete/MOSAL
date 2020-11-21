/****************************************************************************

 It computes the set of non dominated scores for the "Substitution score"-#Gaps problem
 without traceback (only final scores). It is a Dynamic Programming
 (DP) with pruning approach, by keeping three DP matrices: Q, S and T
 (with dimensions (M+1)*(N+1), where M and N are the sizes of the 1st and
 2nd sequence respectively).
 Each entry S[i,j] and T[i,j] will store the set of states corresponding 
 to Pareto optimal alignments of subsequences ending with ('-',b_j) and
 (a_i,'-'), respectively.
 Each entry Q[i,j] will store the states corresponding to Pareto optimal
 alignments of subsequences ending with (a_i,b_j) computed in S[i,j] and T[i,j].

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

/**
 * For the sake of memory, the Dynamic Programming (DP) tables used
 * have only two lines (just like on "Subs score"-#Indels problem), because there is no
 * traceback step (the traceback pointers do not need to be stored). Instead of two lines
 * table S has only one line, because, in each iteration the only cells that
 * are needed are the ones immediately left of the current cell. The other DP tables
 * need cells that are in the line above the current line.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include "pareto_set.h"
#include "bounds.h"

#define MAX_LENGTH		(4000)

void read_sequence(char seq[] , char *filename);
void init_dynamic_tables( );
void init_subs_table( char * filename );
void align( );
int get_score( int i, int j );
void remove_dynamic_tables();

/* read sequences */
char seq1[MAX_LENGTH];
char seq2[MAX_LENGTH];

/* DP main table */
Pareto_set ** Q;
/* DP table with first sequences ending with '-'		('-',bj ) */
Pareto_set * S;			
/* DP table with second sequences ending with '-'		(ai ,'-') */
Pareto_set ** T;

/* substitution score table */
int SS[LEN_ALPHABET][LEN_ALPHABET];

/* sizes of the DP table (M-#lines, N-#colunms) */
int M, N;

/* maximum #states for each DP table cell (not using linked lists) */
int MAX_STATES;

int main( int argc, char *argv[]) {

	if(argc != 5){
		printf("Usage: %s <seq1_file> <seq2_file> <subs_file> <mid_bounds_number>\n", argv[0]);
		return 0;
	}

	read_sequence(seq1, argv[1]);
	read_sequence(seq2, argv[2]);

	M = strlen(seq1);		N = strlen(seq2);

	init_subs_table( argv[3] );
	init_bounds( atoi(argv[4]) );
	init_dynamic_tables();
	align();

	//print scores
	int i, line;
	line = M % 2;	// determine from which line to start printing	
	for( i = 0; i < Q[line][N].num ; ++i )
		printf("%d %d\n", Q[line][N].scores[i].matches, Q[line][N].scores[i].gaps);

	remove_dynamic_tables();
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

void init_dynamic_tables( ){
	/**
	 * Allocates memory for the all the DP tables and
	 * initialize the its base cases.
	 */
	
	int i, j;

	Q = (Pareto_set **) malloc( 2 * sizeof(Pareto_set *) );
	T = (Pareto_set **) malloc( 2 * sizeof(Pareto_set *) );

	for(i=0; i<2 ; i++){
		Q[i] = (Pareto_set *) malloc( (N+1) * sizeof(Pareto_set) );
		T[i] = (Pareto_set *) malloc( (N+1) * sizeof(Pareto_set) );

		/* allocate memory */
		for( j = 0; j < N+1; ++j ){
			Q[i][j].scores = (VMG *) malloc(MAX_STATES * sizeof(VMG));
			Q[i][j].num = 0;
			T[i][j].scores = (VMG *) malloc(MAX_STATES * sizeof(VMG));
			T[i][j].num = 0;			
		}
	}
	S = (Pareto_set *) malloc( (N+1) * sizeof(Pareto_set) );
	for( i = 0; i < N+1; ++i ){
		S[i].scores = (VMG *)malloc( MAX_STATES * sizeof(VMG));
		S[i].num = 0;
	}


	/* base cases */
	append_node( &Q[0][0], 0, 0 );
	for( i = 1; i <= N; ++i ){
		append_node( &T[0][i], 0, INT_MAX );	/* against DP table definition (wrong position of gap if traceback is done) */

		append_node( &Q[0][i], 0, 1 );				
		prune( &Q[0][i], 0, i, 'l' );
	}

	append_node( &S[0], 0, INT_MAX); 	/* against table definition	(wrong position of gap if traceback is done) */
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
				SS[ alphabet[i]-'*' ][ alphabet[j]-'*' ] = atoi(strtok(NULL , sep));
		}
		free(f_contents);
		fclose(f);
	}
	else{
		fprintf(stderr, "Failed to open subsitution score file: %s\n", strerror(errno));
		exit(-1);		
	}
}

void align( ){
	/**
	 * Filling of the main DP table Q by filling the other tables too.
	 * For simplicity, two pointers were created:
	 * 'curr' and 'prev' that point to the 
	 * current and previous line of the DP tables,
	 * because there are only two lines (or one line).
	 */	

	int i, j, k, subs_score;
	int prev, curr;		/* curr - current line in the matrix */
						/* prev - previous line */

	for( i = 1; i <= M; ++i ){

		/* exchange lines for tables T and Q */
		prev = !( i % 2);
		curr = (i % 2);

		/* first column of Q is a base case */
		append_node( &Q[curr][0], 0, 1);
		prune( &Q[curr][0], i, 0, 'u' );
		/* first column of T is out of reach (does not need base case) */
		/* first column of S is a base case but was calculated before once in init_dynamic_tables() ) */


		for( j = 1; j <= N; ++j ){

			subs_score = get_score(i-1,j-1);

			pareto_merge2( &S[j], S[j-1], Q[curr][j-1], i, j, 'l' );
			pareto_merge2( &T[curr][j], T[prev][j], Q[prev][j], i, j, 'u' );
			pareto_merge3( &Q[curr][j], S[j], T[curr][j], Q[prev][j-1], subs_score, i, j );					
		}

		/* reset previous lines of Q and T to store new scores */
		for( k = 0; k < N+1; ++k ){
			Q[prev][k].num = 0;
			T[prev][k].num = 0;
		}
		/* reset previous lines of S to store new scores */
		for( k = 0; k < N+1; ++k ){
			S[k].num = 0;
		}
	}
}

int get_score( int i, int j ){
	/**
	 * Returns the substitution score associated with the letters
	 * at position i in seq1 and j in seq2.
	 */
	return SS[ seq1[i]-'*' ][ seq2[j]-'*' ];
}

void remove_dynamic_tables(){
	/**
	 * Free all DP tables.
	 */
	int i, j;

	/* Q */
	for( i = 0; i < 2; ++i )
		for( j = 0; j < N+1; ++j )
			free( Q[i][j].scores );

	for( i = 0; i < 2; ++i )
		free(Q[i]);
	free(Q);


	/* T */
	for( i = 0; i < 2; ++i )
		for( j = 1; j < N+1; ++j )
			free( T[i][j].scores );	
			
	for( i = 0; i < 2; ++i )
		free(T[i]);
	free(T);


	/* S */
	for( i = 0; i < N+1; ++i )
		free( S[i].scores );
	free(S);

}
