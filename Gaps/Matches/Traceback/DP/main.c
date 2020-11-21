/****************************************************************************

 It computes the set of non dominated scores for the #Matches-#Gaps problem.
 It is a Dynamic Programming (DP) without pruning approach, by keeping three
 DP matrices: Q, S and T (with dimensions (M+1)*(N+1), where M and N are the
 sizes of the 1st and 2nd sequence respectively).
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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include "pareto_set.h"

#define MAX_LENGTH		(4000)					/* maximum length of the sequences */

void read_sequence(char seq[] , char *filename);
void init_tables( );
void align( );
void traceback( );
void reverse(char s[]);
void remove_tables();

/* read sequences */
char seq1[MAX_LENGTH];
char seq2[MAX_LENGTH];

/* DP main table */
Pareto_set ** Q;
/* DP table with first sequences ending with '-'		('-',bj ) */
Pareto_set ** S;
/* DP table with second sequences ending with '-'		(ai ,'-') */
Pareto_set ** T;

/* sizes of the DP table (M-#lines, N-#colunms) */
int M, N;

int main( int argc, char *argv[]) {
	
	srand( time(NULL) );	/* random seed to choose one score between equal score vectores */

	if(argc != 3){
		printf("Usage: %s <seq1_file> <seq2_file>\n", argv[0]);
		return 0;
	}

	read_sequence(seq1, argv[1]);
	read_sequence(seq2, argv[2]);

	M = strlen(seq1);		N = strlen(seq2);

	init_tables();
	align();
	traceback();

	remove_tables();

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

void init_tables( ){
	/**
	 * Allocates memory for the DP table and
	 * initialize the base cases.
	 */	

	int i;

	Q = (Pareto_set **) malloc( (M+1) * sizeof(Pareto_set *) );
	S = (Pareto_set **) malloc( (M+1) * sizeof(Pareto_set *) );
	T = (Pareto_set **) malloc( (M+1) * sizeof(Pareto_set *) );

	for(i=0; i<M+1 ; i++){
		Q[i] = (Pareto_set *) malloc( (N+1) * sizeof(Pareto_set) );
		S[i] = (Pareto_set *) malloc( (N+1) * sizeof(Pareto_set) );
		T[i] = (Pareto_set *) malloc( (N+1) * sizeof(Pareto_set) );
	}

	Q[0][0] = create_list( );
	append_node( Q[0][0], 0, 0, NULL, '\0' );
	/* base cases */
	for( i = 1; i <= M; ++i ){
		S[i][0] = create_list(); append_node( S[i][0], 0, INT_MAX, NULL, '-' );			/* against DP table definition (wrong position of gap if traceback is done) */
		Q[i][0] = create_list(); append_node( Q[i][0], 0, 1, Q[i-1][0]->next, 'u' );
	}
	for( i = 1; i <= N; ++i ){
		T[0][i] = create_list(); append_node( T[0][i], 0, INT_MAX, NULL, '-' );			/* against DP table definition (wrong position of gap if traceback is done) */
		Q[0][i] = create_list(); append_node( Q[0][i], 0, 1, Q[0][i-1]->next, 'l' );				
	}
}

void align( ){
	/**
	 * Filling of the main DP table Q by filling the other tables too.
	 */	
	
	int i, j, match;

	for( i = 1; i <= M; ++i ){
		for( j = 1; j <= N; ++j ){
			match = seq1[i-1] == seq2[j-1];

			S[i][j] = pareto_merge2( S[i][j-1], Q[i][j-1], 'l' );
			T[i][j] = pareto_merge2( T[i-1][j], Q[i-1][j], 'u' );
			Q[i][j] = pareto_merge3( S[i][j], T[i][j], Q[i-1][j-1], match );					
		}
	}
}

void traceback( ){
	/**
	 * Traceback (based on the direction and pointer
	 * stored in each Node struct).
	 * Outputs the two aligned sequences to each solution.
	 */

	int i, j;
	/* length of the each solution to each sequence */
	int len1, len2;
	int matches, gaps;							/* auxiliar variables */
	Node *solution, *tb_ptr;
	char res1[MAX_LENGTH*2], res2[MAX_LENGTH*2];

	solution = Q[M][N]->next;
	while( solution ){

		tb_ptr = solution;
		len1 = len2 = 0;
		i = M; 	j = N;							/* start at the bottom right corner of the DP table */

		matches = solution->score.matches;
		gaps = solution->score.gaps;

		while( i>0 ||  j>0 ){					/* while the upper left corner of the DP table isn't reached */

			switch( tb_ptr->traceback_dir ){
				case 'u':						/* going up */
					i--;
					res2[len2++] = '-';
					res1[len1++] = seq1[i];			
				break;
				case 'l':
					j--;
					res1[len1++] = '-';
					res2[len2++] = seq2[j];				
				break;
				case 'd':						/* going up and left */
					i--; j--;
					res1[len1++] = seq1[i];
					res2[len2++] = seq2[j];		
				break;						
			}		

			tb_ptr = tb_ptr->traceback_ptr;					
		}

		res1[len1] = '\0';
		res2[len2] = '\0';		
		/* sequences must be reversed */
		reverse(res1);
		reverse(res2);	
		printf("matches=%d gaps=%d\n%s\n%s\n\n", matches, gaps, res1, res2);			
		
		solution = solution -> next;			/* next solution */
	}	
}

void reverse(char s[]){
	int i, len;
	char tmp;
	len = strlen(s);

	for(i=0; i<len/2; i++){
		tmp = s[i];
		s[i] = s[len-i-1];
		s[len-i-1] = tmp;
	}
}

void remove_tables(){
	/**
	 * Free all DP tables.
	 */

	int i, j;

	/* Q */
	for( i = 0; i < M+1; ++i )
		for( j = 0; j < N+1; ++j )
			remove_list( Q[i][j] );

	for( i = 0; i < M+1; ++i )
		free(Q[i]);
	free(Q);

	/* T */
	for( i = 0; i < M+1; ++i )
		for( j = 1; j < N+1; ++j )
			remove_list( T[i][j] );
	
	for( i = 0; i < M+1; ++i )
		free(T[i]);
	free(T);		

	/* S */
	for( i = 1; i < M+1; ++i )
		for( j = 0; j < N+1; ++j )
			remove_list( S[i][j] );
	
	for( i = 0; i < M+1; ++i )
		free(S[i]);
	free(S);	

}