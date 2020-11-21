/****************************************************************************

 It computes the set of non dominated scores for the #Matches-#Indels problem.
 It is a Dynamic Programming (DP) without pruning approach, by filling a DP
 matrix - P (with dimensions (M+1)*(N+1), where M and N are the sizes of the
 1st and 2nd sequence respectively).
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
#include <time.h>
#include "pareto_set.h"

#define MAX_LENGTH		(4000)					/* maximum length of the sequences */

void read_sequence(char seq[] , char *filename);
void init_dynamic_table( );
void align( );
void traceback( );
void reverse(char s[]);
void remove_dynamic_table();


/* read sequences */
char seq1[MAX_LENGTH];
char seq2[MAX_LENGTH];

/* DP table */
Pareto_set ** P;

/* sizes of the DP table (M-#lines, N-#colunms) */
int M, N;

int main( int argc, char *argv[] ) {

	srand( time(NULL) );	/* random seed to choose one score between equal score vectores */

	if(argc != 3){
		printf("Usage: %s <seq1_file> <seq2_file>\n", argv[0]);
		return 0;
	}

	read_sequence(seq1, argv[1]);
	read_sequence(seq2, argv[2]);

	M = strlen(seq1);		N = strlen(seq2);
	
	init_dynamic_table(  );	
	align( );
	traceback( );

	remove_dynamic_table();

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
	 * Allocates memory for the all the DP tables and
	 * initialize the its base cases.
	 */

	int i;

	P = (Pareto_set **) malloc( (M+1) * sizeof(Pareto_set *) );

	for(i=0; i<M+1 ; i++)
		P[i] = (Pareto_set *) malloc( (N+1) * sizeof(Pareto_set) );
	
	P[0][0] = create_list( );
	append_node( P[0][0], 0, 0, NULL, '\0' );
	/* base cases */
	for( i = 1; i <= M ; i++ ){
		P[i][0] = create_list();
		append_node( P[i][0], 0, i, P[i-1][0]->next, 'u' );		/* pointing upwards */
	}

	for( i = 1; i <= N ; i++ ){
		P[0][i] = create_list( );
		append_node(P[0][i], 0, i, P[0][i-1]->next, 'l' );		/* pointing leftwards */
	}	
}

void align(  ){
	/**
	 * Filling of the DP table.
	 */	
	
	int i, j, match;

	for( i = 1; i <= M; ++i ){				
		for( j = 1; j <= N; ++j ){
			match = seq1[i-1] == seq2[j-1];
			P[i][j] = pareto_merge( P[i-1][j], P[i-1][j-1], P[i][j-1], match );
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
	int matches, indels;						/* auxiliar variables */
	Node * solution, *tb_ptr;
	char res1[MAX_LENGTH*2], res2[MAX_LENGTH*2];

	solution = P[M][N]->next;
	while( solution ){

		tb_ptr = solution;
		len1 = len2 = 0;
		i = M; 	j = N;							/* start at the bottom right corner of the DP table */

		matches = solution->score.matches;
		indels = solution->score.indels;

		while( i!=0 ||  j!=0 ){					/* while the upper left corner of the DP table isn't reached */

			switch( tb_ptr->traceback_dir ){
				case 'u':						/* going up */
					i--;
					res2[len2++] = '-';
					res1[len1++] = seq1[i];			
				break;
				case 'l':						/* going left */
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
		printf("matches=%d indels=%d\n%s\n%s\n\n", matches, indels, res1, res2);	
		
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

void remove_dynamic_table(){
	/**
	 * Free DP table.
	 */	
	int i, j;

	for( i = 0; i < M+1; ++i )
		for( j = 0; j < N+1; ++j )
			remove_list( P[i][j] );
	

	for( i = 0; i < M+1; ++i ){
		free(P[i]);
	}
	free(P);
}