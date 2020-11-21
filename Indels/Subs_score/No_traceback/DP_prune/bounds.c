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

/**
 * Implementation of the functions in 'bounds.h'.
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "bounds.h"

void init_bounds( int mid_bounds ){
	MAX_MID_BOUNDS_NUM = mid_bounds;
	BOUNDS_NUM     = mid_bounds + 2;	

	init_lower_bound( );
	calculate_lower_bound( );

	init_upper_bound( );
	calculate_upper_bound( );
}

void init_lower_bound( ){
	int i, j;

	LB = (VMD *) malloc( BOUNDS_NUM * sizeof(VMD) );

	/* L malloc (MAX and MIN) */		
	L = (VMD **) malloc( (2)*sizeof(VMD *) );
	for( i = 0; i < 2; ++i )
		L[i] = (VMD *) malloc( (N+1)*sizeof(VMD) );

	/* base cases */
	L[0][0].matches = 0;
	L[0][0].indels  = 0;

	for( j = 1; j < N+1; ++j ){
		L[0][j].matches = 0;
		L[0][j].indels  = j;
	}

	/* M's mallocs (MID's) */
	M_ = (VMD **) malloc( (2)*sizeof(VMD *) );
	for( i = 0; i < 2; ++i )
		M_[i] = (VMD *) malloc( (N+1)*sizeof(VMD) );
	
}

void calculate_lower_bound( ){
	int i;

	calc_MIN();
	calc_MID();
	calc_MAX();

	/* remove tables (not needed anymore) */
	for( i = 0; i < 2; ++i )
		free( L[i] );
	
	free(L);

	for( i = 0; i < 2; ++i )
		free( M_[i] );

	free( M_ );	

	/*  calculates the maximum number of scores vectors per cell
		(used to allocate memory for the main DP table)
	*/
	MAX_STATES = min2( MAX->matches-MIN->matches, MAX->indels-MIN->indels ) + 1;	

}

void calc_MAX(){
	/* LexMDP */

	int i, j, score;
	int prev, curr;		/* curr - current line in the matrix */
						/* prev - previous line */

	/* base cases */
	L[0][0].matches = 0;
	L[0][0].indels  = 0;
	for( j = 1; j < N+1; ++j ){
		L[0][j].matches = 0;
		L[0][j].indels  = j;
	}

	for( i = 1; i < M+1; ++i ){

		/* exchange lines */
		prev = !( i % 2);
		curr = (i % 2);

		/* base cases */
		L[curr][0].matches = 0;
		L[curr][0].indels  = i;

		for( j = 1; j < N+1; ++j ){
			score = get_score(i-1,j-1);
			L[curr][j] = lexmaxM( L[prev][j], L[curr][j-1], L[prev][j-1], score );
		}
	}
	LB[ BOUNDS_NUM-1 ] = L[curr][N];

	if( LB[ BOUNDS_NUM-1 ].matches == LB[ BOUNDS_NUM-2 ].matches )		/* if MAX is equal to last MID point */
		BOUNDS_NUM--;


	MAX = &LB[ BOUNDS_NUM-1 ] ;
	// printf("MAX=(%d,%d)\n", MAX->matches, MAX->indels );
}

void calc_MIN(){
	/* LexDMP */

	int i, j, score;
	int prev, curr;		/* curr - current line in the matrix */
						/* prev - previous line */


	for( i = 1; i < M+1; ++i ){

		/* exchange lines */
		prev = !( i % 2);
		curr = (i % 2);

		/* base cases */
		L[curr][0].matches = 0;
		L[curr][0].indels  = i;

		for( j = 1; j < N+1; ++j ){
			score = get_score(i-1,j-1);
			L[curr][j] = lexminD( L[prev][j], L[curr][j-1], L[prev][j-1], score );
		}
	}
	LB[0] = L[curr][N];
	MIN = &LB[0];
	//printf("MIN=(%d,%d)\n", MIN->matches, MIN->indels);
}

void calc_MID(){
	int m, i, j, ws, wd, score;

	int prev, curr;		/* curr - current line in the matrix */
						/* prev - previous line */

	int last_mid_matches = MIN->matches;		 /* discard same MID points */
	int mid_bounds_num = 0;	

	for( m = 0; m < MAX_MID_BOUNDS_NUM; ++m ){

		/* base cases */
		M_[0][0].matches = 0;
		M_[0][0].indels  = 0;

		for( j = 1; j < N+1; ++j ){
			M_[0][j].matches = 0;
			M_[0][j].indels  = j;
		}		


		for( i = 1; i < M+1; ++i ){

			/* exchange lines for tables M_ and MT */
			prev = !( i % 2);
			curr = (i % 2);

			/* base cases */
			M_[curr][0].matches = 0;
			M_[curr][0].indels  = i;


			for( j = 1; j < N+1; ++j ){
				
				score = get_score(i-1,j-1);
				ws = m + 1;
				wd = MAX_MID_BOUNDS_NUM - m;

				M_[curr][j] = scalar_max( M_[prev][j], M_[curr][j-1], M_[prev][j-1], score, ws, wd );

			}
		}


		if( last_mid_matches != M_[curr][N].matches ){		/* if MID does not exist yet */
			LB[1 + mid_bounds_num++] = M_[curr][N];
		}
		last_mid_matches = M_[curr][N].matches;
		/* printf("MID=(%d,%d)\n", M_[curr][N].matches, M_[curr][N].indels ); */
	}

	/* update value */
	BOUNDS_NUM = 2 + mid_bounds_num;

}

void init_upper_bound( ){
	int i;
	
	U = (VMD **) malloc( (M+1)*sizeof(VMD *) );
	for( i = 0; i <= M; ++i )
		U[i] = (VMD *) malloc( (N+1)*sizeof(VMD) );
}

void calculate_upper_bound( ){
	/**
	 * LCS variation adapted to substitution scores
	 * in the reversed sequences.
	 */
	
	int i, j ;
	/* base cases */
	U[M][N].matches = 0;
	U[M][N].indels  = 0;
	for( i = 0; i < M; ++i ){
		U[i][N].matches = 0;
		U[i][N].indels  = M-i ;
	}
	for( j = 0; j < N; ++j ){
		U[M][j].matches = 0;
		U[M][j].indels  = (N-j) ;
	}	
	for( i = M-1; i >= 0; --i ){			/* begins in the oposite corner of the table (reversed sequences) */
		for( j = N-1; j >= 0; --j ){
			U[i][j].matches = max3(U[i+1][j+1].matches+get_score(i,j), U[i+1][j].matches , U[i][j+1].matches );
			U[i][j].indels  = abs( (M-i)-(N-j) ) ;
		}
	}
}

VMD lexmaxM( VMD up, VMD left, VMD diagonal, int score ){
	/* match priority */
	VMD result;
	int mu, ml, md, iu, il, id;			/* auxiliar variables */

	mu = up.matches;
	ml = left.matches;
	md = diagonal.matches + score;
	
	iu = up.indels + 1;
	il = left.indels + 1;
	id = diagonal.indels;

	if( mu > ml ){
		if( mu > md ){
			result.matches = mu;
			result.indels  = iu;
		}
		else if( mu == md ){
			if( iu <= id ){
				result.matches = mu;
				result.indels  = iu;	
			}
			else{
				result.matches = md;
				result.indels  = id;				
			}
		}
		else{
			result.matches = md;
			result.indels  = id;			
		}
	}
	else if( mu == ml ){
		if( mu > md ){
			if( iu <= il ){
				result.matches = mu;
				result.indels  = iu;
			}
			else{
				result.matches = ml;
				result.indels  = il;				
			}
		}
		else if( mu == md ){
			if( iu <= il ){
				if( iu <= id ){
					result.matches = mu;
					result.indels  = iu;					
				}
				else{
					result.matches = md;
					result.indels  = id;					
				}
			}
			else{
				if( il <= id ){
					result.matches = ml;
					result.indels  = il;					
				}
				else{
					result.matches = md;
					result.indels  = id;					
				}
			}
		}
		else{			/* mu < md */
			result.matches = md;
			result.indels  = id;			
		}
	}
	else{				/* mu < ml */
		if( ml > md ){
			result.matches = ml;
			result.indels  = il;			
		}
		else if( ml == md){
			if( il <= id ){
				result.matches = ml;
				result.indels  = il;				
			}
			else{
				result.matches = md;
				result.indels  = id;			
			}
		}
		else{
			result.matches = md;
			result.indels  = id;						
		}
	}

	return result;
}


VMD lexminD( VMD up, VMD left, VMD diagonal, int score ){
	/* indels priority */

	VMD result;
	int mu, ml, md, iu, il, id;			/* auxiliar variables */

	mu = up.matches;
	ml = left.matches;
	md = diagonal.matches + score;
	
	iu = up.indels + 1;
	il = left.indels + 1;
	id = diagonal.indels;

	if( iu < id ){
		if( iu < il ){
			result.matches = mu;
			result.indels  = iu;			
		}
		else if( iu == il && mu >= ml ){
			result.matches = mu;
			result.indels  = iu;			
		}
		else{
			result.matches = ml;
			result.indels  = il;							
		}
	}

	else if( iu == id ){

		if( il < iu ){
			result.matches = ml;
			result.indels  = il;
		}
		else if( iu == il ){
			if( mu >= md ){
				if( mu >= ml ){
					result.matches = mu;
					result.indels  = iu;
				}
				else{
					result.matches = ml;
					result.indels  = il;
				}
			}
			else{
				if( md >= ml ){
					result.matches = md;
					result.indels  = id;
				}
				else{
					result.matches = ml;
					result.indels  = il;
				}
			}
		}
		else{			/* iu==id && iu < il */
			if( mu >= md ){
				result.matches = mu;
				result.indels  = iu;
			}
			else{
				result.matches = md;
				result.indels  = id;
			}
		}
	}

	else{

		if( id < il ){
			result.matches = md;
			result.indels  = id;
		}
		else if( id == il && md >= ml ){
			result.matches = md;
			result.indels  = id;			
		}
		else{
			result.matches = ml;
			result.indels  = il;								
		}
	}	


	return result;
}

VMD scalar_max( VMD up, VMD left, VMD diagonal, int score, int ws, int wd ){
	VMD result;

	/* auxiliar variables */
	int iu, il, id;
	int mu, ml, md;
	int ru, rl, rd;		/* temporary results */

	iu = up.indels + 1;
	il = left.indels + 1;
	id = diagonal.indels;

	mu = up.matches;
	ml = left.matches;
	md = diagonal.matches + score;

	ru = ws * mu - wd * iu;
	rl = ws * ml - wd * il;
	rd = ws * md - wd * id;

	if( ru >= rl ){
		if( ru >= rd ){
			result.matches = mu;
			result.indels  = iu;
		}
		else{
			result.matches = md;
			result.indels  = id;
		}
	}
	else{
		if( rl >= rd ){
			result.matches = ml;
			result.indels  = il;
		}
		else{
			result.matches = md;
			result.indels  = id;
		}
	}

	return result;
}

void prune( Pareto_set *set, int i, int j ){
	/**
	 * If the current score vector plus the corresponding upper bound
	 * score vector is dominated by any of the lower bound scores, than
	 * it is pruned.
	 */
	int k;

	/* auxiliar variables 'm', 'd' */
	int m = set->scores[ set->num-1 ].matches + U[i][j].matches;
	int d = set->scores[ set->num-1 ].indels + U[i][j].indels;

	k = 0;
    while( k < BOUNDS_NUM && LB[k].matches < m ) k++;	/* advance in bound scores */

    if( k < BOUNDS_NUM ){
        if( m == LB[k].matches ){
            if( d > LB[k].indels ){
                set->num--;        
            }    
        }
        else{
            if( d >= LB[k].indels ){
                set->num--;        
            }
        }
    }
    else{	/* advanced all bounds (this should never happen) */
        set->num--;
	}
}

int min2(int n1, int n2){
	return n1 <= n2 ? n1 : n2;
}

int max3(int n1, int n2, int n3){
	if(n1 >= n2)
		return n1 >= n3 ? n1 : n3;
	return n2 >= n3 ? n2 : n3;
}

void remove_bounds_tables(){
	int i;

	free(LB);

	for( i = 0; i < M+1; ++i ){
		free( U[i] );
	}
	free(U);
}