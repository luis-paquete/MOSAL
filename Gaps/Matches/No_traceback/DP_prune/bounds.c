/****************************************************************************

 It computes the set of non dominated scores for the #Matches-#Gaps problem
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
 * Implementation of the functions in 'bounds.h'.
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "bounds.h"


void init_bounds( int mid_bounds ){
	MAX_MID_BOUNDS_NUM = mid_bounds;
	BOUNDS_NUM         = mid_bounds + 2;

	init_lower_bound( );
	calculate_lower_bound( );

	init_upper_bound( );
	calculate_upper_bound( );
}

void init_lower_bound( ){
	int i, j;

	LB = (VMG *) malloc( BOUNDS_NUM * sizeof(VMG) );

	/* L and TL malloc (MAX and MIN) */		
	L  = (VMG **) malloc( (2)*sizeof(VMG *) );			/* just 2 lines (just left, upper and upper-left cells are needed) */
	LT  = (VMG **) malloc( (2)*sizeof(VMG *) );
	for( i = 0; i < 2; ++i ){
		L[i]  = (VMG *) malloc( (N+1)*sizeof(VMG) );
		LT[i]  = (VMG *) malloc( (N+1)*sizeof(VMG) );
	}
	LS = (VMG *) malloc( (N+1)*sizeof(VMG) );
	
	TL = (char **) malloc( (M+1)*sizeof(char *) );
	for( i = 0; i < M+1; ++i ){
		TL[i] = (char *) malloc( (N+1)*sizeof(char) );
	}
	
	/* base cases */
	L[0][0].matches = 0;
	L[0][0].gaps    = 0;
	
	L[1][0].matches = 0;
	L[1][0].gaps    = 1;
	
	for( i = 0; i < M+1; ++i ){
		TL[i][0]        = 'u';
	}
	for( j = 1; j < N+1; ++j ){
		L[0][j].matches = 0;
		L[0][j].gaps    = 1;

		LT[0][j].matches = 0;
		LT[0][j].gaps    = INT_MAX;		

		TL[0][j]        = 'l';
	}

	LS[0].matches = 0;
	LS[0].gaps    = INT_MAX;

	/* M's mallocs (MID's) */
	M_ = (VMG **) malloc( 2*sizeof(VMG *) );
	MT = (VMG **) malloc( 2*sizeof(VMG *) );
	for( i = 0; i < 2; ++i ){
		M_[i] = (VMG *) malloc( (N+1)*sizeof(VMG) );
		MT[i] = (VMG *) malloc( (N+1)*sizeof(VMG) );
	}

	MS = (VMG *) malloc( (N+1)*sizeof(VMG) );
	MS[0].matches = 0;
    MS[0].gaps    = INT_MAX;
}

void calculate_lower_bound( ){
	int i;

	calc_MIN();
	calc_MID();
	calc_MAX();

	/* remove tables (not needed anymore) */
	for( i = 0; i < 2; ++i ){
		free( L[i] );
		free( LT[i] );
	}
	free(L);
	free(LT);
	free(LS);
	for( i = 0; i < M+1; ++i ){
		free( TL[i] );
	}
	free(TL);

	for( i = 0; i < 2; ++i ){
		free( M_[i] );
		free( MT[i] );
	}
	free( M_ );
	free( MS );
	free( MT );

	/*  calculates the maximum number of scores vectors per cell
		(used to allocate memory for the main DP table)
	*/
	MAX_STATES = min2( MAX->matches-MIN->matches, MAX->gaps-MIN->gaps ) + 1;		

}

void calc_MAX(){
	/* LexMGP */
	
	int i, j, match;
	int prev, curr;		/* curr - current line in the matrix */
						/* prev - previous line */

	for( i = 1; i < M+1; ++i ){

		/* exchange lines */
		prev = !( i % 2);
		curr = (i % 2);

		L[curr][0].matches = 0;
		L[curr][0].gaps    = 1;
		/* first column of LS is a base case but was calculated before once in init_lower_bound() ) */
		/* first column of LT is out of reach (does not need base case) */

		for( j = 1; j < N+1; ++j ){
			match = seq1[i-1] == seq2[j-1];

			LS[j]       = lexmaxMS( LS[j-1], L[curr][j-1], i, j );
			LT[curr][j] = lexmaxMT( LT[prev][j], L[prev][j], i, j );
			L[curr][j]  = lexmaxM3( LS[j], LT[curr][j], L[prev][j-1], match, i, j );

		}
	}
	
	LB[ BOUNDS_NUM-1 ] = L[curr][N];

	if( LB[ BOUNDS_NUM-1 ].matches == LB[ BOUNDS_NUM-2 ].matches )		/* if MAX is equal to last MID point */
		BOUNDS_NUM--;

	MAX = &LB[ BOUNDS_NUM-1 ] ;
	/* printf("MAX=(%d,%d)\n", MAX->matches, MAX->gaps ); */

}
void calc_MIN(){
	/**
	 * Instead of filling DP tables, the MIN score vector can be determined
	 * in a easier way:
	 * 	  - the minimum #gaps is either 0 (M=N) or 1 (M!=N, indels all together)
	 * 	  - the corresponding #matches can then be easily determined based on
	 * 	  	all possibilities of the gap position in the smaller sequence
	 */
	VMG min;
	int i, match, gap_size;
	char * greater_seq, *smaller_seq;

	match = 0;
	min.gaps = M == N ? 0 : 1;

	if( min.gaps == 0 ){
		for( i = 0; i < N; ++i )
			match += seq1[i] == seq2[i];		
		min.matches = match;
	}
	else{
		gap_size = abs( M-N );

		/* determine greater sequence */
		if( M > N )	{
			greater_seq = seq1;
			smaller_seq = seq2;
		}
		else {
			greater_seq = seq2;
			smaller_seq = seq1;
		}

		/* gap is at the beginning */
		for( i = gap_size; i < strlen(greater_seq); ++i )
			match += greater_seq[i] == smaller_seq[i-gap_size];		

		min.matches = match;

		/* 	iterate through all gap positions that are possible
			by advancing the gap one letter in each iteration in the 
			smaller sequence */
		for( i = 0; i < strlen(smaller_seq); ++i ){
			
			match += greater_seq[i] == smaller_seq[i];
			match -= greater_seq[i+gap_size] == smaller_seq[i];

			/* update minimum #matches */
			if( match > min.matches)
				min.matches = match;
		}
	}
	/* printf("MIN=(%d,%d)\n", min.matches, min.gaps); */
	LB[0] = min;
	MIN = &LB[0];
}

void calc_MID(){
	int m, i, j, ws, wg, match;

	int prev, curr;		/* curr - current line in the matrix */
						/* prev - previous line */

	int last_mid_matches = MIN->matches;		 /* discard same MID points */
	int mid_bounds_num = 0;

	for( m = 0; m < MAX_MID_BOUNDS_NUM; ++m ){

		/* base cases */
		M_[0][0].matches = 0;
		M_[0][0].gaps    = 0;

		M_[1][0].matches = 0;
		M_[1][0].gaps    = 1;

		for( j = 1; j < N+1; ++j ){
			M_[0][j].matches = 0;
			M_[0][j].gaps    = 1;

			MT[0][j].matches = 0;
			MT[0][j].gaps    = INT_MAX;
		}		

		for( i = 1; i < M+1; ++i ){

			/* exchange lines for tables M_ and MT */
			prev = !( i % 2);
			curr = (i % 2);

			/* base cases */
			M_[curr][0].matches = 0;
			M_[curr][0].gaps    = 1;
			/* first column of MS is a base case but was calculated before once in init_lower_bound() ) */
			/* first column of MT is out of reach (does not need base case) */

			for( j = 1; j < N+1; ++j ){
				
				match = seq1[i-1] == seq2[j-1];
				ws = m + 1;
				wg = MAX_MID_BOUNDS_NUM - m;

				MS[j] = scalar_max2( MS[j-1], M_[curr][j-1], ws, wg );
				MT[curr][j] = scalar_max2( MT[prev][j], M_[prev][j], ws, wg );
				M_[curr][j] = scalar_max3( MS[j], MT[curr][j], M_[prev][j-1], match, ws, wg );

			}
		}

		if( last_mid_matches != M_[curr][N].matches ){		/* if MID does not exist yet */
			LB[1 + mid_bounds_num++] = M_[curr][N];
		}
		last_mid_matches = M_[curr][N].matches;
		/* printf("MID=(%d,%d)\n", M_[curr][N].matches, M_[curr][N].gaps ); */
	}

	/* update value */
	BOUNDS_NUM = 2 + mid_bounds_num;

}

void init_upper_bound( ){
	int i;

	U = (int **) malloc( (M+1)*sizeof(int *) );
	for( i = 0; i <= M; ++i )
		U[i] = (int *) malloc( (N+1)*sizeof(int) );
	
}

void calculate_upper_bound( ){
	/**
	 * Based on the classical DP algorithm for computing the LCS
	 * in the reversed sequences.
	 */
	int i, j ;
	/* base cases */
	U[M][N] = 0;
	for( i = 0; i < M; ++i )
		U[i][N] = 0;
	
	for( j = 0; j < N; ++j )
		U[M][j] = 0;
	
	for( i = M-1; i >= 0; --i )				/* begins in the oposite corner of the table (reversed sequences) */
		for( j = N-1; j >= 0; --j )
			if( seq1[i] == seq2[j] )
				U[i][j] = U[i+1][j+1] + 1;			
			else
				U[i][j] = max2(U[i+1][j] , U[i][j+1]);
}

VMG lexmaxMS( VMG p1, VMG p2, int i, int j ){
	/* match priority */
	
	VMG result;

	/* auxiliar variables */
	int m1, m2, g1, g2;

	m1 = p1.matches;
	m2 = p2.matches;

	g1 = p1.gaps + IND_PENALTY;				/* add additional indel */
	g2 = p2.gaps + (TL[i][j-1] != 'l') ;	/* adds gaps only if it doesn't already exist */	

	if( m1 > m2 || (m1==m2 && g1 <= g2) ){
		result.matches = m1;
		result.gaps    = g1;
	}
	else{
		result.matches = m2;
		result.gaps    = g2;
	}

	return result;
}

VMG lexmaxMT( VMG p1, VMG p2, int i, int j ){
	/* match priority */
	
	VMG result;

	/* auxiliar variables */
	int m1, m2, g1, g2;

	m1 = p1.matches;
	m2 = p2.matches;

	g1 = p1.gaps + IND_PENALTY;				/* add additional indel */
	g2 = p2.gaps + (TL[i-1][j] != 'u');		/* adds gaps only if it doesn't already exist */	

	if( m1 > m2 || (m1==m2 && g1 <= g2) ){
		result.matches = m1;
		result.gaps    = g1;
	}
	else{
		result.matches = m2;
		result.gaps    = g2;
	}

	return result;
}


VMG lexmaxM3( VMG S, VMG T, VMG Q, int match, int i, int j){
	/* match priority */

	VMG result;
	int mt, ms, mq, gt, gs, gq;

	mt = T.matches;
	ms = S.matches;
	mq = Q.matches + match;
	
	gt = T.gaps;
	gs = S.gaps;
	gq = Q.gaps;

	if( mt > ms ){
		if( mt > mq ){
			result.matches = mt;
			result.gaps    = gt;
			TL[i][j]       = 'u';
		}
		else if( mt == mq ){
			if( gt <= gq ){
				result.matches = mt;
				result.gaps    = gt;	
				TL[i][j]       = 'u';
			}
			else{
				result.matches = mq;
				result.gaps    = gq;	
				TL[i][j]       = 'd';
			}
		}
		else{
			result.matches = mq;
			result.gaps    = gq;	
			TL[i][j]       = 'd';		
		}
	}
	else if( mt == ms ){
		if( mt > mq ){
			if( gt <= gs ){
				result.matches = mt;
				result.gaps    = gt;
				TL[i][j]       = 'u';
			}
			else{
				result.matches = ms;
				result.gaps    = gs;	
				TL[i][j]       = 'l';
			}
		}
		else if( mt == mq ){
			if( gt <= gs ){
				if( gt <= gq ){
					result.matches = mt;
					result.gaps    = gt;	
					TL[i][j]       = 'u';				
				}
				else{
					result.matches = mq;
					result.gaps    = gq;
					TL[i][j]       = 'd';					
				}
			}
			else{
				if( gs <= gq ){
					result.matches = ms;
					result.gaps    = gs;	
					TL[i][j]       = 'l';				
				}
				else{
					result.matches = mq;
					result.gaps    = gq;	
					TL[i][j]       = 'd';				
				}
			}
		}
		else{			/* mt < mq */
			result.matches = mq;
			result.gaps    = gq;			
			TL[i][j]       = 'd';
		}
	}
	else{				/* mt < ms */
		if( ms > mq ){
			result.matches = ms;
			result.gaps    = gs;		
			TL[i][j]       = 'l';	
		}
		else if( ms == mq){
			if( gs <= gq ){
				result.matches = ms;
				result.gaps    = gs;		
				TL[i][j]       = 'l';		
			}
			else{
				result.matches = mq;
				result.gaps    = gq;	
				TL[i][j]       = 'd';		
			}
		}
		else{
			result.matches = mq;
			result.gaps    = gq;	
			TL[i][j]       = 'd';					
		}
	}

	return result;
}

VMG scalar_max2( VMG p1, VMG p2, int ws, int wg ){
	VMG result;

	/* auxiliar variables */
	int g1, g2;
	int m1, m2;
	int r1, r2;		/* temporary results */

	g1 = p1.gaps + IND_PENALTY;		/* add additional indel */
	g2 = p2.gaps + 1;				/* open new gap */

	m1 = p1.matches;
	m2 = p2.matches;

	r1 = g1 == INT_MAX ? INT_MIN : ws * m1 - wg * g1;
	r2 = ws * m2 - wg * g2;                      

	if( r1 >= r2 ){
		result.matches = m1;
		result.gaps    = g1;
	}
	else{
		result.matches = m2;
		result.gaps    = g2;		
	}

	return result;
}

VMG scalar_max3( VMG S, VMG T, VMG Q, int match, int ws, int wg ){
	VMG result;

	/* auxiliar variables */
	int gs, gt, gq;
	int ms, mt, mq;
	int rs, rt, rq;		/* temporary results */

	gs = S.gaps;
	gt = T.gaps;
	gq = Q.gaps;

	ms = S.matches;
	mt = T.matches;
	mq = Q.matches + match;

	rs = ws * ms - wg * gs;
	rt = ws * mt - wg * gt;
	rq = ws * mq - wg * gq;

	if( rs >= rt ){
		if( rs >= rq ){
			result.matches = ms;
			result.gaps    = gs;
		}
		else{
			result.matches = mq;
			result.gaps    = gq;
		}
	}
	else{
		if( rt >= rq ){
			result.matches = mt;
			result.gaps    = gt;
		}
		else{
			result.matches = mq;
			result.gaps    = gq;
		}
	}

	return result;
}

void prune( Pareto_set * set, int i, int j, char dir ){
	/**
	 * If the current score vector plus the corresponding upper bound
	 * score vector is dominated by any of the lower bound scores, than
	 * it is pruned.
	 */
	int add_gap, m, g, h, l, k;
	
	h = M-i;			
	l = N-j;

	/* checks if a gap already exists */
	if( h == l )		add_gap = 0;
	else if( h > l )	add_gap = dir == 'u' ? 0 : 1;
	else				add_gap = dir == 'l' ? 0 : 1;
	
	/* auxiliar variables 'm', 'g' */
	m = set->scores[ set->num-1 ].matches + U[i][j];	
	g = set->scores[ set->num-1 ].gaps + add_gap;


	k = 0;
    while( k < BOUNDS_NUM && LB[k].matches < m ) k++;	/* advance in bound scores */

    if( k < BOUNDS_NUM ){
        if( m == LB[k].matches ){
            if( g > LB[k].gaps ){
                set->num--;  
            }      
        }
        else{
            if( g >= LB[k].gaps ){
                set->num--;       
            }     
        }
    }
    else{	/* advanced all bounds (this should never happen) */
        set->num--;
    }
}

int max2(int n1, int n2){
	return n1 >= n2 ? n1 : n2;
}
int min2(int n1, int n2){
	return n1 <= n2 ? n1 : n2;
}

void remove_bounds_tables(){
	int i;

	free(LB);

	for( i = 0; i < M+1; ++i ){
		free( U[i] );
	}
	free(U);
}
