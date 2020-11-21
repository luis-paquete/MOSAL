/****************************************************************************

 It computes the set of non dominated scores for the "Substitution score"-#Gaps
 problem.
 It is a Dynamic Programming (DP) with pruning approach, by keeping three
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

/**
 * Declaration of variables and function prototypes
 * related with the establishment and initialization
 * of all bounds: lower (MIN, MID's and MAX) and upper bounds.
 * MIN is the score vector with minimum #gaps and minimum substitution score
 * MAX is the score vector with maximum #gaps and maximum substitution score
 * MID's are all the score vectors between.
 * It also includes the function responsible for pruning.
 */

#ifndef _BOUNDS_H
#define _BOUNDS_H

#include "pareto_set.h"

/* 	length of the alphabet to the substitution score table
	(includes all the letters from '*' to 'Z' to get constant time access)
*/
#define LEN_ALPHABET 	('Z'-'*'+1)

VMG *LB;			/* MIN, MIDs, MAX together (ordered by substitution score) */
VMG *MIN, *MAX;		/* MIN and MAX are calculated in a different way */

VMG ** L;			/* Dynamic Programming (DP) table for the main lower bound states table (m,g) */
VMG ** LT;			/* DP table for the lower bound states table (m,g) that keeps solutions ending in (ai,'-') */
VMG * LS;			/* DP table for the lower bound states table (m,g) that keeps solutions ending in ('-',bj) */

char ** TL;  		/* lower bound traceback directions (needed to know if a previous gap already exists to improve puning) */

/* MID bound states tables */
VMG ** M_;		/* main DP table */
VMG * MS;		/* DP table S */
VMG ** MT;		/* DP table T */

int ** U;	/* DP table for the upper bound states table (only substitution score, #gaps are calculated in another simpler way) */

int MAX_MID_BOUNDS_NUM;				/* number of mid bounds */
int BOUNDS_NUM;					/* total number of bounds */

extern int M, N;
extern char seq1[];
extern char seq2[];
extern int SS[][LEN_ALPHABET];

/* bounds initialization and determination (calls all init_*_bound() and calculate_*_bound() functions )*/
void init_bounds( int mid_bounds );
void init_lower_bound( );			/* initialization of the lower bound structures (memory allocations and base cases specification) */
void calculate_lower_bound( );		/* computation of all lower bounds (MIN, MAX and MID's) by filling 'VMG ** L' and 'VMG ** M_' */
void init_upper_bound( );			/* initialization of the upper bound structures (memory allocations and base cases specification) */
void calculate_upper_bound( );		/* computation of all upper bounds by filling 'int ** U' */

/**
 * Returns the lexicographically larger vector between 
 * LS and L current score vectors
 * This is used to determine MAX, together with lexmaxM3() - 'matches' are priority 
 */
VMG lexmaxMS( VMG p1, VMG p2, int i, int j );
/**
 * Returns the lexicographically larger vector between 
 * LT and L current score vectors
 * This is used to determine MAX, together with lexmaxM3() - 'matches' are priority 
 */
VMG lexmaxMT( VMG p1, VMG p2, int i, int j );

/**
 * Returns the lexicographically larger vector between the score
 * vectors of cells 'S', 'T' and 'Q'
 * 'score'- the value of the subs score corresponding the current letters at position (i,j) (used by the 'Q' score vector).
 * This is used to determine MAX - 'matches' are priority
 */
VMG lexmaxM3( VMG S, VMG T, VMG Q, int score, int i, int j);
void calc_MAX();					/* calculates MAX, based on lexmaxM() */

void calc_MIN();					/* calculates MIN */

/**
 * Returns the score vector that maximizes the scalarized score function
 * between the MS, MT and M_('Q') current score vectors.
 * (another version of the bicriteria alignment problem).
 * 'subs_score'- the value of the subs score corresponding the current letters at position (i,j) (used by the 'Q' score vector).
 * 'ws' and 'wd' are the weighting coefficients.
 */
VMG scalar_max3( VMG S, VMG T, VMG Q, int subs_score, int ws, int wg );
/**
 * Returns the score vector that maximizes the scalarized score function
 * between the MS and M_('Q') or MT and M_('Q') current score vectors.
 * (another version of the bicriteria alignment problem).
 * 'ws' and 'wd' are the weighting coefficients.
 */
VMG scalar_max2( VMG p1, VMG p2, int ws, int wg );
void calc_MID();					/* calculates all MID bounds for each MID bound DP table */

/**
 * Determines if the last score vector of the 
 * pareto set 'set' can be pruned and, if so,
 * removes it from the set.
 * Based on the lower and upper bounds previously calculated.
 * The 'dir' improves pruning by knowing if a previous gap already exists
 */
void prune( Node * previous_state, int i, int j, char dir );

void remove_bounds_tables();		/* frees tables still in memory and used during the algorithm */

int get_score( int i, int j );
int get_score2(char *s1, int i,  char *s2, int j );

int max2(int n1, int n2);				/* returns the maximum of two numbers */
int max3(int n1, int n2, int n3);		/* returns the maximum of three numbers */


#endif