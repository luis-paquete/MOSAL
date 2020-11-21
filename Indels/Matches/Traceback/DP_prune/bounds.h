/****************************************************************************

 It computes the set of non dominated scores for the #Matches-#Indels problem.
 It is a Dynamic Programming (DP) with pruning approach, by filling a DP
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
                 
/**
 * Declaration of variables and function prototypes
 * related with the establishment and initialization
 * of all bounds: lower (MIN, MID's and MAX) and upper bounds.
 * MIN is the score vector with minimum #indels and minimum #matches
 * MAX is the score vector with maximum #indels and maximum #matches
 * MID's are all the score vectors between.
 * It also includes the function responsible for pruning.
 */

#ifndef _BOUNDS_H
#define _BOUNDS_H

#include "pareto_set.h"

VMD *LB;			/* MIN, MIDs, MAX together (ordered by #matches) */
VMD *MIN, *MAX;		/* MIN and MAX are calculated in a different way */

VMD ** L;			/* Dynamic Programming (DP) table for the lower bound states table (m,d) */
VMD ** U;			/* DP table for the upper bound states table (m,d) */
VMD ** M_; 			/* MID bound states DP tables */

int MAX_MID_BOUNDS_NUM;			/* number of mid bounds */
int BOUNDS_NUM;					/* total number of bounds */

extern int M, N;
extern char seq1[];
extern char seq2[];

/* bounds initialization and determination (calls all init_*_bound() and calculate_*_bound() functions )*/
void init_bounds( int mid_bounds );
void init_lower_bound( );			/* initialization of the lower bound structures (memory allocations and base cases specification) */
void calculate_lower_bound( );		/* computation of all lower bounds (MIN, MAX and MID's) by filling 'VMD ** L' */
void init_upper_bound( );			/* initialization of the upper bound structures (memory allocations and base cases specification) */
void calculate_upper_bound( );		/* computation of all upper bounds by filling 'VMD ** U' */

/**
 * Returns the lexicographically larger vector between the score
 * vectors of cells 'up', 'left' and 'diagonal' of the current cell
 * in the DP table 'VMD ** L'.
 * 'match' - if the current sequence letters match (used by the diagonal score vector)
 * This is used to determine MAX - 'matches' are priority
 */
VMD lexmaxM( VMD up, VMD left, VMD diagonal, int match );
void calc_MAX();					/* calculates MAX, based on lexmaxM() */
/**
 * Returns the lexicographically smaller vector between the score
 * vectors of cells 'up', 'left' and 'diagonal' of the current cell
 * in the DP table 'VMD ** L'.
 * 'match' - if the current sequence letters match (used by the diagonal score vector)
 * This is used to determine MIN - 'indels' are priority
 */
VMD lexminD( VMD up, VMD left, VMD diagonal, int match );
void calc_MIN();					/* calculates MIN, based on lexminD() */



/**
 * Returns the score vector that maximizes the scalarized score function
 * (another version of the bicriteria alignment problem).
 * 'match' - if the current sequence letters match (used by the diagonal score vector)
 * 'ws' and 'wd' are the weighting coefficients.
 */
VMD scalar_max( VMD up, VMD left, VMD diagonal, int match, int ws, int wd );
void calc_MID();					/* calculates all MID bounds for each MID bound DP table, based on scalar_max() */

/**
 * Determines if the score vector following 'previous_state'
 * can be pruned and, if so, removes it from the linked list.
 * Based on the lower and upper bounds previously calculated.
 */
void prune( Node * previous_state, int i, int j);

void remove_bounds_tables();		/* frees tables still in memory and used during the algorithm */

int max2(int n1, int n2);			/* returns the maximum of two numbers */

#endif