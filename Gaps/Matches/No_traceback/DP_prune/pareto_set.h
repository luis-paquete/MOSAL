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
 * Definition of the structures that compose
 * each Dynamic Programming (DP) table cell and the 
 * functions that fill the DP tables.
 */

#ifndef _PARETO_SET_H
#define _PARETO_SET_H

#define IND_PENALTY			(0)		/* value for each additional indel */

/* score vector */
typedef struct score_vector{
	int matches;
	int gaps;
}VMG;

/* Set of pareto optimal alignment scores.
   Corresponds to a DP table cell.
 */
typedef struct set{
	VMG * scores;
	int num;		/* number of scores */
}Pareto_set;

extern int MAX_STATES;

/**
 * Stores a new score to the array of scores in set->scores.
 * m - #matches
 * g - #gaps
 */
void append_node(Pareto_set *set, int m, int g );

/**
 * Calculates the set of pareto optimal alignments
 * between 'S' and 'Q' or 'T' and 'Q' sets.
 * Called in each DP algorithm iteration.
 * Also, it passes the position on the DP table and also the traceback
 * direction (S or T) in order to prune or not the result to be inserted.
 *
 * Stores the result in the set pointed by 'ptr_res'.
 * 
 */
void pareto_merge2( Pareto_set *ptr_res, Pareto_set p1, Pareto_set p2 , int i, int j, char dir);

/**
 * Calculates the set of pareto optimal alignments
 * between 'Q', 'S' and 'T' sets.
 * Called in each DP algorithm iteration.
 * 'match' - if the current sequence letters match (used by the 'Q' set).
 * Also, it passes the position in order to prune or not the result
 * to be inserted.
 * 
 * Stores the result in the set pointed by 'ptr_res'
 */
void pareto_merge3( Pareto_set *ptr_res, Pareto_set S, Pareto_set T, Pareto_set Q, int match, int i, int j );

#endif