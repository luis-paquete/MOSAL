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
 * Definition of the structures that compose
 * each Dynamic Programming (DP) table cell and the 
 * functions that fill the table.
 */

#ifndef _PARETO_SET_H
#define _PARETO_SET_H

/* score vector */
typedef struct score_vector{
	int matches;
	int indels;
}VMD;

/* Set of pareto optimal alignment scores.
   Corresponds to a DP table cell.
 */
typedef struct set{
	VMD * scores;	
	int num;		/* number of scores */
}Pareto_set;

extern int MAX_STATES;

/**
 * Stores a new score to the array of scores in set->scores.
 * m - #matches
 * i - #indels
 */
void append_node(Pareto_set *set, int m, int i );
/**
 * Calculates the set of pareto optimal alignments alreadt pruned
 * between 'up', 'diagonal' and 'left' sets.
 * Called in each DP algorithm iteration.
 * 'subs_score' - the value of the subs score corresponding the current letters at position (i,j) (used by the diagonal set).
 *
 * Stores the result in the set pointed by 'ptr_res'
 */
void pareto_merge( Pareto_set *ptr_res, Pareto_set up, Pareto_set diagonal, Pareto_set left, int subs_score, int i, int j );

#endif