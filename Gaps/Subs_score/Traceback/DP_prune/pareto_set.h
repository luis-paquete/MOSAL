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
typedef struct n{
	VMG score;						/* score vector */
	struct n* traceback_ptr;		/* traceback pointer */
	char traceback_dir; 			/* 'T'-up   'S'-left   'Q'-diagonal */
	struct n* next;					/* pointer to next score */
}Node;

typedef Node * Pareto_set;

/**
 * Returns a new linked list to store pareto optimal alignments.
 * The list has a header.
 */
Pareto_set create_list( );

/**
 * Appends a new score vector to the linked list pointed by 'ptr'
 * m - #matches
 * g - #gaps
 */
void append_node(Node *ptr, int m, int g, Node *tb_ptr, char tb_dir);

/**
 * Calculates the set of pareto optimal alignments
 * between 'S', 'T' and 'Q' sets.
 * Called in each DP algorithm iteration.
 * Also, it passes the position in order to prune or not the result
 * to be inserted.
 * 
 * Returns a linked list with the result of non dominated score vectors.
 */
Pareto_set pareto_merge2( Pareto_set p1, Pareto_set p2, char tb_dir );

/**
 * Calculates the set of pareto optimal alignments
 * between 'S', 'T' and 'Q' sets.
 * Called in each DP algorithm iteration.
 * 'subs_score' - the value of the subs score corresponding the current letters
 * 				  at position (i,j) (used by the 'Q' set).
 * Also, it passes the position in order to prune or not the result
 * to be inserted.
 * 
 * Returns a linked list with the result of non dominated score vectors.
 */
Pareto_set pareto_merge3( Pareto_set S, Pareto_set T, Pareto_set Q, int subs_score, int i, int j );

/* frees the linked list pointed to by 'ptr' */
void remove_list( Node * ptr );

#endif