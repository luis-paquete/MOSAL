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
                 
/**
 * Definition of the structures that compose
 * each Dynamic Programming (DP) table cell and the 
 * functions that fill the table.
 */

/* score vector */
typedef struct score_vector{
	int matches;
	int indels;
}VMD;

typedef struct n{
	VMD score;						/* score vector */
	struct n* traceback_ptr;		/* traceback pointer */
	char traceback_dir; 			/* 'u'-up   'l'-left   'd'-diagonal */
	struct n* next;					/* pointer to next score */
}Node;

typedef Node * Pareto_set;			/* linked list of nodes 'Node' */

/**
 * Returns a new linked list to store pareto optimal alignments.
 * The list has a header.
 */
Node * create_list( );

/**
 * Appends a new score vector to the linked list pointed by 'ptr'
 * m - #matches
 * i - #indels
 */
void append_node(Node *ptr, int m, int i, Node *tb_ptr, char tb_dir);

/**
 * Calculates the set of pareto optimal alignments
 * between 'up', 'diagonal' and 'left' sets.
 * Called in each DP algorithm iteration.
 * 'match' - if the current sequence letters match (used by the diagonal set).
 *
 * Returns a linked list with the result of non dominated score vectors.
 */
Pareto_set pareto_merge( Pareto_set up, Pareto_set diagonal, Pareto_set left, int match );

/* frees the linked list pointed to by 'ptr' */
void remove_list( Node * ptr );