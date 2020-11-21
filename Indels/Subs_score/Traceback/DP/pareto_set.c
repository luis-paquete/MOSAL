/****************************************************************************

 It computes the set of non dominated scores for the "Substitution score"-#Indels
 problem.
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
 * Implementation of the functions in 'pareto_set.h'.
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "pareto_set.h"

Pareto_set create_list( ){
	Node *ptr = (Node *)malloc(sizeof(Node));
	if(ptr){
		ptr->score.matches = 0;
		ptr->score.indels = 0;
		ptr->traceback_ptr = NULL;
		ptr->traceback_ptr = '\0';
		ptr->next = NULL;
	}

	return ptr;
}

void append_node(Node *ptr, int m, int i, Node *tb_ptr, char tb_dir){
	ptr->next = (Node *)malloc(sizeof(Node));
	ptr->next->score.matches = m;
	ptr->next->score.indels = i;
	ptr->next->traceback_ptr = tb_ptr;
	ptr->next->traceback_dir = tb_dir;
	ptr->next->next = NULL;
}

Pareto_set pareto_merge( Pareto_set up, Pareto_set diagonal, Pareto_set left, int subs_score ){
	/**
	 * Scores are assumed to be ordered incrementally by their #indels.
	 * Until all the iterators reach the end, it is selected the state with
	 * minimum #indels and, in case of equality, with maximum substitution score, from
	 * the three states: up-(mu,iu), diagonal-(md,id) and left-(ml,il).
	 * (If the score vectors are exactly the same, than one is chosen randomly.)
	 * This state is appended to the result as it is always non dominated.
	 * The substitution score of that state (to be stored in SS_MAX) is also used to 
	 * eliminate the (dominated) states following in each set, by advancing 
	 * the iterators until the substitution score of the pointed state is greater than SS_MAX.
	 * The dominance of these states is also assured because the sets are ordered by #indels
	 * and they all have greater #indels compared to the selected state.
	 *
	 */	

	/* allocate memory for the result set */
	Pareto_set res = create_list();
	Node * ptr_res = res;

	/* defines the substitution score of the state to be inserted */
	int SS_MAX;
	/* auxiliar variables for comparisons */
	int iu, id, il;
	int mu, md, ml;

	/* these linked lists have headers */
	up       = up->next;
	diagonal = diagonal->next;
	left     = left->next;

	while(  up || diagonal || left  ){			
		
		/*  Initialization of the auxiliar variables for the up, diagonal and left directions
			In case iterators reach the end of its set (first comparison),
			its values are not inserted, keeping the algorithm structure */

		/* UP */
		if ( up ){
			iu = up->score.indels + 1;
			mu = up->score.matches;
		}
		else{
			iu = INT_MAX;
			mu = INT_MIN;
		}
		/* DIAGONAL */
		if( diagonal ){
			id = diagonal->score.indels;
			md = diagonal->score.matches + subs_score;
		}
		else{
			id = INT_MAX;
			md = INT_MIN;
		}
		/* LEFT */
		if( left ){
			il = left->score.indels + 1;
			ml = left->score.matches;
		}
		else{
			il = INT_MAX;
			ml = INT_MIN;
		}
		/* End of auxiliar variables initialization */	

				
		/*  Select score with minimum #indels and, in case of equality,
			with maximum substitution score */
		if( iu < id ){

			if( iu < il ){
				SS_MAX = mu;		
				append_node( ptr_res, mu, iu , up, 'u' );					
			}
			else if( iu == il ){
				if( mu > ml ){
					SS_MAX = mu;		
					append_node( ptr_res, mu, iu , up, 'u' );										
				}
				else if( mu == ml && (rand()%2) ){
					SS_MAX = mu;		
					append_node( ptr_res, mu, iu , up, 'u' );
				}
				else{
					SS_MAX = ml;		
					append_node( ptr_res, ml, il , left, 'l' );														
				}
			}
			else{
				SS_MAX = ml;		
				append_node( ptr_res, ml, il , left, 'l' );									
			}
		}

		else if( iu == id ){

			if( il < iu ){
				SS_MAX = ml;		
				append_node( ptr_res, ml, il , left, 'l' );	
			}
			else if( iu == il ){
				if( mu > md ){
					if( mu > ml ){
						SS_MAX = mu;		
						append_node( ptr_res, mu, iu , up, 'u' );	
					}
					else if( mu == ml && (rand()%2) ){
						SS_MAX = mu;		
						append_node( ptr_res, mu, iu , up, 'u' );							
					}
					else{
						SS_MAX = ml;		
						append_node( ptr_res, ml, il , left, 'l' );		
					}
				}
				else if( mu == md ){
					if( ml > mu ){
						SS_MAX = ml;		
						append_node( ptr_res, ml, il , left, 'l' );							
					}
					else if( mu == ml ){
						int var = rand()%3;
						if( var == 0){
							SS_MAX = ml;		
							append_node( ptr_res, ml, il , left, 'l' );														
						}
						else if (var == 1){
							SS_MAX = mu;		
							append_node( ptr_res, mu, iu , up, 'u' );							
						}
						else{
							SS_MAX = md;		
							append_node( ptr_res, md, id , diagonal, 'd' );	
						}
					}
					else{
						if( rand()%2 ){
							SS_MAX = mu;		
							append_node( ptr_res, mu, iu , up, 'u' );							
						}
						else{							
							SS_MAX = md;		
							append_node( ptr_res, md, id , diagonal, 'd' );	
						}
					}
				}
				else{		/* iu == il == id && md > mu */
					if( md > ml ){
						SS_MAX = md;		
						append_node( ptr_res, md, id , diagonal, 'd' );	
					}
					else if( md == ml && (rand()%2) ){
						SS_MAX = md;		
						append_node( ptr_res, md, id , diagonal, 'd' );							
					}
					else{
						SS_MAX = ml;		
						append_node( ptr_res, ml, il , left, 'l' );	
					}
				}
			}
			else{			/*  iu==id && iu < il */
				if( mu > md ){
					SS_MAX = mu;		
					append_node( ptr_res, mu, iu , up, 'u' );	
				}
				else if( mu == md && (rand()%2) ){
					SS_MAX = mu;		
					append_node( ptr_res, mu, iu , up, 'u' );	
				}
				else{
					SS_MAX = md;		
					append_node( ptr_res, md, id , diagonal, 'd' );	
				}
			}
		}

		else{

			if( id < il ){
				SS_MAX = md;		
				append_node( ptr_res, md, id , diagonal, 'd' );					
			}
			else if( id == il ){
				if( md > ml ){
					SS_MAX = md;		
					append_node( ptr_res, md, id , diagonal, 'd' );										
				}
				else if( md == ml && (rand()%2) ){
					SS_MAX = md;		
					append_node( ptr_res, md, id , diagonal, 'd' );										
				}
				else{					
					SS_MAX = ml;		
					append_node( ptr_res, ml, il , left, 'l' );														
				}
			}
			else{
				SS_MAX = ml;		
				append_node( ptr_res, ml, il , left, 'l' );									
			}
		}
		/* End of minimum #indels score selection */

		ptr_res = ptr_res->next;

		/* Exclude dominated scores, advancing iterators, based on the value SS_MAX */
		while( up && up->score.matches <= SS_MAX ) up = up->next;
		while( left && left->score.matches <= SS_MAX ) left = left->next;
		while( diagonal && diagonal->score.matches + subs_score <= SS_MAX ) diagonal = diagonal->next;

	}

	return res;

}

void remove_list( Node * ptr ){
	if(ptr){
		Node * tmp;
		while( ptr->next ){		
			tmp = ptr;
			ptr = ptr->next;
			free(tmp);
		}
		free(ptr);
	}
}
