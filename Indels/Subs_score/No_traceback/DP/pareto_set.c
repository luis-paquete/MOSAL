/****************************************************************************

 It computes the set of non dominated scores for the "Substitution score"-#Indels
 problem without traceback (only final scores). It is a Dynamic Programming
 (DP) without pruning approach, by filling a DP matrix - P
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
 * Implementation of the functions in 'pareto_set.h'.
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "pareto_set.h"

void append_node(Pareto_set *set, int m, int i ){
	if( set->num < MAX_STATES ){
		set->scores[ set->num ].matches = m;
		set->scores[ set->num ].indels = i;
		set->num++;
	}
	else{
		printf("Error on inserting node! Exceeded MAX_STATES=%d\n", MAX_STATES);
		exit(-1);
	}
}

void pareto_merge( Pareto_set *ptr_res, Pareto_set up, Pareto_set diagonal, Pareto_set left, int subs_score ){	
	/**
	 * Scores are assumed to be ordered incrementally by their #indels.
	 * Until all the iterators reach the end, it is selected the state with
	 * minimum #indels and, in case of equality, with maximum substitution score, from
	 * the three states: up-(mu,iu), diagonal-(md,id) and left-(ml,il).
	 * This state is appended to the result as it is always non dominated.
	 * The substitution score of that state (to be stored in SS_MAX) is also used to 
	 * eliminate the (dominated) states following in each set, by advancing 
	 * the iterators until the substitution score of the pointed state is greater than SS_MAX.
	 * The dominance of these states is also assured because the sets are ordered by #indels
	 * and they all have greater #indels compared to the selected state.
	 */	

	/* defines substitution score of the state to be inserted */
	int SS_MAX;
	/* auxiliar variables for comparisons */
	int iu, id, il;
	int mu, md, ml;

	int it_up, it_diagonal, it_left;		/* iterators to each set ot pareto optimal alignments */
	it_up = it_diagonal = it_left = 0;

	while(  (it_up < up.num) || (it_diagonal < diagonal.num) || (it_left < left.num)  ){		

		/*  Initialization of the auxiliar variables for the up, diagonal and left directions
			In case iterators reach the end of its set (first comparison),
			its values are not inserted, keeping the algorithm structure */

		/* UP */
		if ( it_up < up.num ){
			iu = up.scores[ it_up ].indels + 1;
			mu = up.scores[ it_up ].matches;
		}
		else{
			iu = INT_MAX;
			mu = INT_MIN;
		}
		/* DIAGONAL */
		if( it_diagonal < diagonal.num ){
			id = diagonal.scores[ it_diagonal ].indels;
			md = diagonal.scores[ it_diagonal ].matches + subs_score;
		}
		else{
			id = INT_MAX;
			md = INT_MIN;
		}
		/* LEFT */
		if( it_left < left.num ){
			il = left.scores[ it_left ].indels + 1;
			ml = left.scores[ it_left ].matches;
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
				append_node( ptr_res, mu, iu );		
			}
			else if( iu == il && mu >= ml ){
				SS_MAX = mu;		
				append_node( ptr_res, mu, iu );					
			}
			else{
				SS_MAX = ml;		
				append_node( ptr_res, ml, il );									
			}
		}

		else if( iu == id ){

			if( il < iu ){
				SS_MAX = ml;		
				append_node( ptr_res, ml, il );	
			}
			else if( iu == il ){
				if( mu >= md ){
					if( mu >= ml ){
						SS_MAX = mu;		
						append_node( ptr_res, mu, iu );	
					}
					else{
						SS_MAX = ml;		
						append_node( ptr_res, ml, il );		
					}
				}
				else{
					if( md >= ml ){
						SS_MAX = md;		
						append_node( ptr_res, md, id );	
					}
					else{
						SS_MAX = ml;		
						append_node( ptr_res, ml, il );	
					}
				}
			}
			else{			/* iu==id && iu < il */
				if( mu >= md ){
					SS_MAX = mu;		
					append_node( ptr_res, mu, iu );	
				}
				else{
					SS_MAX = md;		
					append_node( ptr_res, md, id );	
				}
			}
		}

		else{

			if( id < il ){
				SS_MAX = md;		
				append_node( ptr_res, md, id );					
			}
			else if( id == il && md >= ml ){
				SS_MAX = md;		
				append_node( ptr_res, md, id );					
			}
			else{
				SS_MAX = ml;		
				append_node( ptr_res, ml, il );									
			}
		}
		/* End of minimum #indels score selection */

		/* Exclude dominated scores, advancing iterators, based on the value SS_MAX */
		while( it_up < up.num && up.scores[ it_up ].matches <= SS_MAX ) it_up++;
		while( it_left < left.num && left.scores[ it_left ].matches <= SS_MAX ) it_left++;
		while( it_diagonal < diagonal.num && diagonal.scores[ it_diagonal ].matches + subs_score <= SS_MAX ) it_diagonal++;

	}

}
