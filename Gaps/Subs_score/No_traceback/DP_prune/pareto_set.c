/****************************************************************************

 It computes the set of non dominated scores for the "Substitution score"-#Gaps problem
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
 * Implementation of the functions in 'pareto_set.h'.
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "pareto_set.h"
#include "bounds.h"


void append_node(Pareto_set *set, int m, int g ){
	if( set->num < MAX_STATES ){
		set->scores[ set->num ].matches = m;
		set->scores[ set->num ].gaps = g;
		set->num++;
	}
	else{
		printf("Error on inserting node! Exceeded MAX_STATES=%d\n", MAX_STATES);
		exit(-1);
	}
}


void pareto_merge2( Pareto_set *ptr_res, Pareto_set p1, Pareto_set p2, int i, int j, char dir ){	
	/**
	 * Scores are assumed to be ordered incrementally by their #gaps.
	 * Until all the iterators reach the end, it is selected the state with
	 * minimum #gaps and, in case of equality, with maximum substitution score, from
	 * the two states: p1-(ss_p1,gap_p1), p2-(ss_p2,gap_p2).
	 * This state is appended to the result as it is always non dominated.
	 * The substitution score of that state (to be stored in SS_MAX) is also used to 
	 * eliminate the (dominated) states following in each set, by advancing 
	 * the iterators until the substitution score of the pointed state is greater than SS_MAX.
	 * The dominance of these states is also assured because the sets are ordered by #gaps
	 * and they all have greater #gaps compared to the selected state.
	 */
	
	/* defines the substitution score of the state to be inserted */
	int SS_MAX;
	/* auxiliar variables for comparisons */
	int gap_p1, gap_p2;
	int ss_p1, ss_p2;

	int it_p1, it_p2;		/* iterators to each set ot pareto optimal alignments */
	it_p1 = it_p2 = 0;

	while( (it_p1 < p1.num) || (it_p2 < p2.num) ){

		/*  Initialization of the auxiliar variables for the p1 and p2 cells.
			In case iterators reach the end of its set (first comparison),
			its values are not inserted, keeping the algorithm structure */

		/* P1 */
		if( it_p1 < p1.num ){
			gap_p1 = p1.scores[ it_p1 ].gaps + IND_PENALTY;		/* additional indel penalty */
			ss_p1 = p1.scores[ it_p1 ].matches;
		}
		else{
			gap_p1 = INT_MAX;
			ss_p1 = INT_MIN;			
		}

		/* P2 */
		if( it_p2 < p2.num ){
			gap_p2 = p2.scores[ it_p2 ].gaps + 1;				/* open new gap */
			ss_p2 = p2.scores[ it_p2 ].matches;
		}
		else{
			gap_p2 = INT_MAX;
			ss_p2 = INT_MIN;			
		}		
		/* End of auxiliar variables initialization */


		/*  Select score with minimum #gaps and, in case of equality,
			with maximum substitution score */
		if( gap_p1 < gap_p2 || (gap_p1==gap_p2 && ss_p1 >= ss_p2) ){
			SS_MAX = ss_p1;
			append_node( ptr_res, ss_p1, gap_p1);
		}
		else{
			SS_MAX = ss_p2;
			append_node( ptr_res, ss_p2, gap_p2);	
		}
		/* End of minimum #gaps score selection */

		prune(ptr_res, i, j, dir);

		/* Exclude dominated scores, advancing iterators, based on the value SS_MAX */
		while( it_p1 < p1.num && p1.scores[ it_p1 ].matches <= SS_MAX )	it_p1++;
		while( it_p2 < p2.num && p2.scores[ it_p2 ].matches <= SS_MAX )	it_p2++;

	}
}

void pareto_merge3( Pareto_set *ptr_res, Pareto_set S, Pareto_set T, Pareto_set Q, int subs_score, int i, int j ){
	/**
	 * Scores are assumed to be ordered incrementally by their #gaps.
	 * Until all the iterators reach the end, it is selected the state with
	 * minimum #gaps and, in case of equality, with maximum substitution score, from
	 * the three states: S-(ms,gs), Q-(mq,gq) and T-(mt,gt).
	 * This state is appended to the result as it is always non dominated.
	 * The substitution score of that state (to be stored in SS_MAX) is also used to 
	 * eliminate the (dominated) states following in each set, by advancing 
	 * the iterators until the substitution score of the pointed state is greater than SS_MAX.
	 * The dominance of these states is also assured because the sets are ordered by #gaps
	 * and they all have greater #gaps compared to the selected state.
	 */
	
	char dir; 		/* in order to improve the pruning */

	/* defines the substitution score of the state to be inserted */
	int SS_MAX;
	/* auxiliar variables for comparisons */
	int gt, gq, gs;
	int mt, mq, ms;

	int it_s, it_t, it_q;			/* iterators to each set ot pareto optimal alignments */
	it_s = it_t = it_q = 0;

	while(  (it_t < T.num) || (it_q < Q.num) || (it_s < S.num)  ){				

		/*  Initialization of the auxiliar variables for the Q, S and T cells.
			In case iterators reach the end of its set (first comparison),
			its values are not inserted, keeping the algorithm structure */
		
		/* Q */
		if( it_q < Q.num ){
			gq = Q.scores[ it_q ].gaps;
			mq = Q.scores[ it_q ].matches + subs_score;
		}
		else{
			gq = INT_MAX;
			mq = INT_MIN;
		}
		/* S */
		if( it_s < S.num){
			gs = S.scores[ it_s ].gaps;
			ms = S.scores[ it_s ].matches;
		}
		else{
			gs = INT_MAX;
			ms = INT_MIN;
		}
		/* T */
		if( it_t < T.num ){
			gt = T.scores[ it_t ].gaps;
			mt = T.scores[ it_t ].matches;
		}
		else{
			gt = INT_MAX;
			mt = INT_MIN;
		}
		/* End of auxiliar variables initialization */


		/*  Select score with minimum #gaps and, in case of equality,
			with maximum substitution score */
		if( gt < gq ){

			if( gt < gs ){
				SS_MAX = mt;						
				append_node( ptr_res, mt, gt );		
				dir = 'u';						
			}
			else if( gt == gs && mt >= ms ){
				SS_MAX = mt;		
				append_node( ptr_res, mt, gt );	
				dir = 'u';							
			}
			else{
				SS_MAX = ms;		
				append_node( ptr_res, ms, gs );		
				dir = 'l';										
			}
		}

		else if( gt == gq ){

			if( gs < gt ){
				SS_MAX = ms;		
				append_node( ptr_res, ms, gs );
				dir = 'l';				
			}
			else if( gt == gs ){
				if( mt >= mq ){
					if( mt >= ms ){
						SS_MAX = mt;		
						append_node( ptr_res, mt, gt );	
						dir = 'u';
					}
					else{
						SS_MAX = ms;		
						append_node( ptr_res, ms, gs );
						dir = 'l';		
					}
				}
				else{
					if( mq >= ms ){
						SS_MAX = mq;		
						append_node( ptr_res, mq, gq );
						dir = 'd';
					}
					else{
						SS_MAX = ms;		
						append_node( ptr_res, ms, gs );
						dir = 'l';	
					}
				}
			}
			else{			/* gt==gq && gt < gs */
				if( mt >= mq ){
					SS_MAX = mt;		
					append_node( ptr_res, mt, gt );
					dir = 'u';
				}
				else{
					SS_MAX = mq;		
					append_node( ptr_res, mq, gq );	
					dir = 'd';
				}
			}
		}

		else{

			if( gq < gs ){
				SS_MAX = mq;		
				append_node( ptr_res, mq, gq );	
				dir = 'd';				
			}
			else if( gq == gs && mq >= ms ){
				SS_MAX = mq;		
				append_node( ptr_res, mq, gq );	
				dir = 'd';				
			}
			else{
				SS_MAX = ms;		
				append_node( ptr_res, ms, gs );
				dir = 'l';								
			}
		}
		/* End of minimum #gaps score selection */

		/* prune last score appended */
		prune(ptr_res, i, j, dir);

		/* Exclude dominated scores, advancing iterators, based on the value SS_MAX */
		while( it_s < S.num && S.scores[ it_s ].matches <= SS_MAX ) it_s++;
		while( it_t < T.num && T.scores[ it_t ].matches <= SS_MAX ) it_t++;
		while( it_q < Q.num && Q.scores[ it_q ].matches + subs_score <= SS_MAX ) it_q++;

	}

}
