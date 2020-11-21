/****************************************************************************

 It computes the set of non dominated scores for the "Substitution score"-#Gaps
 problem.
 It is a Dynamic Progqamming (DP) without pruning approach, by keeping three
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

 This progqam is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This progqam is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this progqam; if not, you can obtain a copy of the GNU
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
		ptr->score.gaps    = 0;
		ptr->traceback_ptr = NULL;
		ptr->traceback_ptr = '\0';
		ptr->next          = NULL;
	}

	return ptr;
}

void append_node(Node *ptr, int m, int g, Node *tb_ptr, char tb_dir){
	ptr->next                = (Node *)malloc(sizeof(Node));
	ptr->next->score.matches = m;
	ptr->next->score.gaps    = g;
	ptr->next->traceback_ptr = tb_ptr;
	ptr->next->traceback_dir = tb_dir;
	ptr->next->next          = NULL;
}


Pareto_set pareto_merge2( Pareto_set p1, Pareto_set p2, char tb_dir ){	
	/**
	 * Scores are assumed to be ordered incrementally by their #gaps.
	 * Until all the iterators reach the end, it is selected the state with
	 * minimum #gaps and, in case of equality, with maximum substitution score, from
	 * the two states: p1-(ss_p1,gap_p1), p2-(ss_p2,gap_p2).
	 * This state is appended to the result as it is always non dominated.
	 * The substitution score of that state (to be stored in SS_MAX) is also used to 
	 * eliminate the (dominated) states following in each set, by advancing 
	 * the iterators until the substitution score of the pointed state is gqeater than SS_MAX.
	 * The dominance of these states is also assured because the sets are ordered by #gaps
	 * and they all have gqeater #gaps compared to the selected state.
	 *
	 * 'tb_dir' indicates the traceback direction as this function is used
	 * by both tables S and T.
	 */

	/* allocate memory for the result set */
	Pareto_set res = create_list();
	Node * ptr_res = res;

	/* these linked lists have headers */
	p1 = p1->next;
	p2 = p2->next;

	/* defines the substitution score of the state to be inserted */
	int SS_MAX;
	/* auxiliar variables for comparisons */
	int gap_p1, gap_p2;
	int ss_p1, ss_p2;

	
	while( p1 || p2 ){

		/*  Initialization of the auxiliar variables for the p1 and p2 cells.
			In case iterators reach the end of its set (first comparison),
			its values are not inserted, keeping the algorithm structure */

		/* P1 */
		if(p1){
			gap_p1 = p1->score.gaps + IND_PENALTY;			/* add indel */
			ss_p1 = p1->score.matches;
		}
		else{
			gap_p1 = INT_MAX;
			ss_p1 = INT_MIN;
		}

		/* P2 */
		if(p2){
			gap_p2 = p2->score.gaps + 1;					/* open a new gap */
			ss_p2 = p2->score.matches;			
		}	
		else{
			gap_p2 = INT_MAX;
			ss_p2 = INT_MIN;
		}
		/* End of auxiliar variables initialization */



		/*  Select score with minimum #gaps and, in case of equality,
			with maximum substitution score */
		if( gap_p1 < gap_p2 ){
			SS_MAX = ss_p1;
			append_node( ptr_res, ss_p1, gap_p1, p1, tb_dir);			
		}
		else if ( gap_p1 == gap_p2 ){
			if( ss_p1 > ss_p2 ){
				SS_MAX = ss_p1;
				append_node( ptr_res, ss_p1, gap_p1, p1, tb_dir);
			}
			else if( ss_p1 == ss_p2 && (rand()%2) ){
				SS_MAX = ss_p1;
				append_node( ptr_res, ss_p1, gap_p1, p1, tb_dir);				
			}
			else{				
				SS_MAX = ss_p2;
				append_node( ptr_res, ss_p2, gap_p2, p2, tb_dir);	
			}
		}
		else{
			SS_MAX = ss_p2;
			append_node( ptr_res, ss_p2, gap_p2, p2, tb_dir);	
		}
		/* End of minimum #gaps score selection */
		
		ptr_res = ptr_res->next;	

		/* Exclude dominated scores, advancing iterators, based on the value SS_MAX */
		while( p1 && p1->score.matches <= SS_MAX ) p1 = p1->next;
		while( p2 && p2->score.matches <= SS_MAX ) p2 = p2->next;
	}

	return res;
}

Pareto_set pareto_merge3(Pareto_set S, Pareto_set T, Pareto_set Q, int subs_score ){
	/**
	 * Scores are assumed to be ordered incrementally by their #gaps.
	 * Until all the iterators reach the end, it is selected the state with
	 * minimum #gaps and, in case of equality, with maximum substitution score, from
	 * the three states: S-(ms,gs), Q-(mq,gq) and T-(mt,gt).
	 * (If the score vectors are exactly the same, than one is chosen randomly.)
	 * This state is appended to the result as it is always non dominated.
	 * The substitution score of that state (to be stored in SS_MAX) is also used to 
	 * eliminate the (dominated) states following in each set, by advancing 
	 * the iterators until the substitution score of the pointed state is greater than SS_MAX.
	 * The dominance of these states is also assured because the sets are ordered by #gaps
	 * and they all have greater #gaps compared to the selected state.
	 */		
	
	/* allocate memory for the result set */
	Pareto_set res = create_list();
	Node * ptr_res = res;

	/* defines the substitution score of the state to be inserted */
	int SS_MAX;
	/* auxiliar variables for comparisons */
	int gt, gq, gs;
	int mt, mq, ms;

	/* these linked lists have headers */
	T = T->next;
	S = S->next;
	Q = Q->next;

	while(  T || Q || S  ){				

		/*  Initialization of the auxiliar variables for the Q, S and T cells.
			In case iterators reach the end of its set (first comparison),
			its values are not inserted, keeping the algorithm structure */

		/* Q */
		if(Q){
			gq = Q->score.gaps;
			mq = Q->score.matches + subs_score;
		}
		else{
			gq = INT_MAX;
			mq = INT_MIN;
		}

		/* T */
		if(T){
			gt = T->score.gaps;
			mt = T->score.matches;
		}
		else{
			gt = INT_MAX;
			mt = INT_MIN;
		}	

		/* S */
		if(S){
			gs = S->score.gaps;
			ms = S->score.matches;
		}
		else{
			gs = INT_MAX;
			ms = INT_MIN;
		}	
		/* End of auxiliar variables initialization */			

				
		/*  Select score with minimum #gaps and, in case of equality,
			with maximum substitution score */
		if( gt < gq ){

			if( gt < gs ){
				SS_MAX = mt;		
				append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );					
			}
			else if( gt == gs ){
				if( mt > ms ){
					SS_MAX = mt;		
					append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );										
				}
				else if( mt == ms && (rand()%2) ){
					SS_MAX = mt;		
					append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );
				}
				else{
					SS_MAX = ms;		
					append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );														
				}
			}
			else{
				SS_MAX = ms;		
				append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );									
			}
		}

		else if( gt == gq ){

			if( gs < gt ){
				SS_MAX = ms;		
				append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );	
			}
			else if( gt == gs ){
				if( mt > mq ){
					if( mt > ms ){
						SS_MAX = mt;		
						append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );	
					}
					else if( mt == ms && (rand()%2) ){
						SS_MAX = mt;		
						append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );							
					}
					else{
						SS_MAX = ms;		
						append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );		
					}
				}
				else if( mt == mq ){
					if( ms > mt ){
						SS_MAX = ms;		
						append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );							
					}
					else if( mt == ms ){
						int var = rand()%3;
						if( var == 0){
							SS_MAX = ms;		
							append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );														
						}
						else if (var == 1){
							SS_MAX = mt;		
							append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );							
						}
						else{
							SS_MAX = mq;		
							append_node( ptr_res, mq, gq , Q, 'd' );	
						}
					}
					else{
						if( rand()%2 ){
							SS_MAX = mt;		
							append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );							
						}
						else{							
							SS_MAX = mq;		
							append_node( ptr_res, mq, gq , Q, 'd' );	
						}
					}
				}
				else{		// gt == gs == gq && mq > mt
					if( mq > ms ){
						SS_MAX = mq;		
						append_node( ptr_res, mq, gq , Q, 'd' );	
					}
					else if( mq == ms && (rand()%2) ){
						SS_MAX = mq;		
						append_node( ptr_res, mq, gq , Q, 'd' );							
					}
					else{
						SS_MAX = ms;		
						append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );	
					}
				}
			}
			else{			// gt==gq && gt < gs
				if( mt > mq ){
					SS_MAX = mt;		
					append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );	
				}
				else if( mt == mq && (rand()%2) ){
					SS_MAX = mt;		
					append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );	
				}
				else{
					SS_MAX = mq;		
					append_node( ptr_res, mq, gq , Q, 'd' );	
				}
			}
		}

		else{

			if( gq < gs ){
				SS_MAX = mq;		
				append_node( ptr_res, mq, gq , Q, 'd' );					
			}
			else if( gq == gs ){
				if( mq > ms ){
					SS_MAX = mq;		
					append_node( ptr_res, mq, gq , Q, 'd' );										
				}
				else if( mq == ms && (rand()%2) ){
					SS_MAX = mq;		
					append_node( ptr_res, mq, gq , Q, 'd' );										
				}
				else{					
					SS_MAX = ms;		
					append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );														
				}
			}
			else{
				SS_MAX = ms;		
				append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );									
			}
		}
		/* End of minimum #gaps score selection */

		ptr_res = ptr_res->next;

		/* Exclude dominated scores, advancing iterators, based on the value SS_MAX */
		while( T && T->score.matches <= SS_MAX ) T = T->next;
		while( S && S->score.matches <= SS_MAX ) S = S->next;
		while( Q && Q->score.matches + subs_score <= SS_MAX ) Q = Q->next;

	}

	return res;

}

void remove_list( Node * ptr ){
	Node * tmp;
	while( ptr->next ){		
		tmp = ptr;
		ptr = ptr->next;
		free(tmp);
	}
	free(ptr);
}