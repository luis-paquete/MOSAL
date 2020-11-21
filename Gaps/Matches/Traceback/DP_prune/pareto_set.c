/****************************************************************************

 It computes the set of non dominated scores for the #Matches-#Gaps problem.
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
 * Implementation of the functions in 'pareto_set.h'.
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "pareto_set.h"
#include "bounds.h"

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
	 * Until all the iterators reach the end, it is selected the score with
	 * minimum #gaps and, in case of equality, with maximum #matches, from
	 * the two scores: p1-(mat_p1,gap_p1), p2-(mat_p2,gap_p2).
	 * (If the score vectors are exactly the same, than one is chosen randomly.)
	 * This score is appended to the result as it is always non dominated.
	 * The #matches of that score (to be stored in M_MAX) is also used to 
	 * eliminate the (dominated) scores following in each set, by advancing 
	 * the iterators until the #matches of the pointed score is greater than M_MAX.
	 * The dominance of these scores is also assured because the sets are ordered by #gaps
	 * and they all have greater #gaps compared to the selected score.
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

	/* defines the number of matches of the score to be inserted */
	int M_MAX;
	/* auxiliar variables for comparisons */
	int gap_p1, gap_p2;	
	int mat_p1, mat_p2;		

	while( p1 || p2 ){

		/*  Initialization of the auxiliar variables for the p1 and p2 cells.
			In case iterators reach the end of its set (first comparison),
			its values are not inserted, keeping the algorithm structure */
		/* P1 */
		if( p1 ){
			gap_p1 = p1->score.gaps + IND_PENALTY;		/* add indel */
			mat_p1 = p1->score.matches;
		}
		else{
			gap_p1 = INT_MAX;
			mat_p1 = INT_MIN;			
		}

		/* P2 */
		if( p2 ){
			gap_p2 = p2->score.gaps + 1;				/* open a new gap */
			mat_p2 = p2->score.matches;
		}
		else{
			gap_p2 = INT_MAX;
			mat_p2 = INT_MIN;			
		}		
		/* End of auxiliar variables initialization */


		/*  Select score with minimum #gaps and, in case of equality,
			with maximum #matches */
		if( gap_p1 < gap_p2 ){
			M_MAX = mat_p1;
			append_node( ptr_res, mat_p1, gap_p1, p1, tb_dir);			
		}
		else if ( gap_p1 == gap_p2 ){
			if( mat_p1 > mat_p2 ){
				M_MAX = mat_p1;
				append_node( ptr_res, mat_p1, gap_p1, p1, tb_dir);
			}
			else if( mat_p1 == mat_p2 && (rand()%2) ){
				M_MAX = mat_p1;
				append_node( ptr_res, mat_p1, gap_p1, p1, tb_dir);				
			}
			else{				
				M_MAX = mat_p2;
				append_node( ptr_res, mat_p2, gap_p2, p2, tb_dir);	
			}
		}
		else{
			M_MAX = mat_p2;
			append_node( ptr_res, mat_p2, gap_p2, p2, tb_dir);	
		}
		/* End of minimum #gaps score selection */

		ptr_res = ptr_res->next;	
		
		/* Exclude dominated scores, advancing iterators, based on the value M_MAX */
		while( p1 && p1->score.matches <= M_MAX ) p1 = p1->next;
		while( p2 && p2->score.matches <= M_MAX ) p2 = p2->next;
	}

	return res;
}

Pareto_set pareto_merge3(Pareto_set S, Pareto_set T, Pareto_set Q, int match, int i, int j){
	/**
	 * Scores are assumed to be ordered incrementally by their #gaps.
	 * Until all the iterators reach the end, it is selected the score with
	 * minimum #gaps and, in case of equality, with maximum #matches, from
	 * the three scores: S-(ms,gs), Q-(mq,gq) and T-(mt,gt).
	 * (If the score vectors are exactly the same, than one is chosen randomly.)
	 * This score is appended to the result as it is always non dominated.
	 * The #matches of that score (to be stored in M_MAX) is also used to 
	 * eliminate the (dominated) scores following in each set, by advancing 
	 * the iterators until the #matches of the pointed score is greater than M_MAX.
	 * The dominance of these scores is also assured because the sets are ordered by #gaps
	 * and they all have greater #gaps compared to the selected score.
	 */		
	
	/* allocate memory for the result set */
	Pareto_set res = create_list();
	Node * ptr_res = res;

	char dir; 		/* in order to improve the pruning */

	/* defines the number of matches of the score to be inserted */
	int M_MAX;
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
		if( Q ){
			gq = Q->score.gaps;
			mq = Q->score.matches + match;
		}
		else{
			gq = INT_MAX;
			mq = INT_MIN;
		}
		/* S */
		if( S ){
			gs = S->score.gaps;
			ms = S->score.matches;
		}
		else{
			gs = INT_MAX;
			ms = INT_MIN;
		}
		/* T */
		if( T ){
			gt = T->score.gaps;
			mt = T->score.matches;
		}
		else{
			gt = INT_MAX;
			mt = INT_MIN;
		}
		/* End of auxiliar variables initialization */			


		/*  Select score with minimum #gaps and, in case of equality,
			with maximum #matches */
		if( gt < gq ){

			if( gt < gs ){
				M_MAX = mt;		
				append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );					
				dir= 'u';
			}
			else if( gt == gs ){
				if( mt > ms ){
					M_MAX = mt;		
					append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );		
					dir= 'u';								
				}
				else if( mt == ms && (rand()%2) ){
					M_MAX = mt;		
					append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );
					dir= 'u';
				}
				else{
					M_MAX = ms;		
					append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );	
					dir= 'l';													
				}
			}
			else{
				M_MAX = ms;		
				append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );		
				dir= 'l';							
			}
		}

		else if( gt == gq ){

			if( gs < gt ){
				M_MAX = ms;		
				append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );	
				dir= 'l';
			}
			else if( gt == gs ){
				if( mt > mq ){
					if( mt > ms ){
						M_MAX = mt;		
						append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );	
						dir= 'u';
					}
					else if( mt == ms && (rand()%2) ){
						M_MAX = mt;		
						append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );	
						dir= 'u';						
					}
					else{
						M_MAX = ms;		
						append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );	
						dir= 'l';	
					}
				}
				else if( mt == mq ){
					if( ms > mt ){
						M_MAX = ms;		
						append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );	
						dir= 'l';						
					}
					else if( mt == ms ){
						int var = rand()%3;
						if( var == 0){
							M_MAX = ms;		
							append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );
							dir= 'l';														
						}
						else if (var == 1){
							M_MAX = mt;		
							append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );	
							dir= 'u';						
						}
						else{
							M_MAX = mq;		
							append_node( ptr_res, mq, gq , Q, 'd' );	
							dir= 'd';
						}
					}
					else{
						if( rand()%2 ){
							M_MAX = mt;		
							append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );
							dir= 'u';							
						}
						else{							
							M_MAX = mq;		
							append_node( ptr_res, mq, gq , Q, 'd' );	
							dir= 'd';
						}
					}
				}
				else{		/* gt == gs == gq && mq > mt */
					if( mq > ms ){
						M_MAX = mq;		
						append_node( ptr_res, mq, gq , Q, 'd' );	
						dir= 'd';
					}
					else if( mq == ms && (rand()%2) ){
						M_MAX = mq;		
						append_node( ptr_res, mq, gq , Q, 'd' );	
						dir= 'd';						
					}
					else{
						M_MAX = ms;		
						append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );	
						dir= 'l';
					}
				}
			}
			else{			/* gt==gq && gt < gs */
				if( mt > mq ){
					M_MAX = mt;		
					append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );	
					dir= 'u';
				}
				else if( mt == mq && (rand()%2) ){
					M_MAX = mt;		
					append_node( ptr_res, mt, gt , T->traceback_ptr, 'u' );	
					dir= 'u';
				}
				else{
					M_MAX = mq;		
					append_node( ptr_res, mq, gq , Q, 'd' );
					dir= 'd';	
				}
			}
		}

		else{

			if( gq < gs ){
				M_MAX = mq;		
				append_node( ptr_res, mq, gq , Q, 'd' );		
				dir= 'd';			
			}
			else if( gq == gs ){
				if( mq > ms ){
					M_MAX = mq;		
					append_node( ptr_res, mq, gq , Q, 'd' );	
					dir= 'd';									
				}
				else if( mq == ms && (rand()%2) ){
					M_MAX = mq;		
					append_node( ptr_res, mq, gq , Q, 'd' );										
					dir= 'd';
				}
				else{					
					M_MAX = ms;		
					append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );	
					dir= 'l';													
				}
			}
			else{
				M_MAX = ms;		
				append_node( ptr_res, ms, gs , S->traceback_ptr, 'l' );	
				dir= 'l';								
			}
		}
		/* End of minimum #gaps score selection */

		/* prune last score appended */
		prune(ptr_res, i, j, dir);

		if( ptr_res->next )			/* advance ptr only if it was not prunned */
			ptr_res = ptr_res->next;

		/* Exclude dominated scores, advancing iterators, based on the value M_MAX */
		while( T && T->score.matches <= M_MAX ) T = T->next;
		while( S && S->score.matches <= M_MAX ) S = S->next;
		while( Q && Q->score.matches + match <= M_MAX ) Q = Q->next;

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