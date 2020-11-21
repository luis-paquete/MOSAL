gap = Gaps
ind = Indels
mat = Matches
ss = Subs_score
tb = Traceback
ntb = No_traceback
dp = DP
dpp = DP_prune

all:
	make -C $(gap)/$(mat)/$(ntb)/$(dp)
	make -C $(gap)/$(mat)/$(ntb)/$(dpp)

	make -C $(gap)/$(mat)/$(tb)/$(dp)
	make -C $(gap)/$(mat)/$(tb)/$(dpp)

	make -C $(gap)/$(ss)/$(ntb)/$(dp)
	make -C $(gap)/$(ss)/$(ntb)/$(dpp)

	make -C $(gap)/$(ss)/$(tb)/$(dp)
	make -C $(gap)/$(ss)/$(tb)/$(dpp)

	make -C $(ind)/$(mat)/$(ntb)/$(dp)
	make -C $(ind)/$(mat)/$(ntb)/$(dpp)

	make -C $(ind)/$(mat)/$(tb)/$(dp)
	make -C $(ind)/$(mat)/$(tb)/$(dpp)

	make -C $(ind)/$(ss)/$(ntb)/$(dp)
	make -C $(ind)/$(ss)/$(ntb)/$(dpp)

	make -C $(ind)/$(ss)/$(tb)/$(dp)
	make -C $(ind)/$(ss)/$(tb)/$(dpp)

	gcc main.c -Wall -o mosal

clean:
	rm -f $(gap)/$(mat)/$(ntb)/$(dp)/*.o $(gap)/$(mat)/$(ntb)/$(dp)/prog
	rm -f $(gap)/$(mat)/$(ntb)/$(dpp)/*.o $(gap)/$(mat)/$(ntb)/$(dpp)/prog

	rm -f $(gap)/$(mat)/$(tb)/$(dp)/*.o $(gap)/$(mat)/$(tb)/$(dp)/prog
	rm -f $(gap)/$(mat)/$(tb)/$(dpp)/*.o $(gap)/$(mat)/$(tb)/$(dpp)/prog

	rm -f $(gap)/$(ss)/$(ntb)/$(dp)/*.o $(gap)/$(ss)/$(ntb)/$(dp)/prog
	rm -f $(gap)/$(ss)/$(ntb)/$(dpp)/*.o $(gap)/$(ss)/$(ntb)/$(dpp)/prog

	rm -f $(gap)/$(ss)/$(tb)/$(dp)/*.o $(gap)/$(ss)/$(tb)/$(dp)/prog
	rm -f $(gap)/$(ss)/$(tb)/$(dpp)/*.o $(gap)/$(ss)/$(tb)/$(dpp)/prog

	rm -f $(ind)/$(mat)/$(ntb)/$(dp)/*.o $(ind)/$(mat)/$(ntb)/$(dp)/prog
	rm -f $(ind)/$(mat)/$(ntb)/$(dpp)/*.o $(ind)/$(mat)/$(ntb)/$(dpp)/prog

	rm -f $(ind)/$(mat)/$(tb)/$(dp)/*.o $(ind)/$(mat)/$(tb)/$(dp)/prog
	rm -f $(ind)/$(mat)/$(tb)/$(dpp)/*.o $(ind)/$(mat)/$(tb)/$(dpp)/prog

	rm -f $(ind)/$(ss)/$(ntb)/$(dp)/*.o $(ind)/$(ss)/$(ntb)/$(dp)/prog
	rm -f $(ind)/$(ss)/$(ntb)/$(dpp)/*.o $(ind)/$(ss)/$(ntb)/$(dpp)/prog

	rm -f $(ind)/$(ss)/$(tb)/$(dp)/*.o $(ind)/$(ss)/$(tb)/$(dp)/prog
	rm -f $(ind)/$(ss)/$(tb)/$(dpp)/*.o $(ind)/$(ss)/$(tb)/$(dpp)/prog

	rm -f mosal
