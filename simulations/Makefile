parallel_paper: simulations_for_paper.c
	g++ -g -O3 -Wall -o MCSim  simulations_for_paper.c -lgsl -lgslcblas -fopenmp

clean:
	rm r1*
	rm run_i_surviving*
	rm *MCSim
	rm true_params.txt

#run by doing ./MCSim run_inputs.txt

