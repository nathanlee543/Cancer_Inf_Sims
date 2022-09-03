parallel_paper: two_nested_subclones_simulation.c
	g++ -g -O3 -Wall -o MCSim  two_nested_subclones_simulation.c -lgsl -lgslcblas -fopenmp

clean:
	rm r1*
	rm run_i_surviving*
	rm *MCSim
	rm true_params.txt

#run by doing ./MCSim run_inputs.txt

