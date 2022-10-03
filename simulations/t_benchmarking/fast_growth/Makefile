parallel_paper: t_benchmarking_with_size.c
	g++ -g -O3 -Wall -o MCSim t_benchmarking_with_size.c -lgsl -lgslcblas -fopenmp

clean:
	rm r1*
	rm run_i_surviving*
	rm *MCSim
	rm true_params.txt

#run by doing ./MCSim run_inputs.txt

