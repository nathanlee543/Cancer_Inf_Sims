//
//  monte_carlo_tumor_size_benchmarking.c
//  
//  Code to accompany manuscript "Inferring parameters of cancer evolution from sequencing and clinical data" by Lee and Bozic.
//  Created by Nathan Lee, 2022
// Monte Carlo simulation of tumor evolution
// mutations occur over time, not only at divisions
// cell division with rate b, death with rate d, mutation with rate u.
//  
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <omp.h>

using namespace std;

const gsl_rng_type * T; 
gsl_rng * r;

double true_params_out[2];
double surviving_runs;

ofstream myfile;

double expansion(double b, double d, FILE *fp)
{    
    double cells;
    double Time,time_scale;
    double rand;
    double t_estimate, t_true;
    double time_init = 0.0;
    int tumor_size_length = 8;
    double tumor_size[tumor_size_length] = {pow(10.0,3),pow(10.0,4),pow(10.0,5),pow(10.0,6),pow(10.0,7),pow(10.0,8),pow(10.0,9), pow(10.0,10)};
    
    // initializing
    int tumor_size_i = 0;
    cells = 1.0;
    Time = time_init;

    // draw random samples to enter for loop
    time_scale = cells*(b+d);
    Time += gsl_ran_exponential(r,1/time_scale);
    rand=gsl_ran_flat(r,0.0,cells);

    // birth-death process with one type
    // run until we hit max tumor size
    for (; cells < tumor_size[tumor_size_length-1]; )
    {
        if (rand<b*cells/(b+d))
        {
            //birth
            cells+=1.0;
        }
        else
        {
            //death
            cells-=1.0;
        }           

        // if population goes extinct restart simulation run
        if (cells == 0.0)
        {
            // initializing
            tumor_size_i = 0;
            cells = 1.0;
            Time = time_init;
        }

        // if cell count is one of the specified sizes, save estimate and time
        if (cells == tumor_size[tumor_size_i])
        {
            // compute estimate for t and record the true value
            t_estimate = (1./(b-d))*log(cells);
            t_true = Time;
            // write results to output file as csv
            fp = fopen("t_benchmarking_simulation_results.txt","a");    
            fprintf(fp, "%f,%f,%f\n",t_estimate, t_true, tumor_size[tumor_size_i]);
            fclose(fp); 

            tumor_size_i += 1;
        }

        // sample for the next time step and event type
        time_scale = cells*(b+d);
        Time += gsl_ran_exponential(r,1/time_scale);
        rand=gsl_ran_flat(r,0.0,cells);

    }    
    return 1.0;
}


int main(int argc, char *argv[]){
    if (argc < 2){
        printf("didn't supply filename");
        exit(EXIT_FAILURE);
    }
    else{
        
        // processing the input files
        FILE * pFile;
        FILE * fp;
        char * filename = argv[1];
        int param_count = 5;
        char  str[param_count][100];
        char  str2[param_count][100];
        pFile = fopen(filename, "r");

        if (pFile == NULL){
            perror("Error opening file");
            exit(EXIT_FAILURE);
        } 
        else{
            int i;
            for (i=0;i<param_count;i++){
                fscanf(pFile, "%s", str[i]); 
                cout << str[i] << endl;
            }   

        }
        fclose(pFile);
        int i,j;
        // seed random number
        for (i = 0; i<param_count;i++){
            cout << str[i] << endl;
            for (j = 0;j<100;j++){
                
                if (str[i][j]=='#'){
                    cout << str[i][j] << endl;
                    j = 100;
                } else{
                    cout << str[i][j] << endl;
                    str2[i][j] = str[i][j];
                }
            }
        }
        
        double b,d;
        b = stod(str2[0],NULL);
        d = stod(str2[1],NULL);
        int num_threads = strtol(str2[2],NULL,10); 
        int runs = strtol(str2[3],NULL,10);
        int par0_single_run1 = strtol(str2[4],NULL,10);
    
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        
        true_params_out[0] = b;
        true_params_out[1] = d;

        int counter;
        
        // add a header to the output file
        fp = fopen("t_benchmarking_simulation_results.txt","w");    
        fprintf(fp,"t_estimate,t_true,tumor_size\n");
        fclose(fp);
           
        // input file determines if we multithread
        if (par0_single_run1 == 0)
        {
            surviving_runs=0.0;
            int completed = 0;
            omp_set_dynamic(0);     // Explicitly disable dynamic teams
            omp_set_num_threads(num_threads); // set number of threads for multithreading
            #pragma omp parallel for shared(completed)
                for(counter=0; counter<runs;counter++ )
                {
                    surviving_runs+=expansion(b, d, fp);
                    //printf("finished job %f\n",surviving_runs);

                    #pragma omp critical(PRINT)
                    {
                        completed++;
                        cout << completed << " / " << runs << endl; 
                    }
                }
        } else // otherwise, do a single run
        {
            surviving_runs=expansion(b, d, fp);
            printf("finished job %f\n",surviving_runs);
        }
        // save the true parameter values to textfile
        myfile.open("true_b_d.txt");
        for (j=0; j<2; j++) myfile << true_params_out[j] <<" ";
        myfile.close();
    }
    return 0;
}
