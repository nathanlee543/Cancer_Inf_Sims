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

#define write_mut_freq 1
double true_params_out[7];
double surviving_runs;

ofstream myfile;

double expansion(double b[2], double d[2], double t1, double t, double delta, double u, int max_clones)
{
    restart_sim:
        
        int i,j,k,driver;
        int estimates_length = 15;
        double estimates_out[estimates_length];

        // dynamically allocate memory for large arrays. 
        // Since we're compiling with g++, use c++ syntax
        int (*where)[2] = new int[max_clones][2];
        int (*ancestor)[2] = new int[max_clones][2];
        int (*surviving_ancestor)[2] = new int[max_clones][2];
        double (*Clone)[2] = new double[max_clones][2];
        double (*surviving_clone)[2] = new double[max_clones][2];
        int m_int = 0;

        
        double cells[2];
        
        double f1, f2, nf1, nf2, min_SNV_count, alpha1, alpha2,Time,time_scale, extra_m_analytic;
        int num_clones[2] = {0,0};
        int num_surviving_clones[2];
        double rand;
        int t1_yet = 0;
        int t_yet[2] = {0,0};
        double M1 = 0.0;
        double M2 = 0.0;
        double m_clonal = 0.0, m = 0.0, subclonal = 0.0;
        double m_obs;
        int clone_ind;
        double r_map, r1_map, b1_map, b_map, t1_map, t_map, u_estimate;
        double time_init[2] = {0.0,t1};
        int m_mutations[1000];
        double alpha2_min = 0.1, alpha2_max = 0.9; 

        // initializing 
        for (driver=0; driver<2;driver++)
        {
            cells[driver] = 1.0;
            Clone[0][driver]=1.0;
            for(i=1; i<max_clones; i++)
            {
                Clone[i][driver]=0.0;
                surviving_clone[i][driver]=0.0;
                ancestor[i][driver]=0;
                surviving_ancestor[i][driver]=0;
            }
        }

        for (driver=0; driver<2;driver++)
        {
            Time = time_init[driver];
            time_scale = cells[driver]*(b[driver]+d[driver]+u);
            Time += gsl_ran_exponential(r,1/time_scale);

            for (; Time < (t1 + t + delta) && num_clones[driver]<max_clones ; )
            {
                if (Time >= t1 && t1_yet == 0)
                {   
                    //printf("t1\n");
                    t1_yet = 1;
                    cells[driver] -= 1.0;
                    // get number of passengers present in founder type-1 cell
                    // pick which type-0 gets the driver mutation
                    rand=gsl_ran_flat(r,0.0,cells[driver]);
                    
                    for(j=0; j<=num_clones[driver]; j++)
                    {
                        if (rand<Clone[j][driver])
                        {
                            clone_ind = j;
                            m_mutations[m_int] = clone_ind;
                            j=num_clones[driver]+1;
                        }
                        else rand-=Clone[j][driver];
                    }
                    Clone[clone_ind][driver] -= 1.0;
                    // now find how many passengers clone_ind has at t1 (m)
                    for (; clone_ind > 0;)
                    {
                        clone_ind = ancestor[clone_ind][driver];
                        m += 1.0;
                        m_int += 1;
                        m_mutations[m_int] = clone_ind;

                        if (clone_ind == 0) printf("true m: %f \n", m);
                    }
                }

                if (Time >= (t1+t) && t_yet[driver] == 0)
                {

                    t_yet[driver] = 1;
                    if (driver == 1)
                    {
                        alpha1 = cells[1]/(M1+cells[1]);
                    }
                    M1 += cells[driver];
                }

                // determine which subclone event happens to
                // update num cells in that subclone
                rand=gsl_ran_flat(r,0.0,cells[driver]);
                for (j=0; j<=num_clones[driver]; j++)
                {
                    if (rand<Clone[j][driver])
                    {

                        if (rand<u*Clone[j][driver]/(b[driver]+d[driver]+u))
                        {
                            //new clone
                            num_clones[driver]++;
                            Clone[num_clones[driver]][driver]=1.0;
                            ancestor[num_clones[driver]][driver]=j;
                            Clone[j][driver] -= 1.0;
                            j=num_clones[driver]+1;
                            
                        }
                        else
                        {
                            if (rand<(b[driver]+u)*Clone[j][driver]/(b[driver]+d[driver]+u))
                            {
                                //one j is born
                                cells[driver]+=1.0;
                                Clone[j][driver]+=1.0;
                                j=num_clones[driver]+1;
                            }
                            else
                            {
                                //one j dies
                                cells[driver]-=1.0;
                                Clone[j][driver]-=1.0;
                                j=num_clones[driver]+1;
                            }
                        }
                    }
                    else rand-=Clone[j][driver];
                    
                }
                if (cells[driver] == 0.0)
                {
                    cells[driver] = 1.0;
                    Clone[0][driver]=1.0;
                    num_clones[driver]= 0;
                    Time = time_init[driver];

                    for(i=1; i<max_clones; i++)
                    {
                        Clone[i][driver]=0.0;
                        surviving_clone[i][driver]=0.0;
                        ancestor[i][driver]=0;
                        surviving_ancestor[i][driver]=0;
                    }
                }
                time_scale=cells[driver]*(b[driver]+d[driver]+u);
                Time+=gsl_ran_exponential(r,1/time_scale);
            }

        }

        for (driver=0; driver<2;driver++)
        {
            // clone[i] is going to include all the cells that have mutation i
            for(i=num_clones[driver]; i>0; i--) //first add clones then remove non-zero ones
            {
                Clone[ancestor[i][driver]][driver]+=Clone[i][driver];
            }
            
            // add cells[1] to subclone size for the m mutations
            for (i = 0; i <= m_int; i++)
            {
                Clone[m_mutations[i]][0] += cells[1];
            }

            num_surviving_clones[driver]=0;
            
            for(i=0; i<=num_clones[driver]; i++)
            { 
                // if ancester mut is extinct, change ancester to previous ancester in lineage
                if (Clone[ancestor[i][driver]][driver]==0.0)
                {
                    ancestor[i][driver]=ancestor[ancestor[i][driver]][driver];
                }
                // create array of only surviving clones
                // where[i] is index of same clone in re-indexed surviving array
                // survivng_ancester is ancester array, with re-idnexed subclones
                if (Clone[i][driver]>0 || i==0)
                {
                    surviving_clone[num_surviving_clones[driver]][driver]=Clone[i][driver];
                    where[i][driver]=num_surviving_clones[driver];
                    surviving_ancestor[num_surviving_clones[driver]][driver]=where[ancestor[i][driver]][driver];
                    num_surviving_clones[driver]+=1;
                }
                
            }
            num_surviving_clones[driver]-=1;
        }
        // prevent a stack overflow error if you have more clones than you initialized    
        if (num_clones[0]<max_clones && num_clones[1]<max_clones)    
        {

            M2 = cells[0]+cells[1];
            f1 = 0.01;
            f2 = 0.2;
            nf1 = f1*M2;
            nf2 = f2*M2;

            // the minimum frequency of a variant that will be saved to file is 0.1%
            min_SNV_count = 0.001*M2;

            for (driver=0; driver<2;driver++)
            {
                for (i=0; i<num_surviving_clones[driver]; i++)
                {
                    if (driver == 1 && surviving_clone[i][driver] == cells[driver]) m_clonal += 1.0;
                    if (surviving_clone[i][driver]>nf1 && surviving_clone[i][driver]<nf2) subclonal += 1.0;

                }
            }
            alpha2 = cells[1]/M2;
            m_obs = m + m_clonal;
            printf("subclonal = %f bn (%f,%f)\n",subclonal, f1, f2);
            
            // Compute estimates
            // compute estimates of runs for all alpha2 values, and save them to file
            // for purposes relating to parallelization, function will only return if 0.1 < alpha2 < 0.9
            // so that we have 100 runs where 0.1 < alpha2 < 0.9
            r1_map = (1./delta)*log(alpha2*M2/(alpha1*M1));
            r_map =  (1./delta)*log((1.-alpha2)*M2/((1.-alpha1)*M1));
            

            u_estimate = subclonal*pow(((alpha2/r1_map + (1.-alpha2)/r_map)*(1./f1 - 1./f2)),-1);
            extra_m_analytic = u_estimate/r1_map;
            t1_map = (m_obs-extra_m_analytic)/u_estimate;

            t_map = /* TO DO: fill in with estimates */
            b_map = r_map*exp(-r_map*(t_map + t1_map + delta))* (1 - alpha2)*M2;
            b1_map = r1_map*exp(-r1_map*(t_map+delta))*alpha2*M2;

            
            estimates_out[0] = r1_map;
            estimates_out[1] = r_map;
            estimates_out[2] = u_estimate;
            estimates_out[3] = t1_map;
            estimates_out[4] = t_map;
            estimates_out[5] = b_map;
            estimates_out[6] = b1_map;
            estimates_out[7] = M1;
            estimates_out[8] = alpha1;
            estimates_out[9] = M2;
            estimates_out[10] = alpha2;
            estimates_out[11] = subclonal;
            estimates_out[12] = (1/r1_map)*log(M1*alpha1) - (0.5772156649/r1_map);
            estimates_out[13] = m;
            estimates_out[14] = m_obs;
   
            // write the mutation frequencies of all surviving clones with frequency greater than f1
            if (write_mut_freq == 1) 
            {
                #pragma omp critical (mut_freq)
                    myfile.open("run_i_surviving_mut_j_freq_type0.txt",fstream::app);
                    for(i=0; i<num_surviving_clones[0]; i++)
                    {
                        if (surviving_clone[i][0] > min_SNV_count)
                        {
                            myfile << surviving_clone[i][0] <<" ";
                        }
                    }
                    myfile << endl;
                    myfile.close();
                    // start line with the run that it came from
                    myfile.open("run_i_surviving_mut_j_freq_type1.txt",fstream::app);
                    for(i=0; i<num_surviving_clones[1]; i++)
                    {
                        if (surviving_clone[i][1] > min_SNV_count)
                        {
                            myfile << surviving_clone[i][1] <<" ";
                        }

                    }
                    myfile << endl;
                    myfile.close();

                    char filename[sizeof "r1_r_u_t1_t_b_b1_M1_alpha1_M2_alpha2_gamma_1000cells_f1_p01_f2_p2.txt"];
                    sprintf(filename, "r1_r_u_t1_t_b_b1_M1_alpha1_M2_alpha2_gamma_f1_p01_f2_p2.txt");
                    
                    myfile.open(filename,fstream::app);  
                    for (j=0; j<estimates_length; j++)
                    {
                        myfile << estimates_out[j] <<" ";
                    }
                    myfile << endl;
                    myfile.close();
                     
            }
            
            // free memory
            delete[] Clone;
            delete[] surviving_clone;
            delete[] where;
            delete[] ancestor;
            delete[] surviving_ancestor;

            
            
            // only return function if .1 < alpha2 < .9
            if (alpha2 > alpha2_min && alpha2 < alpha2_max)
            {
                return 1.0;
            }
            else
            {   
                goto restart_sim;
            }
        }
        else
        {
            // this will append lines to a file that shows the number of clones 
            // really more of a peace-of-mind check
            // having runs where you restart for having more clones than initialized will bias results
            printf("max clones reached, increase size of max_clones");
            delete[] Clone;
            delete[] surviving_clone;
            delete[] where;
            delete[] ancestor;
            delete[] surviving_ancestor;
            exit(EXIT_FAILURE);
        }
}


int main(int argc, char *argv[]){
    if (argc < 2){
        printf("didn't supply filename");
        exit(EXIT_FAILURE);
    }
    else{
        
        FILE * pFile;
        char * filename = argv[1];
        int param_count = 12;
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
                    printf("comment removed");
                    cout << str[i][j] << endl;
                    j = 100;
                } else{
                    cout << str[i][j] << endl;
                    str2[i][j] = str[i][j];
                }
            }
        }
        printf("print str2[0]:%s\n",str[0]);
        
        double b[2],d[2],t1,t,delta,u;
        b[0] = stod(str2[0],NULL);
        b[1] = stod(str2[1],NULL);
        d[0] = stod(str2[2],NULL);
        d[1] = stod(str2[3],NULL);
        t1 = stod(str2[4],NULL);
        t = stod(str2[5],NULL);
        delta = stod(str2[6],NULL);
        u = stod(str2[7],NULL);

        int num_threads = strtol(str2[8],NULL,10);
        int runs = strtol(str2[9],NULL,10);
        int max_clones = strtol(str2[10],NULL,10);
        int par0_single_run1 = strtol(str2[11],NULL,10);

        printf("argc:%i\n",argc);
    
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        
        true_params_out[0] = (b[1] - d[1]);
        true_params_out[1] = (b[0]-d[0]);
        true_params_out[2] = u;
        true_params_out[3] = t1;
        true_params_out[4] = t;
        true_params_out[5] = b[0];
        true_params_out[6] = b[1];

        int counter;
        
        // add a comment header to the estimates file

        char filename[sizeof "r1_r_u_t1_t_b_b1_M1_alpha1_M2_alpha2_gamma_1000cells_f1_p01_f2_p2.txt"];
        sprintf(filename, "r1_r_u_t1_t_b_b1_M1_alpha1_M2_alpha2_gamma_f1_p01_f2_p2.txt");
        
        myfile.open(filename,fstream::app);  

        myfile << "r1 r u t1 t_ssc b b1 M1 alpha1 M2 alpha2 gamma t_bd m m_obs" <<" ";
        
        myfile << endl;
        myfile.close();
	    


        // input file determines if we multithread
        if (par0_single_run1 == 0)
        {
            surviving_runs=0.0;
            int completed = 0;
            omp_set_dynamic(0);     // Explicitly disable dynamic teams
            omp_set_num_threads(num_threads); // set number of threads for multithreading
            #pragma omp parallel for
                for(counter=0; counter<runs;counter++ )
                {
                    surviving_runs+=expansion(b, d, t1, t, delta, u, max_clones);
                    
                    #pragma omp critical(PRINT)
                    {
                        completed++;
                        cout << completed << " / " << runs << endl; 
                    }
                }
        } else // otherwise, do a single run
        {
            expansion(b, d, t1, t, delta, u, max_clones);
        }
        // save the true parameter values to textfile
        myfile.open("true_params.txt");
        for (j=0; j<7; j++) myfile << true_params_out[j] <<" ";
        myfile.close();
    }

    return 0;
    
}
