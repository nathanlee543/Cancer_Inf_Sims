//
//  two_sibling_subclones_simulation.c
//  
//  Code to accompany manuscript "Inferring parameters of cancer evolution from sequencing and clinical data" by Lee and Bozic.
//  Created by Nathan Lee, 2022
// Monte Carlo simulation of tumor with 2 sibling subclones
// mutations occur over time, not only at divisions 
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <omp.h>

using namespace std;

const gsl_rng_type * T;
gsl_rng * r;

double true_params_out[8];
double surviving_runs;

ofstream myfile;

double expansion(double b[3], double d[3], double t1, double t2, double t, double delta, double u, int max_clones, FILE *fp)
{
    int i,j,driver;

    // dynamically allocate memory for large arrays. 
    // Since we're compiling with g++, use c++ syntax
    int (*where)[3] = new int[max_clones][3];
    int (*ancestor)[3] = new int[max_clones][3];
    int (*surviving_ancestor)[3] = new int[max_clones][3];
    double (*Clone)[3] = new double[max_clones][3];
    double (*surviving_clone)[3] = new double[max_clones][3];
            
    double cells[3];
    double alpha1, alpha2, beta1, beta2;
    double beta2_min = 0.1, beta2_max = 0.9, beta1_min = 0.1, beta1_max = 0.9; 
    double f1, f2, nf1, nf2,Time,time_scale;
    int num_clones[3] = {0,0,0};
    int num_surviving_clones[3];
    double rand;
    int t1_yet = 0;
    int t2_yet = 0;
    int t_yet[3] = {0,0,0};
    double M1 = 0.0, M2 = 0.0;
    double WBC1[3] = {0.0,0.0,0.0};
    double m1_clonal = 0.0, m2_clonal = 0.0, m1 = 0.0, m2 = 0.0, subclonal = 0.0;
    int m1_int = 0, m2_int = 0;
    double m1_obs, m2_obs;
    int clone_ind1, clone_ind2;
    double r_hat, r1_hat, r2_hat, t1_hat, t2_hat, tau1_hat, tau2_hat, u_hat, u_top, u_bottom, m1_hat, m2_hat;
    double time_init[3] = {0.0,t1,t2};
    int m1_mutations[1000], m2_mutations[1000];

    // initializing 
    for (driver=0; driver<3;driver++)
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

    // main loop for growing populations
    for (driver=0; driver<3;driver++)
    {
        Time = time_init[driver];
        time_scale = cells[driver]*(b[driver]+d[driver]+u);
        Time += gsl_ran_exponential(r,1/time_scale);

        for (; Time < (t1 + t + delta) && num_clones[driver]<max_clones ; )
        {
            // if we just hit t1 in last time step, record values
            // this condition will be satisfied during driver = 0 loop
            if (Time >= t1 && t1_yet == 0)
            {   
                t1_yet = 1;
                cells[driver] -= 1.0; // equivalent to cells[0] -= 1.0;

                // get number of passengers present in founder type-1 cell
                // pick which type-0 gets the driver mutation
                rand=gsl_ran_flat(r,0.0,cells[driver]);
                
                for(j=0; j<=num_clones[driver]; j++)
                {
                    if (rand<Clone[j][driver])
                    {
                        clone_ind1 = j;
                        m1_mutations[m1_int] = clone_ind1;
                        j=num_clones[driver]+1;
                    }
                    else rand-=Clone[j][driver];
                }
                Clone[clone_ind1][driver] -= 1.0;
                // now find how many passengers clone_ind1 has at t1 (m)
                for (; clone_ind1 > 0;)
                {
                    clone_ind1 = ancestor[clone_ind1][driver];
                    m1 += 1.0;
                    m1_int += 1;
                    m1_mutations[m1_int] = clone_ind1;
                }
            }

            // if we just hit t2 in last time step, record values
            // this condition will be satisfied during driver = 0 loop
            if (Time >= t2 && t2_yet == 0)
            {   
                t2_yet = 1;
                cells[driver] -= 1.0; // equivalent to cells[0] -= 1.0;

                // get number of passengers present in founder type-2 cell
                // pick which type-0 gets the driver mutation
                rand=gsl_ran_flat(r,0.0,cells[driver]);
                
                for(j=0; j<=num_clones[driver]; j++)
                {
                    if (rand<Clone[j][driver])
                    {
                        clone_ind2 = j;
                        m2_mutations[m2_int] = clone_ind2;
                        j=num_clones[driver]+1;
                    }
                    else rand-=Clone[j][driver];
                }
                Clone[clone_ind2][driver] -= 1.0;
                // now find how many passengers clone_ind2 has at t2 (m)
                for (; clone_ind2 > 0;)
                {
                    clone_ind2 = ancestor[clone_ind2][driver];
                    m2 += 1.0;
                    m2_int += 1;
                    m2_mutations[m2_int] = clone_ind2;
                }
            }

            // if we just hit t1+t in the last time step record values
            if (Time >= (t1+t) && t_yet[driver] == 0)
            {
                t_yet[driver] = 1;
                WBC1[driver] = cells[driver];
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
                t_yet[driver] = 0;
                if (driver == 0)
                {
                    t1_yet = 0;
                    t2_yet = 0;
                }

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

    // mutation and lineage processing
    for (driver=0; driver<3;driver++)
    {
        // clone[i] is going to include all the cells that have mutation i
        for(i=num_clones[driver]; i>0; i--) //first add clones then remove non-zero ones
        {
            Clone[ancestor[i][driver]][driver]+=Clone[i][driver];
        }
        
        // add cells[1] to subclone size for the m1 mutations
        for (i = 0; i <= m1_int; i++)
        {
            Clone[m1_mutations[i]][0] += cells[1];
        }
        // add cells[2] to subclone size for the m2 mutations
        for (i = 0; i <= m2_int; i++)
        {
            Clone[m2_mutations[i]][0] += cells[2];
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
    if (num_clones[0]<max_clones && num_clones[1]<max_clones && num_clones[2]<max_clones)    
    {
        // calculate M1, M2, and clone frequencies
        M1 = WBC1[0] + WBC1[1] + WBC1[2];
        M2 = cells[0] + cells[1] + cells[2];
        alpha1 = WBC1[1]/M1;
        alpha2 = WBC1[2]/M1;
        beta1 = cells[1]/M2;
        beta2 = cells[2]/M2;

        f1 = 0.01;
        f2 = 0.2;
        nf1 = f1*M2;
        nf2 = f2*M2;

        for (driver=0; driver<3;driver++)
        {
            for (i=0; i<num_surviving_clones[driver]; i++)
            {
                if (driver == 1 && surviving_clone[i][driver] == cells[driver]) m1_clonal += 1.0; 
                if (driver == 2 && surviving_clone[i][driver] == cells[driver]) m2_clonal += 1.0; 
                if (surviving_clone[i][driver]>nf1 && surviving_clone[i][driver]<nf2) subclonal += 1.0;
            }
        }
        m1_obs = m1 + m1_clonal;
        m2_obs = m2 + m2_clonal;
        printf("subclonal = %f bn (%f,%f)\n",subclonal, f1, f2);
        
        // Compute estimates
        // note, WBC variable is equivalent to WBC1_i = CCF_i * M1
        r_hat = (1./delta)*log(cells[0]/WBC1[0]);
        r1_hat = (1./delta)*log(cells[1]/WBC1[1]);
        r2_hat = (1./delta)*log(cells[2]/WBC1[2]);            

        u_top = f1*f2*subclonal;
        u_bottom = (f2-f1)*((1-beta1-beta2)/r_hat + beta1/r1_hat + beta2/r2_hat);
        u_hat = u_top/u_bottom;
        
        m1_hat = m1_obs - u_hat/r1_hat;
        m2_hat = m2_obs - u_hat/r2_hat;
        t1_hat = m1_hat/u_hat;
        t2_hat = m2_hat/u_hat;

        tau1_hat = (1./r1_hat)*log(cells[1]);
        tau2_hat = (1./r2_hat)*log(cells[2]);

        #pragma omp critical (write_estimates)
        {
            fp = fopen("sibling_subclones_simulation_results.txt","a");    
            fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",r_hat, r1_hat, r2_hat, u_hat, t1_hat, t2_hat, tau1_hat,tau2_hat, m1, m2, m1_obs, m2_obs, subclonal, alpha1, alpha2, beta1, beta2);
            fclose(fp);
        }
        // free memory
        delete[] Clone;
        delete[] surviving_clone;
        delete[] where;
        delete[] ancestor;
        delete[] surviving_ancestor;
        
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
    return 1.0;
}


int main(int argc, char *argv[]){
    if (argc < 2){
        printf("didn't supply filename");
        exit(EXIT_FAILURE);
    }
    else{
        
        FILE * pFile;
        FILE * fp;
        FILE * ft;
        char * filename = argv[1];
        int param_count = 15;
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
        
        double b[3],d[3],t1,t2,t,delta,u;
        b[0] = stod(str2[0],NULL);
        b[1] = stod(str2[1],NULL);
        b[2] = stod(str2[2],NULL);
        d[0] = stod(str2[3],NULL);
        d[1] = stod(str2[4],NULL);
        d[2] = stod(str2[5],NULL);
        t1 = stod(str2[6],NULL);
        t2 = stod(str2[7],NULL);
        t = stod(str2[8],NULL);
        delta = stod(str2[9],NULL);
        u = stod(str2[10],NULL);

        int num_threads = strtol(str2[11],NULL,10);
        int runs = strtol(str2[12],NULL,10);
        int max_clones = strtol(str2[13],NULL,10);
        int par0_single_run1 = strtol(str2[14],NULL,10);

        printf("argc:%i\n",argc);
    
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        
        true_params_out[0] = (b[0]-d[0]); // r
        true_params_out[1] = (b[1] - d[1]); // r1
        true_params_out[2] = (b[2] - d[2]); // r2
        true_params_out[3] = u; 
        true_params_out[4] = t1;
        true_params_out[5] = t2;
        true_params_out[6] = t + delta; // tau1
        true_params_out[7] = t + delta; // tau2 (have t1 = t2)

        // save the true parameter values to textfile
        // add a header to the output file
        ft = fopen("true_params.txt","w");    
        fprintf(ft,"r,r1,r2,u,t1,t2,tau1,tau2\n");
        fclose(ft);

        ft = fopen("true_params.txt","a");
        fprintf(ft,"%f,%f,%f,%f,%f,%f,%f,%f\n",true_params_out[0],true_params_out[1],true_params_out[2],true_params_out[3],true_params_out[4],true_params_out[5],true_params_out[6],true_params_out[7]);
        fclose(ft);

        int counter;
        
        // add a header to the output file
        fp = fopen("sibling_subclones_simulation_results.txt","w");    
        fprintf(fp,"r_hat,r1_hat,r2_hat,u_hat,t1_hat,t2_hat,tau1_hat,tau2_hat,m1,m2,m1_obs,m2_obs,gamma,alpha1,alpha2,beta1,beta2\n");
        fclose(fp);

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
                    surviving_runs+=expansion(b, d, t1, t2, t, delta, u, max_clones, fp);
                    
                    #pragma omp critical (status)
                    {
                        completed++;
                        cout << completed << " / " << runs << endl; 
                    }
                }
        } else // otherwise, do a single run
        {
            expansion(b, d, t1, t2, t, delta, u, max_clones, fp);
        }
    }

    return 0;
    
}
