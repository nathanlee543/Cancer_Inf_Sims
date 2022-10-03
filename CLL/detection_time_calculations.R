library(dplyr)
library(MASS)

WBC_t <- function(r,y0,WBC,x){
  if (sign(sum((exp(r*0+y0))) - WBC) == sign(sum((exp(r*100+y0))) - WBC)){
    print(y0)
    print(r)
    print(exp(r*0+y0))
  }
  out <- sum((exp(r*x+y0))) - WBC
  return(out)  
}  

generate_estimates_pt21 <- function(pt = "21", timepoints = c("1","2"),
                                    n_draws = 10^4){
  clones_nums_my_phylogic_runs <- list("21" = 2)
  inputs_df <- read.csv("data/CLL_inputs.csv", header = TRUE, fill = TRUE, na.strings = c("NA"," ",""))
  #create data frame for saving the parameter values
  clone_cols <- c("alpha1", "beta1","r","r1","t_cll_subclones","t_M1_subclones","t_M2_subclones")
  
  param_df <- data.frame(data = matrix(0,nrow = n_draws, ncol = length(clone_cols)))
  colnames(param_df) <- clone_cols
  
  # draw CCFs
  CCF_post_df <- read.table("data/pt21.cluster_ccfs.txt",
                            header = TRUE, fill = TRUE, na.strings = c("NA"," ",""), sep = "\t")
  CCF_vals <- seq(from = 0, to = 1, by = 0.01)
  
  
  CCF_roots <- c("alpha", "beta")
  for (i in 1:length(timepoints)){
    for (j in 1:length(clones_nums_my_phylogic_runs[[pt]]))
    {
      t <- paste0("sample_0",timepoints[i])
      id <- clones_nums_my_phylogic_runs[[pt]][j]
      
      CCF_post <- filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id))[(ncol(filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id)))-100):ncol(filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id)))]
      CCF <- sample(CCF_vals, size = n_draws, replace = TRUE, prob = CCF_post)
      
      # will be alpha_1, alpha2, beta1, or beta2, depending on patient
      colname <- paste0(CCF_roots[i],j)
      param_df[colname] <- CCF
    }
  }
  
  
  # compute estimates
  delta <- inputs_df[inputs_df$pt == pt,"delta"]
  M1 <- inputs_df[inputs_df$pt == pt,"M1"]
  M2 <- inputs_df[inputs_df$pt == pt,"M2"]
  tp1 <- inputs_df[inputs_df$pt == pt,"tp1"]
  tp2 <- inputs_df[inputs_df$pt == pt,"tp2"]

  age <- inputs_df[inputs_df$pt == pt,"age"]
  times <- c(tp1,tp2) + age
  
  for (i in 1:nrow(param_df)){
    WBC_df <- data.frame(t = times, 
                         logWBC_cl0 = log(5*10^9*c(M1*(1-param_df[i,"alpha1"]),M2*(1-param_df[i,"beta1"]))),
                         logWBC_cl1 = log(5*10^9*c(M1*param_df[i,"alpha1"],M2*param_df[i,"beta1"])),
                         logWBC = log(5*10^9*c(M1,M2)))
    curve_fit_cl0 <- lm(logWBC_cl0~t, data = WBC_df)
    curve_fit_cl1 <- lm(logWBC_cl1~t, data = WBC_df)
    param_df[i,"r"] <- coef(curve_fit_cl0)[[2]]
    param_df[i,"r1"] <- coef(curve_fit_cl1)[[2]]

    # calculate lifetime growth rates for the other parameter calculations
    # clone 0
    r_min_tp1 <- log(5*10^9*M1*(1-param_df[i,"alpha1"]))/(age+tp1)
    param_df[i,'r_tp1'] <- max(param_df[i,'r'],r_min_tp1)
    
    lifetime_r <- param_df[i,'r']<r_min_tp1
    param_df[i,"lifetime_r"] <- lifetime_r
    
    # clone 1
    r1_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha1"]))/(age+tp1)
    r1_min_tp2 <- log(5*10^9*M2*(param_df[i,"beta1"]))/(age+tp2)
    param_df[i,'r1_tp1'] <- max(param_df[i,'r1'],r1_min_tp1)
    lifetime_r1 <- param_df[i,'r1']<r1_min_tp1
    param_df[i,"lifetime_r1"] <- lifetime_r1
    
    y0 <- vector(length = 2)
    if (lifetime_r){
      y0[1] <- 1
    }else{
      y0[1] <- coef(curve_fit_cl0)[[1]]
    }
    
    if (lifetime_r1){
      y0[2] <- 1
    }else{
      y0[2] <- coef(curve_fit_cl1)[[1]]
    }
    
    r <- c(param_df[i,'r_tp1'],param_df[i,'r1_tp1'] )


    param_df[i,"t_cll_subclones"] <- uniroot(function(x) WBC_t(r,y0,5.75*10^10,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M1_subclones"] <- uniroot(function(x) WBC_t(r,y0,M1*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M2_subclones"] <- uniroot(function(x) WBC_t(r,y0,M2*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M1_subclones_PE"] <- 100*(param_df[i,"t_M1_subclones"]-(times[1]))/times[1]
    param_df[i,"t_M2_subclones_PE"] <- 100*(param_df[i,"t_M2_subclones"]-times[2])/times[2]

  }
  
  row_names <- c("t_cll_subclones","t_M1_subclones","t_M2_subclones","t_M1_subclones_PE","t_M2_subclones_PE")
  # save confidence intervals
  CI_df <- data.frame(matrix(0, nrow = length(row_names), ncol = 3))
  rownames(CI_df) <- row_names
  colnames(CI_df) <- c("2.5", "50", "97.5")
  
  for (param in rownames(CI_df)){
    CI_df[param,] <- quantile(param_df[[param]],c(.025,0.5,.975)) %>%  round(digits = 3)
  }
  
  return(list(CI_df,param_df))
  
}

# function to generate estimates for patient 3
generate_estimates_pt3 <- function(pt = "3", timepoints = c("1","2","3"),
                                   n_draws = 10000){
  
  clones_nums_my_phylogic_runs <- list("3" = c(2))
  inputs_df <- read.csv("data/CLL_inputs.csv", header = TRUE, 
                        fill = TRUE, na.strings = c("NA"," ",""))
  #create data frame for saving the parameter values
  
  clone_cols <- c("alpha1", "beta1","ccf_last1","r","r1","u","t_cll_subclones","t_M1_subclones","t_M2_subclones","t_M3_subclones")
  
  param_df <- data.frame(data = matrix(0,nrow = n_draws, ncol = length(clone_cols)))
  colnames(param_df) <- clone_cols
  
  # draw CCFs
  CCF_post_df <- read.table("data/pt3.cluster_ccfs.txt",
                            header = TRUE, fill = TRUE, na.strings = c("NA"," ",""), sep = "\t")
  CCF_vals <- seq(from = 0, to = 1, by = 0.01)
  
  CCF_roots <- c("alpha", "beta","ccf_last")
  for (i in 1:length(timepoints)){
    for (j in 1:length(clones_nums_my_phylogic_runs[[pt]]))
    {
      t <- paste0("sample_0",timepoints[i])
      id <- clones_nums_my_phylogic_runs[[pt]][j]
      
      CCF_post <- filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id))[(ncol(filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id)))-100):ncol(filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id)))]
      CCF <- sample(CCF_vals, size = n_draws, replace = TRUE, prob = CCF_post)
      
      # will be alpha_1, alpha2, beta1, or beta2, depending on patient
      colname <- paste0(CCF_roots[i],j)
      param_df[colname] <- CCF
    }
  }
  
  # compute estimates
  delta <- inputs_df[inputs_df$pt == pt,"delta"]
  M1 <- inputs_df[inputs_df$pt == pt,"M1"]
  M2 <- inputs_df[inputs_df$pt == pt,"M2"]
  M3 <- inputs_df[inputs_df$pt == pt,"M3"]
  tp1 <- inputs_df[inputs_df$pt == pt,"tp1"]
  tp2 <- inputs_df[inputs_df$pt == pt,"tp2"]
  tp3 <- inputs_df[inputs_df$pt == pt,"tp3"]
  age <- inputs_df[inputs_df$pt == pt,"age"]
  times <- c(tp1,tp2,tp3) + age

  for (i in 1:nrow(param_df)){
    WBC_df <- data.frame(t = times, 
                         logWBC_cl0 = log(5*10^9*c(M1*(1-param_df[i,"alpha1"]),M2*(1-param_df[i,"beta1"]),M3*(1-param_df[i,"ccf_last1"]))),
                         logWBC_cl1 = log(5*10^9*c(M1*param_df[i,"alpha1"],M2*param_df[i,"beta1"],M3*param_df[i,"ccf_last1"])),
                         logWBC = log(5*10^9*c(M1,M2,M3)))
    curve_fit_cl0 <- lm(logWBC_cl0~t, data = WBC_df)
    curve_fit_cl1 <- lm(logWBC_cl1~t, data = WBC_df)
    
    sampled_params_cl0 <- mvrnorm(1, coef(curve_fit_cl0), vcov(curve_fit_cl0))
    sampled_params_cl1 <- mvrnorm(1, coef(curve_fit_cl1), vcov(curve_fit_cl1))
    
    param_df[i,"r"] <- sampled_params_cl0[[2]]
    param_df[i,"r1"] <- sampled_params_cl1[[2]]
    
    # lifetime growth rates to generate clone size
    # clone 0
    r_min_tp1 <- log(5*10^9*M1*(1-param_df[i,"alpha1"]))/(age+tp1)
    param_df[i,'r_tp1'] <- max(param_df[i,'r'],r_min_tp1)
    lifetime_r <- param_df[i,'r']<r_min_tp1
    param_df[i,"lifetime_r"] <- lifetime_r
    
    # clone 1
    r1_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha1"]))/(age+tp1)
    param_df[i,'r1_tp1'] <- max(param_df[i,'r1'],r1_min_tp1)
    lifetime_r1 <- param_df[i,'r1']<r1_min_tp1
    param_df[i,"lifetime_r1"] <- lifetime_r1

    y0 <- vector(length = 2)
    if (lifetime_r){
      y0[1] <- 1
    }else{
      y0[1] <- sampled_params_cl0[[1]]
    }
    
    if (lifetime_r1){
      y0[2] <- 1
    }else{
      y0[2] <- sampled_params_cl1[[1]]
    }
    
    r <- c(param_df[i,'r_tp1'],param_df[i,'r1_tp1'] )


    param_df[i,"t_cll_subclones"] <- uniroot(function(x) WBC_t(r,y0,5.75*10^10,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M1_subclones"] <- uniroot(function(x) WBC_t(r,y0,M1*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M2_subclones"] <- uniroot(function(x) WBC_t(r,y0,M2*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M3_subclones"] <- uniroot(function(x) WBC_t(r,y0,M3*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M1_subclones_PE"] <- 100*(param_df[i,"t_M1_subclones"]-(times[1]))/times[1]
    param_df[i,"t_M2_subclones_PE"] <- 100*(param_df[i,"t_M2_subclones"]-times[2])/times[2]
    param_df[i,"t_M3_subclones_PE"] <- 100*(param_df[i,"t_M3_subclones"]-times[3])/times[3]
  }
  
  row_names <- c("t_cll_subclones","t_M1_subclones","t_M2_subclones","t_M3_subclones","t_M1_subclones_PE","t_M2_subclones_PE","t_M3_subclones_PE")
  CI_df <- data.frame(matrix(0, nrow = length(row_names), ncol = 3))
  rownames(CI_df) <- row_names
  colnames(CI_df) <- c("2.5", "50", "97.5")
  for (param in rownames(CI_df)){
    CI_df[param,] <- quantile(param_df[[param]],c(.025,0.5,.975)) %>%  round(digits = 3)
  }

  return(list(CI_df,param_df))
}


generate_estimates_pt6 <- function(pt = "6", timepoints = c("1","2","3"),
                                   n_draws = 10000){
  clones_nums_my_phylogic_runs <- list("6"=c(2,3,4))
  inputs_df <- read.csv("data/CLL_inputs.csv", header = TRUE, fill = TRUE, na.strings = c("NA"," ",""))
  #create data frame for saving the parameter values
  
  clone_cols <- c("alpha1", "alpha2", "alpha3","beta1", "beta2","beta3","ccf_last1","ccf_last2", "ccf_last3","r","r1","r2","r3","t_cll_subclones","t_M1_subclones","t_M2_subclones","t_M3_subclones")
  
  param_df <- data.frame(data = matrix(0,nrow = n_draws, ncol = length(clone_cols)))
  colnames(param_df) <- clone_cols
  
  # draw CCFs from phylogicNDT posteriors
  CCF_post_df <- read.table("data/pt6.cluster_ccfs.txt",
                            header = TRUE, fill = TRUE, na.strings = c("NA"," ",""), sep = "\t")
  CCF_vals <- seq(from = 0, to = 1, by = 0.01)
  CCF_draw_fun <- function(support = NA, timepoint = NA, cluster = NA, CCF_post_df = NA){
    t <- paste0("sample_0",timepoint)
    id <- cluster
    
    CCF_post <- filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id))[(ncol(filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id)))-100):ncol(filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id)))]
    CCF <- sample(support, size = 1, replace = TRUE, prob = CCF_post)
    
    return(CCF)
  }
  for (i in 1:n_draws){
    need_another_loop <- 1
    while (need_another_loop == 1) {
      param_df[i,"alpha1"] <- CCF_draw_fun(support = CCF_vals, 
                                           timepoint = 1, cluster = 2, 
                                           CCF_post_df = CCF_post_df)
      
      param_df[i,"alpha2"] <- CCF_draw_fun(support = CCF_vals, 
                                           timepoint = 1, cluster = 3, 
                                           CCF_post_df = CCF_post_df)
      
      param_df[i,"alpha3"] <- CCF_draw_fun(support = CCF_vals, 
                                           timepoint = 1, cluster = 4, 
                                           CCF_post_df = CCF_post_df)
      param_df[i,"beta1"] <- CCF_draw_fun(support = CCF_vals, 
                                          timepoint = 2, cluster = 2, 
                                          CCF_post_df = CCF_post_df)
      
      param_df[i,"beta2"] <- CCF_draw_fun(support = CCF_vals, 
                                          timepoint = 2, cluster = 3, 
                                          CCF_post_df = CCF_post_df)
      
      param_df[i,"beta3"] <- CCF_draw_fun(support = CCF_vals, 
                                          timepoint = 2, cluster = 4, 
                                          CCF_post_df = CCF_post_df)
      
      param_df[i,"ccf_last1"] <- CCF_draw_fun(support = CCF_vals, 
                                              timepoint = 3, cluster = 2, 
                                              CCF_post_df = CCF_post_df)
      
      param_df[i,"ccf_last2"] <- CCF_draw_fun(support = CCF_vals, 
                                              timepoint = 3, cluster = 3, 
                                              CCF_post_df = CCF_post_df)
      
      param_df[i,"ccf_last3"] <- CCF_draw_fun(support = CCF_vals, 
                                              timepoint = 3, cluster = 4, 
                                              CCF_post_df = CCF_post_df)
      
      if ((param_df[i,"alpha3"] >= param_df[i,"alpha1"])| 
          (param_df[i,"beta3"] >=param_df[i,"beta1"])|
          (param_df[i,"ccf_last3"] >= param_df[i,"ccf_last1"])|
          (1-param_df[i,"alpha1"] - param_df[i,"alpha2"] <= 0)|
          (1-param_df[i,"beta1"] - param_df[i,"beta2"] <= 0)|
          (1-param_df[i,"ccf_last1"] - param_df[i,"ccf_last2"] <= 0)
      ){
        need_another_loop <- 1
      }else{
        need_another_loop <- 0
      }
    }
  }
  
  
  # compute estimates
  delta <- inputs_df[inputs_df$pt == pt,"delta"]
  M1 <- inputs_df[inputs_df$pt == pt,"M1"]
  M2 <- inputs_df[inputs_df$pt == pt,"M2"]
  M3 <- inputs_df[inputs_df$pt == pt,"M3"]
  tp1 <- inputs_df[inputs_df$pt == pt,"tp1"]
  tp2 <- inputs_df[inputs_df$pt == pt,"tp2"]
  tp3 <- inputs_df[inputs_df$pt == pt,"tp3"]
  age <- inputs_df[inputs_df$pt == pt,"age"]
  
  times <- c(tp1,tp2,tp3) + age
  for (i in 1:nrow(param_df)){
    WBC_df <- data.frame(t = times, 
                         logWBC_cl0 = log(5*10^9*c(M1*(1-param_df[i,"alpha1"] - param_df[i,"alpha2"]),M2*(1-param_df[i,"beta1"] - param_df[i,"beta2"]),M3*(1-param_df[i,"ccf_last1"] - param_df[i,"ccf_last2"]))),
                         logWBC_cl1 = log(5*10^9*c(M1*(param_df[i,"alpha1"]-param_df[i,"alpha3"]),M2*(param_df[i,"beta1"]-param_df[i,"beta3"]),M3*(param_df[i,"ccf_last1"]-param_df[i,"ccf_last3"]))),
                         logWBC_cl2 = log(5*10^9*c(M1*(param_df[i,"alpha2"]),M2*(param_df[i,"beta2"]),M3*(param_df[i,"ccf_last2"]))),
                         logWBC_cl3 = log(5*10^9*c(M1*param_df[i,"alpha3"],M2*param_df[i,"beta3"],M3*param_df[i,"ccf_last3"])),
                         logWBC = log(5*10^9*c(M1,M2,M3)))
    curve_fit_cl0 <- lm(logWBC_cl0~t, data = WBC_df)
    curve_fit_cl1 <- lm(logWBC_cl1~t, data = WBC_df)
    curve_fit_cl2 <- lm(logWBC_cl2~t, data = WBC_df)
    curve_fit_cl3 <- lm(logWBC_cl3~t, data = WBC_df)
    
    # sample from mv normal of fit parameters
    sampled_params_cl0 <- mvrnorm(1, coef(curve_fit_cl0), vcov(curve_fit_cl0))
    sampled_params_cl1 <- mvrnorm(1, coef(curve_fit_cl1), vcov(curve_fit_cl1))
    sampled_params_cl2 <- mvrnorm(1, coef(curve_fit_cl2), vcov(curve_fit_cl2))
    sampled_params_cl3 <- mvrnorm(1, coef(curve_fit_cl3), vcov(curve_fit_cl3))
    
    param_df[i,"r"] <- sampled_params_cl0[[2]]
    param_df[i,"r1"] <- sampled_params_cl1[[2]]
    param_df[i,"r2"] <- sampled_params_cl2[[2]]
    param_df[i,"r3"] <- sampled_params_cl3[[2]]
    
    # clone 0
    r_min_tp1 <- log(5*10^9*M1*(1-param_df[i,"alpha1"]-param_df[i,"alpha2"]))/(age+tp1)
    param_df[i,'r_tp1'] <- max(param_df[i,'r'],r_min_tp1)
    lifetime_r <- param_df[i,'r']<r_min_tp1
    param_df[i,"lifetime_r"] <- lifetime_r

    # clone 1
    r1_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha1"]-param_df[i,"alpha3"]))/(age+tp1)
    param_df[i,'r1_tp1'] <- max(param_df[i,'r1'],r1_min_tp1)
    lifetime_r1 <- param_df[i,'r1']<r1_min_tp1
    param_df[i,"lifetime_r1"] <- lifetime_r1

    # clone 2
    r2_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha2"]))/(age+tp1)
    param_df[i,'r2_tp1'] <- max(param_df[i,'r2'],r2_min_tp1,na.rm = T)
    lifetime_r2 <- param_df[i,'r2']<r2_min_tp1
    param_df[i,"lifetime_r2"] <- lifetime_r2    
    
    # clone 3 
    r3_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha3"]))/(age+tp1)
    param_df[i,'r3_tp1'] <- max(param_df[i,'r3'],r3_min_tp1,na.rm = T)
    lifetime_r3 <- param_df[i,'r3']<r3_min_tp1
    param_df[i,"lifetime_r3"] <- lifetime_r3    

    y0 <- vector(length = (length(timepoints)+1))
    
    if (lifetime_r){
      y0[1] <- 1
    }else{
      y0[1] <- sampled_params_cl0[[1]]
    }
    
    if (lifetime_r1){
      y0[2] <- 1
    }else{
      y0[2] <- sampled_params_cl1[[1]]
    }

    if (lifetime_r2){
      y0[3] <- 1
    }else{
      y0[3] <- sampled_params_cl2[[1]]
    }

    if (lifetime_r3){
      y0[4] <- 1
    }else{
      y0[4] <- sampled_params_cl3[[1]]
    }
    
    r <- c(param_df[i,'r_tp1'],param_df[i,'r1_tp1'],param_df[i,'r2_tp1'],param_df[i,'r3_tp1'] )
    
    param_df[i,"t_cll_subclones"] <- uniroot(function(x) WBC_t(r,y0,5.75*10^10,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M1_subclones"] <- uniroot(function(x) WBC_t(r,y0,M1*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M2_subclones"] <- uniroot(function(x) WBC_t(r,y0,M2*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M3_subclones"] <- uniroot(function(x) WBC_t(r,y0,M3*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M1_subclones_PE"] <- 100*(param_df[i,"t_M1_subclones"]-(times[1]))/times[1]
    param_df[i,"t_M2_subclones_PE"] <- 100*(param_df[i,"t_M2_subclones"]-times[2])/times[2]
    param_df[i,"t_M3_subclones_PE"] <- 100*(param_df[i,"t_M3_subclones"]-times[3])/times[3]
  }
  
  row_names <- c("t_cll_subclones","t_M1_subclones","t_M2_subclones","t_M3_subclones","t_M1_subclones_PE","t_M2_subclones_PE","t_M3_subclones_PE")
  CI_df <- data.frame(matrix(0, nrow = length(row_names), ncol = 3))
  rownames(CI_df) <- row_names
  colnames(CI_df) <- c("2.5", "50", "97.5")
  for (param in rownames(CI_df)){
    CI_df[param,] <- quantile(param_df[[param]],c(.025,0.5,.975)) %>%  round(digits = 3)
  }
  
  return(list(CI_df,param_df))
}

# function to generate estimates for patient 9
generate_estimates_pt9 <- function(pt = "9", timepoints = c("1","2","3"),
                                   n_draws = 10000){
  clones_nums_my_phylogic_runs <- list("9"=c(3,2,4))
  inputs_df <- read.csv("data/CLL_inputs.csv", header = TRUE, fill = TRUE, na.strings = c("NA"," ",""))
  #create data frame for saving the parameter values
  
  clone_cols <- c("alpha1", "alpha2", "alpha3","beta1", "beta2","beta3","ccf_last1","ccf_last2", "ccf_last3","r","r1","r2","r3","t_cll_subclones","t_M1_subclones","t_M2_subclones","t_M3_subclones")
  
  param_df <- data.frame(data = matrix(0,nrow = n_draws, ncol = length(clone_cols)))
  colnames(param_df) <- clone_cols
  
  # draw CCFs
  CCF_post_df <- read.table("data/pt9.cluster_ccfs.txt",
                            header = TRUE, fill = TRUE, na.strings = c("NA"," ",""), sep = "\t")
  CCF_vals <- seq(from = 0, to = 1, by = 0.01)
  CCF_draw_fun <- function(support = NA, timepoint = NA, cluster = NA, CCF_post_df = NA){
    t <- paste0("sample_0",timepoint)
    id <- cluster
    
    CCF_post <- filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id))[(ncol(filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id)))-100):ncol(filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id)))]
    CCF <- sample(support, size = 1, replace = TRUE, prob = CCF_post)
    
    return(CCF)
  }
  for (i in 1:n_draws){
    need_another_loop <- 1
    while (need_another_loop == 1) {
      param_df[i,"alpha1"] <- CCF_draw_fun(support = CCF_vals, 
                                           timepoint = 1, cluster = 3, 
                                           CCF_post_df = CCF_post_df)
      
      param_df[i,"alpha2"] <- CCF_draw_fun(support = CCF_vals, 
                                           timepoint = 1, cluster = 2, 
                                           CCF_post_df = CCF_post_df)
      
      param_df[i,"alpha3"] <- CCF_draw_fun(support = CCF_vals, 
                                           timepoint = 1, cluster = 4, 
                                           CCF_post_df = CCF_post_df)
      param_df[i,"beta1"] <- CCF_draw_fun(support = CCF_vals, 
                                          timepoint = 2, cluster = 3, 
                                          CCF_post_df = CCF_post_df)
      
      param_df[i,"beta2"] <- CCF_draw_fun(support = CCF_vals, 
                                          timepoint = 2, cluster = 2, 
                                          CCF_post_df = CCF_post_df)
      
      param_df[i,"beta3"] <- CCF_draw_fun(support = CCF_vals, 
                                          timepoint = 2, cluster = 4, 
                                          CCF_post_df = CCF_post_df)
      
      param_df[i,"ccf_last1"] <- CCF_draw_fun(support = CCF_vals, 
                                              timepoint = 3, cluster = 3, 
                                              CCF_post_df = CCF_post_df)
      
      param_df[i,"ccf_last2"] <- CCF_draw_fun(support = CCF_vals, 
                                              timepoint = 3, cluster = 2, 
                                              CCF_post_df = CCF_post_df)
      
      param_df[i,"ccf_last3"] <- CCF_draw_fun(support = CCF_vals, 
                                              timepoint = 3, cluster = 4, 
                                              CCF_post_df = CCF_post_df)
      
      if (param_df[i,"alpha3"] >= param_df[i,"alpha2"] | 
          param_df[i,"beta3"] >= param_df[i,"beta2"]|
          param_df[i,"ccf_last3"] >= param_df[i,"ccf_last2"]|
          ((1-param_df[i,"alpha1"]-param_df[i,"alpha2"])<=0)|
          ((1-param_df[i,"beta1"]-param_df[i,"beta2"])<=0)|
          ((1-param_df[i,"ccf_last1"]-param_df[i,"ccf_last2"])<=0)
      ){
        need_another_loop <- 1
      }else{
        need_another_loop <- 0
      }
    }
  }
  
  
  # compute estimates
  delta <- inputs_df[inputs_df$pt == pt,"delta"]
  M1 <- inputs_df[inputs_df$pt == pt,"M1"]
  M2 <- inputs_df[inputs_df$pt == pt,"M2"]
  M3 <- inputs_df[inputs_df$pt == pt,"M3"]
  tp1 <- inputs_df[inputs_df$pt == pt,"tp1"]
  tp2 <- inputs_df[inputs_df$pt == pt,"tp2"]
  tp3 <- inputs_df[inputs_df$pt == pt,"tp3"]
  
  age <- inputs_df[inputs_df$pt == pt,"age"]
  
  times <- c(tp1,tp2,tp3) + age
  for (i in 1:nrow(param_df)){
    WBC_df <- data.frame(t = times, 
                         logWBC_cl0 = log(5*10^9*c(M1*(1-param_df[i,"alpha1"] - param_df[i,"alpha2"]),M2*(1-param_df[i,"beta1"] - param_df[i,"beta2"]),M3*(1-param_df[i,"ccf_last1"] - param_df[i,"ccf_last2"]))),
                         logWBC_cl1 = log(5*10^9*c(M1*param_df[i,"alpha1"],M2*param_df[i,"beta1"],M3*param_df[i,"ccf_last1"])),
                         logWBC_cl2 = log(5*10^9*c(M1*(param_df[i,"alpha2"]-param_df[i,"alpha3"]),M2*(param_df[i,"beta2"]-param_df[i,"beta3"]),M3*(param_df[i,"ccf_last2"]-param_df[i,"ccf_last3"]))),
                         logWBC_cl3 = log(5*10^9*c(M1*param_df[i,"alpha3"],M2*param_df[i,"beta3"],M3*param_df[i,"ccf_last3"])),
                         logWBC = log(5*10^9*c(M1,M2,M3)))
    curve_fit_cl0 <- lm(logWBC_cl0~t, data = WBC_df)
    curve_fit_cl1 <- lm(logWBC_cl1~t, data = WBC_df)
    curve_fit_cl2 <- lm(logWBC_cl2~t, data = WBC_df)
    curve_fit_cl3 <- lm(logWBC_cl3~t, data = WBC_df)
    
    # sample from parameter of the fit
    sampled_params_cl0 <- mvrnorm(1, coef(curve_fit_cl0), vcov(curve_fit_cl0))
    sampled_params_cl1 <- mvrnorm(1, coef(curve_fit_cl1), vcov(curve_fit_cl1))
    sampled_params_cl2 <- mvrnorm(1, coef(curve_fit_cl2), vcov(curve_fit_cl2))
    sampled_params_cl3 <- mvrnorm(1, coef(curve_fit_cl3), vcov(curve_fit_cl3))
    
    param_df[i,"r"] <- sampled_params_cl0[[2]]
    param_df[i,"r1"] <- sampled_params_cl1[[2]]
    param_df[i,"r2"] <- sampled_params_cl2[[2]]
    param_df[i,"r3"] <- sampled_params_cl3[[2]]
    
    # clone 0
    r_min_tp1 <- log(5*10^9*M1*(1-param_df[i,"alpha1"]-param_df[i,"alpha2"]))/(age+tp1)
    param_df[i,'r_tp1'] <- max(param_df[i,'r'],r_min_tp1)
    lifetime_r <- param_df[i,'r']<r_min_tp1
    param_df[i,"lifetime_r"] <- lifetime_r

    # clone 1
    r1_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha1"]))/(age+tp1)
    param_df[i,'r1_tp1'] <- max(param_df[i,'r1'],r1_min_tp1)
    lifetime_r1 <- param_df[i,'r1']<r1_min_tp1
    param_df[i,"lifetime_r1"] <- lifetime_r1
    
    # clone 2
    r2_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha2"]-param_df[i,"alpha3"]))/(age+tp1)
    param_df[i,'r2_tp1'] <- max(param_df[i,'r2'],r2_min_tp1,na.rm = T)
    lifetime_r2 <- param_df[i,'r2']<r2_min_tp1
    param_df[i,"lifetime_r2"] <- lifetime_r2    
    
    # clone 3 
    r3_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha3"]))/(age+tp1)
    param_df[i,'r3_tp1'] <- max(param_df[i,'r3'],r3_min_tp1,na.rm = T)
    lifetime_r3 <- param_df[i,'r3']<r3_min_tp1
    param_df[i,"lifetime_r3"] <- lifetime_r3 
    y0 <- vector(length = (length(timepoints)+1))
    
    if (lifetime_r){
      y0[1] <- 1
    }else{
      y0[1] <- sampled_params_cl0[[1]]
    }
    
    if (lifetime_r1){
      y0[2] <- 1
    }else{
      y0[2] <- sampled_params_cl1[[1]]
    }

    if (lifetime_r2){
      y0[3] <- 1
    }else{
      y0[3] <- sampled_params_cl2[[1]]
    }

    if (lifetime_r3){
      y0[4] <- 1
    }else{
      y0[4] <- sampled_params_cl3[[1]]
    }
    
    r <- c(param_df[i,'r_tp1'],param_df[i,'r1_tp1'],param_df[i,'r2_tp1'],param_df[i,'r3_tp1'] )

    param_df[i,"t_cll_subclones"] <- uniroot(function(x) WBC_t(r,y0,5.75*10^10,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M1_subclones"] <- uniroot(function(x) WBC_t(r,y0,M1*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M2_subclones"] <- uniroot(function(x) WBC_t(r,y0,M2*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M3_subclones"] <- uniroot(function(x) WBC_t(r,y0,M3*5*10^9,x), lower = 0, upper = 100)[[1]]
    param_df[i,"t_M1_subclones_PE"] <- 100*(param_df[i,"t_M1_subclones"]-(times[1]))/times[1]
    param_df[i,"t_M2_subclones_PE"] <- 100*(param_df[i,"t_M2_subclones"]-times[2])/times[2]
    param_df[i,"t_M3_subclones_PE"] <- 100*(param_df[i,"t_M3_subclones"]-times[3])/times[3]
  }
  
  row_names <- c("t_cll_subclones","t_M1_subclones","t_M2_subclones","t_M3_subclones","t_M1_subclones_PE","t_M2_subclones_PE","t_M3_subclones_PE")
  CI_df <- data.frame(matrix(0, nrow = length(row_names), ncol = 3))
  rownames(CI_df) <- row_names
  colnames(CI_df) <- c("2.5", "50", "97.5")
  
  for (param in rownames(CI_df)){
    CI_df[param,] <- quantile(param_df[[param]],c(.025,0.5,.975)) %>%  round(digits = 3)
  }
  
  return(list(CI_df,param_df))
}

# function for generating and saving the posterior distributions of multiple patients
generate_CIs_multiple_pts <- function(pt_list){
  for (pt in pt_list){
    print(paste0("generating confidence intervals for pt.", pt))
    if (pt == '3'){
      out <- generate_estimates_pt3()
    }else if (pt == '6'){
      out <- generate_estimates_pt6()
    }else if (pt == '9'){
      out <- generate_estimates_pt9()
    }else if (pt == '21'){
      out <- generate_estimates_pt21()
    }
    
    filename_CIs <- paste0("data/detection_time/pt",pt,"_CIs.csv")
    filename_params <- paste0("data/detection_time/pt",pt,"_params.csv")
    write.csv(out[[1]], file = filename_CIs)
    write.csv(out[[2]], file = filename_params)
    
  }
}

generate_CIs_multiple_pts(list("3","21","6","9"))

