library(dplyr)
library(pracma)
library(MASS)
library(Rmpfr)

# functions for the distributions for the different clonal structures
sample_t1_from_post <- function(m, age, u){
  p <- runif(1, min = 0, max = 1)
  CDF <- function(x){
    out <- (1- gammainc(u*x,1+m)[[2]]/gamma(1+m))/(1 - gammainc(age*u,1+m)[[2]]/gamma(1+m)) - p
    return(out)
  }
  return(uniroot(CDF, lower = 0, upper = age)[[1]])
}

v_gammainc <- Vectorize(gammainc)

marginal_t1 <- function(t1,m1,m2,age,u){
  out <- exp(-t1*u)*(t1*u)^m1*(gamma(1+m2) - v_gammainc((age-t1)*u,1+m2)["uppinc",])/u
  return(out)
}

marginal_t2 <- function(t2,m1,m2,age,u){
  out <- exp(-t2*u)*(t2*u)^m2*(gamma(1+m1) - v_gammainc((age-t2)*u,1+m1)["uppinc",])/u
  return(out)
}

const_fun <- function(m1,m2,age,u){
  out <- integrate(marginal_t1,lower = 0, upper = age, m1 = m1, m2 = m2,age = age, u = u)[[1]]
  return(out)
}

analytic_integrate_amt2_dt3dt2dt1 <- function(t3, t2, m1,m2,m3,age,u){
  out <- exp(-t3*u)*(t3*u)^m3*exp(-t2*u)*(t2*u)^m2*(gamma(1+m1) - v_gammainc((age-t2)*u,1+m1)["uppinc",])/u
  return(out)
}

analytic_integrate_amt3_dt3dt2dt1 <- function(t3, t2, m1,m2,m3,age,u){
  out <- exp(-t3*u)*(t3*u)^m3*exp(-t2*u)*(t2*u)^m2*(gamma(1+m1) - v_gammainc((age-t3)*u,1+m1)["uppinc",])/u
  return(out)
}

analytic_integrate_amt2_dt2dt3dt1 <- function(t2, t3, m1,m2,m3,age,u){
  out <- exp(-t3*u)*(t3*u)^m3*exp(-t2*u)*(t2*u)^m2*(gamma(1+m1) - v_gammainc((age-t2)*u,1+m1)["uppinc",])/u
  return(out)
}

analytic_integrate_amt3_dt2dt3dt1 <- function(t2, t3, m1,m2,m3,age,u){
  out <- exp(-t3*u)*(t3*u)^m3*exp(-t2*u)*(t2*u)^m2*(gamma(1+m1) - v_gammainc((age-t3)*u,1+m1)["uppinc",])/u
  return(out)
}

conditional_t1 <- function(t1,t2,t3,m1,m2,m3,age,u){
  if (t3 >= t2){
    out <- u*exp(-t1*u)*(t1*u)^m1/(gamma(1+m1) - v_gammainc((age-t3)*u,1+m1)["uppinc",]) 
  }else{
    out <- u*exp(-t1*u)*(t1*u)^m1/(gamma(1+m1) - v_gammainc((age-t2)*u,1+m1)["uppinc",])
  }
  return(out)
}

t1_t2_joint <- function(t1,t2,m1,m2,age,u){
  out <- (u*t1)^m1*exp(-u*t1)*(u*t2)^m2*exp(-u*t2)
  return(out)
}

conditional_t2 <- function(t1,t2,m1,m2,age,u){
  out <- exp(-u*t2)*u*(u*t2)^m2/(gamma(1+m2) - v_gammainc((age-t1)*u,1+m2)["uppinc",])
  return(out)
}

sample_t1_t2_from_post <- function(m1, m2, age, u){
  p <- runif(1, min = 0, max = 1)
  
  const <- integrate(marginal_t1,lower = 0, upper = age, m1 = m1, m2 = m2,age = age, u = u)[[1]]
  
  CDF_t1 <- function(x){
    out <- integrate(marginal_t1,lower = 0, upper = x, m1 = m1, m2 = m2,age = age, u = u)[[1]]/const - p
    return(out)
  }
  
  t1_out <- uniroot(CDF_t1, lower = 0, upper = age)[[1]]
  t1_CDF_check <- CDF_t1(age)+p
  p <- runif(1, min = 0, max = 1)
  CDF_t2 <- function(x){
    out <- integrate(conditional_t2,lower = 0, upper = x, t1 = t1_out, m1 = m1, m2 = m2,age = age, u = u)[[1]] - p
    return(out)
  }
  t2_out <- uniroot(CDF_t2, lower = 0, upper = age - t1_out)[[1]]

  t2_CDF_check <- CDF_t2(age - t1_out)+p
  return(list(t1_out,t2_out, t1_CDF_check, t2_CDF_check))
}

sample_ti <- function(m = NA, age = NA, u = NA){
  ti_out <- try(sample_t1_from_post(m, age, u))
  if (class(ti_out) == "try-error"){
    return(NA)
  }else{
    return(ti_out)
  }
}


generate_estimates_pt21 <- function(pt = "21", timepoints = c("1","2"),
                                    add_poi_err = TRUE, r_processing = 'lifetime',
                                    n_draws = 10000){
  clones_nums_my_phylogic_runs <- list("21" = 2)
  inputs_df <- read.csv("data/CLL_inputs.csv", header = TRUE, fill = TRUE, na.strings = c("NA"," ",""))
  #create data frame for saving the parameter values
  clone_cols <- c("alpha1", "beta1","gamma1_obs", "gamma2_obs", "gamma1", "gamma2", "m","r","r1","u","MRCA", "t1","t","t1_in_years")
  
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
  f1 <- inputs_df[inputs_df$pt == pt,"f1"]
  f2 <- inputs_df[inputs_df$pt == pt,"f2"]
  L <- inputs_df[inputs_df$pt == pt,"L"]
  age <- inputs_df[inputs_df$pt == pt,"age"]
  gamma1_obs <- inputs_df[inputs_df$pt == pt,"gamma1_obs"]
  gamma2_obs <- inputs_df[inputs_df$pt == pt,"gamma2_obs"]
  rho1 <- inputs_df[inputs_df$pt == pt,"coverage1"]
  rho2 <- inputs_df[inputs_df$pt == pt,"coverage2"]
  purity1 <- inputs_df[inputs_df$pt == pt,"purity1"]
  purity2 <- inputs_df[inputs_df$pt == pt,"purity2"]

  for (i in 1:nrow(param_df)){
    WBC_df <- data.frame(t = c(tp1,tp2), 
                         logWBC_cl0 = log(5*10^9*c(M1*(1-param_df[i,"alpha1"]),M2*(1-param_df[i,"beta1"]))),
                         logWBC_cl1 = log(5*10^9*c(M1*param_df[i,"alpha1"],M2*param_df[i,"beta1"])))
    curve_fit_cl0 <- lm(logWBC_cl0~t, data = WBC_df)
    curve_fit_cl1 <- lm(logWBC_cl1~t, data = WBC_df)
    param_df[i,"r"] <- coef(curve_fit_cl0)[[2]]
    param_df[i,"r1"] <- coef(curve_fit_cl1)[[2]]
    #when the time points are in units of years after diagnosis this t will be
    
    # calculate lifetime growth rates for the other parameter calculations
    # clone 0
    r_min_tp1 <- log(5*10^9*M1*(1-param_df[i,"alpha1"]))/(age+tp1)
    r_min_tp2 <- log(5*10^9*M2*(1-param_df[i,"beta1"]))/(age+tp2)
    param_df[i,'r_tp1'] <- max(param_df[i,'r'],r_min_tp1)
    param_df[i,'r_tp2'] <- max(param_df[i,'r'],r_min_tp2)
    
    # clone 1
    r1_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha1"]))/(age+tp1)
    r1_min_tp2 <- log(5*10^9*M2*(param_df[i,"beta1"]))/(age+tp2)
    param_df[i,'r1_tp1'] <- max(param_df[i,'r1'],r1_min_tp1)
    param_df[i,'r1_tp2'] <- max(param_df[i,'r1'],r1_min_tp2)
    
    # the time between start of cl1 and diagnosis
    param_df[i,"t"] <- coef(curve_fit_cl1)[[1]]/param_df[i,'r1_tp1']
  }

  m_obs <- inputs_df[inputs_df$pt == pt,"m1"]

  #calculated corrected gamma at each timepoint
  param_df["gamma1"] <- gamma1_obs*((1/f1 - 1/f2)/(purity1*rho1/(2*L) - 1/f2))
  param_df["gamma2"] <- gamma2_obs*((1/f1 - 1/f2)/(purity2*rho2/(2*L) - 1/f2))
  param_df["gamma1_obs"] <- gamma1_obs
  param_df["gamma2_obs"] <- gamma2_obs

  # mutation rate is averaged across timepoints
  param_df["u"] <- ((f1*f2*param_df[["gamma1"]])/((f2-f1)*((1-param_df[["alpha1"]])/param_df[["r_tp1"]] + param_df[["alpha1"]]/param_df[['r1_tp1']]))+
                      (f1*f2*param_df[["gamma2"]])/((f2-f1)*((1-param_df[["beta1"]])/param_df[["r_tp2"]] + param_df[["beta1"]]/param_df[['r1_tp2']])))/2
  param_df["t_BD"] <- (1/param_df[['r1_tp1']])*log(M1*(param_df[["alpha1"]])*5*10^9) - tp1

  m <- m_obs - param_df[["u"]]/param_df[["r_tp1"]]

  t_greater_age_counter <- 0
  for (i in 1:nrow(param_df)){
    if (age - param_df[i,"t"] <= 0){
      t_greater_age_counter <- t_greater_age_counter + 1
      param_df[i,"t1"] <- 0
      param_df[i,"MRCA"] <- 0
    }else{
      param_df[i,"t1"] <- sample_t1_from_post(m[[i]], age - param_df[i,"t"], param_df[i,"u"])
      param_df[i,"MRCA"] <- age - param_df[i,"t"] - param_df[i,"t1"]
    }
    
  }
  
  param_df["t1_in_years"] <- param_df[["MRCA"]] + param_df[["t1"]]
  
  # save confidence intervals
  CI_df <- data.frame(matrix(0, nrow = 7, ncol = 3))
  rownames(CI_df) <- c("r","r1","u","MRCA (years)", "t1 (years)", "t (y0/r)", "t (BD proc)")
  colnames(CI_df) <- c("2.5", "50", "97.5")
  CI_df["t1 (years)",] <- quantile(param_df$t1_in_years, probs = c(0.025, 0.5, 0.975)) %>%  round(digits = 1)
  CI_df["MRCA (years)",] <-  quantile(param_df[["MRCA"]],c(.025,0.5,.975)) %>%  round(digits = 1)
  CI_df["u",] <- quantile(param_df[["u"]],c(.025,0.5,.975)) %>%  round(digits = 3)
  CI_df["r",] <- quantile(param_df[["r"]],c(.025,0.5,.975)) %>%  round(digits = 3)
  CI_df["r1",] <- quantile(param_df[["r1"]],c(.025,0.5,.975)) %>%  round(digits = 3)
  CI_df["t (y0/r)",] <- quantile(param_df[["t"]],c(.025,0.5,.975)) %>%  round(digits = 1)
  CI_df["t (BD proc)",] <- quantile(param_df[["t_BD"]],c(.025,0.5,.975)) %>%  round(digits = 1)
  return(list(CI_df,param_df))
  
}

# function to generate estimates for patient 3
generate_estimates_pt3 <- function(pt = "3", timepoints = c("1","2","3"),
                                              add_poi_err = TRUE,
                                              r_processing = 'lifetime',
                                              n_draws = 10000){

  clones_nums_my_phylogic_runs <- list("3" = c(2))
  inputs_df <- read.csv("data/CLL_inputs.csv", header = TRUE, 
                        fill = TRUE, na.strings = c("NA"," ",""))
  #create data frame for saving the parameter values
  
  clone_cols <- c("alpha1", "beta1","ccf_last1","gamma1_obs", "gamma2_obs","gamma3_obs", "gamma1", "gamma2","gamma3", "m","r","r1","u","MRCA", "t1","t","t1_in_years")
  
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
  f1 <- inputs_df[inputs_df$pt == pt,"f1"]
  f2 <- inputs_df[inputs_df$pt == pt,"f2"]
  L <- inputs_df[inputs_df$pt == pt,"L"]
  age <- inputs_df[inputs_df$pt == pt,"age"]
  gamma1_obs <- inputs_df[inputs_df$pt == pt,"gamma1_obs"]
  gamma2_obs <- inputs_df[inputs_df$pt == pt,"gamma2_obs"]
  gamma3_obs <- inputs_df[inputs_df$pt == pt,"gamma3_obs"]
  rho1 <- inputs_df[inputs_df$pt == pt,"coverage1"]
  rho2 <- inputs_df[inputs_df$pt == pt,"coverage2"]
  rho3 <- inputs_df[inputs_df$pt == pt,"coverage3"]
  purity1 <- inputs_df[inputs_df$pt == pt,"purity1"]
  purity2 <- inputs_df[inputs_df$pt == pt,"purity2"]
  purity3 <- inputs_df[inputs_df$pt == pt,"purity3"]
  
  for (i in 1:nrow(param_df)){
    WBC_df <- data.frame(t = c(tp1,tp2,tp3), 
                         logWBC_cl0 = log(5*10^9*c(M1*(1-param_df[i,"alpha1"]),M2*(1-param_df[i,"beta1"]),M3*(1-param_df[i,"ccf_last1"]))),
                         logWBC_cl1 = log(5*10^9*c(M1*param_df[i,"alpha1"],M2*param_df[i,"beta1"],M3*param_df[i,"ccf_last1"])))
    curve_fit_cl0 <- lm(logWBC_cl0~t, data = WBC_df)
    curve_fit_cl1 <- lm(logWBC_cl1~t, data = WBC_df)
    
    sampled_params_cl0 <- mvrnorm(1, coef(curve_fit_cl0), vcov(curve_fit_cl0))
    sampled_params_cl1 <- mvrnorm(1, coef(curve_fit_cl1), vcov(curve_fit_cl1))
    
    param_df[i,"r"] <- sampled_params_cl0[[2]]
    param_df[i,"r1"] <- sampled_params_cl1[[2]]
    
    # lifetime growth rates to generate clone size
    # clone 0
    r_min_tp1 <- log(5*10^9*M1*(1-param_df[i,"alpha1"]))/(age+tp1)
    r_min_tp2 <- log(5*10^9*M2*(1-param_df[i,"beta1"]))/(age+tp2)
    r_min_tp3 <- log(5*10^9*M3*(1-param_df[i,"ccf_last1"]))/(age+tp3)
    param_df[i,'r_tp1'] <- max(param_df[i,'r'],r_min_tp1)
    param_df[i,'r_tp2'] <- max(param_df[i,'r'],r_min_tp2)
    param_df[i,'r_tp3'] <- max(param_df[i,'r'],r_min_tp3)
    
    # clone 1
    r1_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha1"]))/(age+tp1)
    r1_min_tp2 <- log(5*10^9*M2*(param_df[i,"beta1"]))/(age+tp2)
    r1_min_tp3 <- log(5*10^9*M3*(param_df[i,"ccf_last1"]))/(age+tp3)
    param_df[i,'r1_tp1'] <- max(param_df[i,'r1'],r1_min_tp1)
    param_df[i,'r1_tp2'] <- max(param_df[i,'r1'],r1_min_tp2)
    param_df[i,'r1_tp3'] <- max(param_df[i,'r1'],r1_min_tp3)
    
    param_df[i,"t"] <- sampled_params_cl1[[1]]/param_df[i,'r1_tp1']
    param_df[i,"t_BD"] <- (1/param_df[i,"r1_tp1"])*log(M1*(param_df[i,"alpha1"])*5*10^9) - tp1
  }
  
  m_obs <- inputs_df[inputs_df$pt == pt,"m1"]
  
  #calculated corrected gamma
  param_df["gamma1"] <- gamma1_obs*((1/f1 - 1/f2)/(purity1*rho1/(2*L) - 1/f2))
  param_df["gamma2"] <- gamma2_obs*((1/f1 - 1/f2)/(purity2*rho2/(2*L) - 1/f2))
  param_df["gamma3"] <- gamma3_obs*((1/f1 - 1/f2)/(purity3*rho3/(2*L) - 1/f2))
  param_df["gamma1_obs"] <- gamma1_obs
  param_df["gamma2_obs"] <- gamma2_obs
  param_df["gamma3_obs"] <- gamma3_obs
  
  param_df["u"] <- ((f1*f2*param_df[["r1_tp1"]]*param_df[["r_tp1"]]*param_df[["gamma1"]])/((f2-f1)*((param_df[["alpha1"]])*param_df[["r_tp1"]] + param_df[["r1_tp1"]]*(1 - param_df[["alpha1"]])))+
                      (f1*f2*param_df[["r1_tp2"]]*param_df[["r_tp2"]]*param_df[["gamma2"]])/((f2-f1)*((param_df[["beta1"]])*param_df[["r_tp2"]] + param_df[["r1_tp2"]]*(1 - param_df[["beta1"]])))+
                      (f1*f2*param_df[["r1_tp3"]]*param_df[["r_tp3"]]*param_df[["gamma3"]])/((f2-f1)*((param_df[["ccf_last1"]])*param_df[["r_tp3"]] + param_df[["r1_tp3"]]*(1 - param_df[["ccf_last1"]]))))/3
  m <- m_obs - param_df[["u"]]/param_df[["r1_tp1"]]
  if (sum(param_df["u"] < 0)>0){
    warning("u < 0")
  }
  
  #add in uncertainty from the poisson process for the time t1
  t_greater_age_counter <- 0
  for (i in 1:nrow(param_df)){
    if (age - param_df[i,"t"] <= 0){
      t_greater_age_counter <- t_greater_age_counter + 1
      param_df[i,"t1"] <- 0
      param_df[i,"MRCA"] <- 0
    }else{
      param_df[i,"t1"] <- sample_t1_from_post(m[[i]], age - param_df[i,"t"], param_df[i,"u"])
      param_df[i,"MRCA"] <- age - param_df[i,"t"] - param_df[i,"t1"]
    }
    
  }
  
  param_df["t1_in_years"] <- param_df[["MRCA"]] + param_df[["t1"]]
  
  
  row_names <- c("r","r1","u","MRCA (years)", "t1 (years)","t (y0/r)","t (BD proc)")
  CI_df <- data.frame(matrix(0, nrow = length(row_names), ncol = 3))
  rownames(CI_df) <- row_names
  colnames(CI_df) <- c("2.5", "50", "97.5")
  CI_df["t1 (years)",] <- quantile(param_df$t1_in_years, probs = c(0.025, 0.5, 0.975)) %>%  round(digits = 1)
  CI_df["MRCA (years)",] <-  quantile(param_df[["MRCA"]],c(.025,0.5,.975))  %>%  round(digits = 1)
  CI_df["u",] <- quantile(param_df[["u"]],c(.025,0.5,.975)) %>%  round(digits = 3)
  CI_df["r",] <- quantile(param_df[["r"]],c(.025,0.5,.975)) %>%  round(digits = 3)
  CI_df["r1",] <- quantile(param_df[["r1"]],c(.025,0.5,.975)) %>%  round(digits = 3)
  CI_df["t (y0/r)",] <- quantile(param_df[["t"]],c(.025,0.5,.975)) %>%  round(digits = 1)
  CI_df["t (BD proc)",] <- quantile(param_df[["t_BD"]],c(.025,0.5,.975)) %>%  round(digits = 1)

  return(list(CI_df,param_df))
  }
  

generate_estimates_pt6 <- function(pt = "6", timepoints = c("1","2","3"),
                                                  add_poi_err = TRUE,
                                                  r_processing = 'lifetime',
                                                  n_draws = 10000){
  clones_nums_my_phylogic_runs <- list("6"=c(2,3,4))
  inputs_df <- read.csv("data/CLL_inputs.csv", header = TRUE, fill = TRUE, na.strings = c("NA"," ",""))
  #create data frame for saving the parameter values
  
  clone_cols <- c("alpha1", "alpha2", "alpha3","beta1", "beta2","beta3","ccf_last1","ccf_last2", "ccf_last3","gamma1_obs", "gamma2_obs","gamma3_obs", "gamma1", "gamma2","gamma3","m1", "m2","m3","r","r1","r2","r3","u","MRCA", "t1","t2","t3","t")
  
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
  f1 <- inputs_df[inputs_df$pt == pt,"f1"]
  f2 <- inputs_df[inputs_df$pt == pt,"f2"]
  L <- inputs_df[inputs_df$pt == pt,"L"]
  age <- inputs_df[inputs_df$pt == pt,"age"]
  gamma1_obs <- inputs_df[inputs_df$pt == pt,"gamma1_obs"]
  gamma2_obs <- inputs_df[inputs_df$pt == pt,"gamma2_obs"]
  gamma3_obs <- inputs_df[inputs_df$pt == pt,"gamma3_obs"]
  param_df["gamma1_obs"] <- gamma1_obs
  param_df["gamma2_obs"] <- gamma2_obs
  param_df["gamma3_obs"] <- gamma3_obs
  rho1 <- inputs_df[inputs_df$pt == pt,"coverage1"]
  rho2 <- inputs_df[inputs_df$pt == pt,"coverage2"]
  rho3 <- inputs_df[inputs_df$pt == pt,"coverage3"]
  purity1 <- inputs_df[inputs_df$pt == pt,"purity1"]
  purity2 <- inputs_df[inputs_df$pt == pt,"purity2"]
  purity3 <- inputs_df[inputs_df$pt == pt,"purity3"]
  
  m1_obs <- inputs_df[inputs_df$pt == pt,"m1"]
  m2_obs <- inputs_df[inputs_df$pt == pt,"m2"]
  m3_obs <- inputs_df[inputs_df$pt == pt,"m3"]
  
  #calculated corrected gamma at each timepoint
  param_df["gamma1"] <- gamma1_obs*((1/f1 - 1/f2)/(purity1*rho1/(2*L) - 1/f2))
  param_df["gamma2"] <- gamma2_obs*((1/f1 - 1/f2)/(purity2*rho2/(2*L) - 1/f2))
  param_df["gamma3"] <- gamma3_obs*((1/f1 - 1/f2)/(purity3*rho3/(2*L) - 1/f2))
  
  # fit WBC counts over multiple timepoints
  t_after_driver1 <- vector(length = nrow(param_df))
  for (i in 1:nrow(param_df)){
    WBC_df <- data.frame(t = c(tp1,tp2,tp3), 
                         logWBC_cl0 = log(5*10^9*c(M1*(1-param_df[i,"alpha1"] - param_df[i,"alpha2"]),M2*(1-param_df[i,"beta1"] - param_df[i,"beta2"]),M3*(1-param_df[i,"ccf_last1"] - param_df[i,"ccf_last2"]))),
                         logWBC_cl1 = log(5*10^9*c(M1*(param_df[i,"alpha1"]-param_df[i,"alpha3"]),M2*(param_df[i,"beta1"]-param_df[i,"beta3"]),M3*(param_df[i,"ccf_last1"]-param_df[i,"ccf_last3"]))),
                         logWBC_cl2 = log(5*10^9*c(M1*(param_df[i,"alpha2"]),M2*(param_df[i,"beta2"]),M3*(param_df[i,"ccf_last2"]))),
                         logWBC_cl3 = log(5*10^9*c(M1*param_df[i,"alpha3"],M2*param_df[i,"beta3"],M3*param_df[i,"ccf_last3"])))
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
    r_min_tp2 <- log(5*10^9*M2*(1-param_df[i,"beta1"]-param_df[i,"beta2"]))/(age+tp2)
    r_min_tp3 <- log(5*10^9*M3*(1-param_df[i,"ccf_last1"]-param_df[i,"ccf_last2"]))/(age+tp3)
    param_df[i,'r_tp1'] <- max(param_df[i,'r'],r_min_tp1)
    param_df[i,'r_tp2'] <- max(param_df[i,'r'],r_min_tp2)
    param_df[i,'r_tp3'] <- max(param_df[i,'r'],r_min_tp3)
    
    # clone 1
    r1_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha1"]-param_df[i,"alpha3"]))/(age+tp1)
    r1_min_tp2 <- log(5*10^9*M2*(param_df[i,"beta1"]-param_df[i,"beta3"]))/(age+tp2)
    r1_min_tp3 <- log(5*10^9*M3*(param_df[i,"ccf_last1"]-param_df[i,"ccf_last3"]))/(age+tp3)
    param_df[i,'r1_tp1'] <- max(param_df[i,'r1'],r1_min_tp1)
    param_df[i,'r1_tp2'] <- max(param_df[i,'r1'],r1_min_tp2)
    param_df[i,'r1_tp3'] <- max(param_df[i,'r1'],r1_min_tp3)
    
    # clone 2
    r2_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha2"]))/(age+tp1)
    r2_min_tp2 <- log(5*10^9*M2*(param_df[i,"beta2"]))/(age+tp2)
    r2_min_tp3 <- log(5*10^9*M3*(param_df[i,"ccf_last2"]))/(age+tp3)
    param_df[i,'r2_tp1'] <- max(param_df[i,'r2'],r2_min_tp1,na.rm = T)
    param_df[i,'r2_tp2'] <- max(param_df[i,'r2'],r2_min_tp2,na.rm = T)
    param_df[i,'r2_tp3'] <- max(param_df[i,'r2'],r2_min_tp3,na.rm = T)
    
    # clone 3 
    r3_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha3"]))/(age+tp1)
    r3_min_tp2 <- log(5*10^9*M2*(param_df[i,"beta3"]))/(age+tp2)
    r3_min_tp3 <- log(5*10^9*M3*(param_df[i,"ccf_last3"]))/(age+tp3)
    param_df[i,'r3_tp1'] <- max(param_df[i,'r3'],r3_min_tp1,na.rm = T)
    param_df[i,'r3_tp2'] <- max(param_df[i,'r3'],r3_min_tp2,na.rm = T)
    param_df[i,'r3_tp3'] <- max(param_df[i,'r3'],r3_min_tp3,na.rm = T)
    
    param_df[i,"t"] <- sampled_params_cl3[[1]]/param_df[i,'r3_tp1']
    
    t_after_driver1[[i]] <- sampled_params_cl1[[1]]/sampled_params_cl1[[2]] 
    
  }
  
  # calculate mutation rate, averaged across timepoints
  param_df["u"] <- ((f1*f2*param_df[["gamma1"]])/((f2-f1)*((1-param_df[["alpha1"]] - param_df[["alpha2"]])/param_df[["r_tp1"]] + (param_df[["alpha1"]]-param_df[["alpha3"]])/param_df[["r1_tp1"]] + (param_df[["alpha2"]])/param_df[["r2_tp1"]]+param_df[["alpha3"]]/param_df[["r3_tp1"]]))+
                      (f1*f2*param_df[["gamma2"]])/((f2-f1)*((1-param_df[["beta1"]] - param_df[["beta2"]])/param_df[["r_tp2"]] + (param_df[["beta1"]]-param_df[["beta3"]])/param_df[["r1_tp2"]] + (param_df[["beta2"]])/param_df[["r2_tp2"]] + param_df[["beta3"]]/param_df[["r3_tp2"]])) + 
                      (f1*f2*param_df[["gamma3"]])/((f2-f1)*((1-param_df[["ccf_last1"]] - param_df[["ccf_last2"]])/param_df[["r_tp3"]] + (param_df[["ccf_last1"]]-param_df[["ccf_last3"]])/param_df[["r1_tp3"]] + (param_df[["ccf_last2"]])/param_df[["r2_tp3"]]+param_df[["ccf_last3"]]/param_df[["r3_tp3"]])))/3
  m1 <- m1_obs - param_df[["u"]]/param_df[["r1_tp1"]]
  m2 <- m2_obs - param_df[["u"]]/param_df[["r2_tp1"]]
  m3 <- m3_obs - param_df[["u"]]/param_df[["r3_tp1"]]+param_df[["u"]]/param_df[["r1_tp1"]]
  
  # compare the intercept method for t with the approximation using branching processes
  param_df["t_BD"] <- (1/param_df[["r3_tp1"]])*log(M1*(param_df[["alpha3"]])*5*10^9) - tp1 
  
  
  param_df["m1"] <- m1
  param_df["m2"] <- m2
  param_df["m3"] <- m3
  
  if (sum(param_df[["u"]] < 0)>0){
    warning("u < 0 for:")
    print(param_df[param_df[["u"]]< 0,])
  }
  
  t_greater_age_counter <- 0
  for (i in 1:nrow(param_df)){
    if (age <= param_df[i,"t"]){
      t_greater_age_counter <- t_greater_age_counter + 1
      param_df[i,"MRCA"] <- 0
      param_df[i,"t1"] <- 0
      param_df[i,"t1_in_years"] <- 0
      param_df[i,"t3"] <- 0
      param_df[i,"t3_in_years"] <- 0
      param_df[i,"t2"] <- sample_ti(m2[[i]], age, param_df[i,"u"])
      param_df[i,"t2_in_years"] <- param_df[i,"t2"]
    }else{
      t_out <- sample_t1_t2_from_post(m1[[i]], m3[[i]], age-param_df[i,"t"], param_df[i,"u"])
      param_df[i,"t1"] <- t_out[[1]]
      param_df[i,"t3"] <- t_out[[2]]
      param_df[i,"MRCA"] <- age - param_df[i,"t3"]- param_df[i,"t1"]  - param_df[i,"t"]
      param_df[i,"t1_in_years"] <- param_df[i,"MRCA"] + param_df[i,"t1"]
      param_df[i,"t3_in_years"] <- param_df[i,"MRCA"] + param_df[i,"t1"]+ param_df[i,"t3"]
      param_df[i,"t2"] <- sample_ti(m2[[i]], param_df[i,"t1"]+param_df[i,"t3"]+param_df[i,"t"], param_df[i,"u"])
      param_df[i,"t2_in_years"] <- param_df[i,"MRCA"] + param_df[i,"t2"]
    }
  }
  row_names <- c("r","r1","r2","r3","u","MRCA (years)", 
                 "t1 (years)", "t2 (years)","t3 (years)","t (y0/r)","t (BD proc)")
  CI_df <- data.frame(matrix(0, nrow = length(row_names), ncol = 3))
  rownames(CI_df) <- row_names
  colnames(CI_df) <- c("2.5", "50", "97.5")

  CI_df["t1 (years)",] <- quantile(param_df$t1_in_years, probs = c(0.025, 0.5, 0.975), na.rm = T) %>%  round(digits = 1)
  CI_df["t2 (years)",] <- quantile(param_df$t2_in_years, probs = c(0.025, 0.5, 0.975), na.rm = T) %>%  round(digits = 1)
  CI_df["t3 (years)",] <- quantile(param_df$t3_in_years, probs = c(0.025, 0.5, 0.975), na.rm = T) %>%  round(digits = 1)
  CI_df["MRCA (years)",] <-  quantile(param_df[["MRCA"]],c(.025,0.5,.975), na.rm = T)  %>%  round(digits = 1)
  CI_df["u",] <- quantile(param_df[["u"]],c(.025,0.5,.975), na.rm = T) %>%  round(digits = 3)
  CI_df["r",] <- quantile(param_df[["r"]],c(.025,0.5,.975), na.rm = T) %>%  round(digits = 3)
  CI_df["r1",] <- quantile(param_df[["r1"]],c(.025,0.5,.975), na.rm = T) %>%  round(digits = 3)
  CI_df["r2",] <- quantile(param_df[["r2"]],c(.025,0.5,.975), na.rm = T) %>%  round(digits = 3)
  CI_df["r3",] <- quantile(param_df[["r3"]],c(.025,0.5,.975), na.rm = T) %>%  round(digits = 3)
  CI_df["t (y0/r)",] <- quantile(param_df[["t"]],c(.025,0.5,.975)) %>%  round(digits = 1)
  CI_df["t (BD proc)",] <- quantile(param_df[["t_BD"]],c(.025,0.5,.975)) %>%  round(digits = 1)
  CI_df["t1_na_count",] <- sum(is.na(param_df$t1_in_years))
  CI_df["t2_na_count",] <- sum(is.na(param_df$t2_in_years))
  CI_df["t3_na_count",] <- sum(is.na(param_df$t3_in_years))
  
  return(list(CI_df,param_df))
}

# function to generate estimates for patient 9
generate_estimates_pt9 <- function(pt = "9", timepoints = c("1","2","3"),
                                              add_poi_err = TRUE,
                                              r_processing = 'lifetime',
                                              n_draws = 10000){
  clones_nums_my_phylogic_runs <- list("9"=c(3,2,4))
  inputs_df <- read.csv("data/CLL_inputs.csv", header = TRUE, fill = TRUE, na.strings = c("NA"," ",""))
  #create data frame for saving the parameter values
  
  clone_cols <- c("alpha1", "alpha2", "alpha3","beta1", "beta2","beta3","ccf_last1","ccf_last2", "ccf_last3","gamma1_obs", "gamma2_obs","gamma3_obs", "gamma1", "gamma2","gamma3","m1", "m2","m3","r","r1","r2","r3","u","MRCA", "t1","t2","t3","t")
  
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
  f1 <- inputs_df[inputs_df$pt == pt,"f1"]
  f2 <- inputs_df[inputs_df$pt == pt,"f2"]
  L <- inputs_df[inputs_df$pt == pt,"L"]
  age <- inputs_df[inputs_df$pt == pt,"age"]
  gamma1_obs <- inputs_df[inputs_df$pt == pt,"gamma1_obs"]
  gamma2_obs <- inputs_df[inputs_df$pt == pt,"gamma2_obs"]
  gamma3_obs <- inputs_df[inputs_df$pt == pt,"gamma3_obs"]
  param_df["gamma1_obs"] <- gamma1_obs
  param_df["gamma2_obs"] <- gamma2_obs
  param_df["gamma3_obs"] <- gamma3_obs
  rho1 <- inputs_df[inputs_df$pt == pt,"coverage1"]
  rho2 <- inputs_df[inputs_df$pt == pt,"coverage2"]
  rho3 <- inputs_df[inputs_df$pt == pt,"coverage3"]
  purity1 <- inputs_df[inputs_df$pt == pt,"purity1"]
  purity2 <- inputs_df[inputs_df$pt == pt,"purity2"]
  purity3 <- inputs_df[inputs_df$pt == pt,"purity3"]
  
  m1_obs <- inputs_df[inputs_df$pt == pt,"m1"]
  m2_obs <- inputs_df[inputs_df$pt == pt,"m2"]
  m3_obs <- inputs_df[inputs_df$pt == pt,"m3"]
  
  #calculated corrected gamma for each timepoint
  param_df["gamma1"] <- gamma1_obs*((1/f1 - 1/f2)/(purity1*rho1/(2*L) - 1/f2))
  param_df["gamma2"] <- gamma2_obs*((1/f1 - 1/f2)/(purity2*rho2/(2*L) - 1/f2))
  param_df["gamma3"] <- gamma3_obs*((1/f1 - 1/f2)/(purity3*rho3/(2*L) - 1/f2))
  
  # fit WBC count across the 3 timepoints
  t_after_driver1 <- vector(length = nrow(param_df))
  for (i in 1:nrow(param_df)){
    WBC_df <- data.frame(t = c(tp1,tp2,tp3), 
                         logWBC_cl0 = log(5*10^9*c(M1*(1-param_df[i,"alpha1"] - param_df[i,"alpha2"]),M2*(1-param_df[i,"beta1"] - param_df[i,"beta2"]),M3*(1-param_df[i,"ccf_last1"] - param_df[i,"ccf_last2"]))),
                         logWBC_cl1 = log(5*10^9*c(M1*param_df[i,"alpha1"],M2*param_df[i,"beta1"],M3*param_df[i,"ccf_last1"])),
                         logWBC_cl2 = log(5*10^9*c(M1*(param_df[i,"alpha2"]-param_df[i,"alpha3"]),M2*(param_df[i,"beta2"]-param_df[i,"beta3"]),M3*(param_df[i,"ccf_last2"]-param_df[i,"ccf_last3"]))),
                         logWBC_cl3 = log(5*10^9*c(M1*param_df[i,"alpha3"],M2*param_df[i,"beta3"],M3*param_df[i,"ccf_last3"])))
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
    r_min_tp2 <- log(5*10^9*M2*(1-param_df[i,"beta1"]-param_df[i,"beta2"]))/(age+tp2)
    r_min_tp3 <- log(5*10^9*M3*(1-param_df[i,"ccf_last1"]-param_df[i,"ccf_last2"]))/(age+tp3)
    param_df[i,'r_tp1'] <- max(param_df[i,'r'],r_min_tp1)
    param_df[i,'r_tp2'] <- max(param_df[i,'r'],r_min_tp2)
    param_df[i,'r_tp3'] <- max(param_df[i,'r'],r_min_tp3)
    
    # clone 1
    r1_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha1"]))/(age+tp1)
    r1_min_tp2 <- log(5*10^9*M2*(param_df[i,"beta1"]))/(age+tp2)
    r1_min_tp3 <- log(5*10^9*M3*(param_df[i,"ccf_last1"]))/(age+tp3)
    param_df[i,'r1_tp1'] <- max(param_df[i,'r1'],r1_min_tp1)
    param_df[i,'r1_tp2'] <- max(param_df[i,'r1'],r1_min_tp2)
    param_df[i,'r1_tp3'] <- max(param_df[i,'r1'],r1_min_tp3)
    
    # clone 2
    r2_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha2"]-param_df[i,"alpha3"]))/(age+tp1)
    r2_min_tp2 <- log(5*10^9*M2*(param_df[i,"beta2"]-param_df[i,"beta3"]))/(age+tp2)
    r2_min_tp3 <- log(5*10^9*M3*(param_df[i,"ccf_last2"]-param_df[i,"ccf_last3"]))/(age+tp3)
    param_df[i,'r2_tp1'] <- max(param_df[i,'r2'],r2_min_tp1,na.rm = T)
    param_df[i,'r2_tp2'] <- max(param_df[i,'r2'],r2_min_tp2,na.rm = T)
    param_df[i,'r2_tp3'] <- max(param_df[i,'r2'],r2_min_tp3,na.rm = T)
    
    # clone 3 
    r3_min_tp1 <- log(5*10^9*M1*(param_df[i,"alpha3"]))/(age+tp1)
    r3_min_tp2 <- log(5*10^9*M2*(param_df[i,"beta3"]))/(age+tp2)
    r3_min_tp3 <- log(5*10^9*M3*(param_df[i,"ccf_last3"]))/(age+tp3)
    param_df[i,'r3_tp1'] <- max(param_df[i,'r3'],r3_min_tp1,na.rm = T)
    param_df[i,'r3_tp2'] <- max(param_df[i,'r3'],r3_min_tp2,na.rm = T)
    param_df[i,'r3_tp3'] <- max(param_df[i,'r3'],r3_min_tp3,na.rm = T)
    
    param_df[i,"t"] <- sampled_params_cl2[[1]]/param_df[i,'r2_tp1']
    
    t_after_driver1[[i]] <- sampled_params_cl1[[1]]/sampled_params_cl1[[2]] 
    
  }
  
  # calculate mutation rate, averaged across 3 timepoints
  param_df["u"] <- ((f1*f2*param_df[["gamma1"]])/((f2-f1)*((1-param_df[["alpha1"]] - param_df[["alpha2"]])/param_df[["r_tp1"]] + (param_df[["alpha1"]])/param_df[["r1_tp1"]] + (param_df[["alpha2"]]-param_df[["alpha3"]])/param_df[["r2_tp1"]]+param_df[["alpha3"]]/param_df[["r3_tp1"]]))+
                      (f1*f2*param_df[["gamma2"]])/((f2-f1)*((1-param_df[["beta1"]] - param_df[["beta2"]])/param_df[["r_tp2"]] + (param_df[["beta1"]])/param_df[["r1_tp2"]] + (param_df[["beta2"]]-param_df[["beta3"]])/param_df[["r2_tp2"]] + param_df[["beta3"]]/param_df[["r3_tp2"]])) + 
                      (f1*f2*param_df[["gamma3"]])/((f2-f1)*((1-param_df[["ccf_last1"]] - param_df[["ccf_last2"]])/param_df[["r_tp3"]] + (param_df[["ccf_last1"]])/param_df[["r1_tp3"]] + (param_df[["ccf_last2"]]-param_df[["ccf_last3"]])/param_df[["r2_tp3"]]+param_df[["ccf_last3"]]/param_df[["r3_tp3"]])))/3
  m1 <- m1_obs - param_df[["u"]]/param_df[["r1_tp1"]]
  m2 <- m2_obs - param_df[["u"]]/param_df[["r2_tp1"]]
  m3 <- m3_obs - param_df[["u"]]/param_df[["r3_tp1"]]+param_df[["u"]]/param_df[["r2_tp1"]]
  
  # compare the intercept method for t with the approximation using branching processes
  param_df["t_BD"] <- (1/param_df[["r2_tp1"]])*log(M1*(param_df[["alpha2"]]-param_df[["alpha3"]])*5*10^9) - tp1 
  
  
  param_df["m1"] <- m1
  param_df["m2"] <- m2
  param_df["m3"] <- m3
  
  if (sum(param_df[["u"]] < 0)>0){
    warning("u < 0 for:")
    print(param_df[param_df[["u"]]< 0,])
  }
  
  t_greater_age_counter <- 0
  for (i in 1:nrow(param_df)){
    if (age <= param_df[i,"t"]){
      t_greater_age_counter <- t_greater_age_counter + 1
      param_df[i,"MRCA"] <- 0
      param_df[i,"t2"] <- 0
      param_df[i,"t2_in_years"] <- 0
      param_df[i,"t1"] <- sample_ti(m1[[i]], age, param_df[i,"u"])
      param_df[i,"t3"] <- sample_ti(m3[[i]], age, param_df[i,"u"])
      param_df[i,"t3_in_years"] <- param_df[i,"MRCA"] + param_df[i,"t2"]+ param_df[i,"t3"]
    }else{
      t_out <- sample_t1_t2_from_post(m2[[i]], m3[[i]], age-param_df[i,"t"], param_df[i,"u"])
      param_df[i,"t2"] <- t_out[[1]]
      param_df[i,"t3"] <- t_out[[2]]
      param_df[i,"MRCA"] <- age - param_df[i,"t2"] - param_df[i,"t"]
      param_df[i,"t2_in_years"] <- param_df[i,"MRCA"] + param_df[i,"t2"]
      param_df[i,"t3_in_years"] <- param_df[i,"MRCA"] + param_df[i,"t2"]+ param_df[i,"t3"]
      param_df[i,"t1"] <- sample_ti(m1[[i]], param_df[i,"t2"]+param_df[i,"t"], param_df[i,"u"])
    }
    
    
    param_df[i,"t1_in_years"] <- param_df[i,"MRCA"] + param_df[i,"t1"]
  }
  row_names <- c("r","r1","r2","r3","u","MRCA (years)", 
                 "t1 (years)", "t2 (years)","t3 (years)","t (y0/r)","t (BD proc)")
  CI_df <- data.frame(matrix(0, nrow = length(row_names), ncol = 3))
  rownames(CI_df) <- row_names
  colnames(CI_df) <- c("2.5", "50", "97.5")

  CI_df["t1 (years)",] <- quantile(param_df$t1_in_years, probs = c(0.025, 0.5, 0.975), na.rm = T) %>%  round(digits = 1)
  CI_df["t2 (years)",] <- quantile(param_df$t2_in_years, probs = c(0.025, 0.5, 0.975), na.rm = T) %>%  round(digits = 1)
  CI_df["t3 (years)",] <- quantile(param_df$t3_in_years, probs = c(0.025, 0.5, 0.975), na.rm = T) %>%  round(digits = 1)
  CI_df["MRCA (years)",] <-  quantile(param_df[["MRCA"]],c(.025,0.5,.975), na.rm = T)  %>%  round(digits = 1)
  CI_df["u",] <- quantile(param_df[["u"]],c(.025,0.5,.975), na.rm = T) %>%  round(digits = 3)
  CI_df["r",] <- quantile(param_df[["r"]],c(.025,0.5,.975), na.rm = T) %>%  round(digits = 3)
  CI_df["r1",] <- quantile(param_df[["r1"]],c(.025,0.5,.975), na.rm = T) %>%  round(digits = 3)
  CI_df["r2",] <- quantile(param_df[["r2"]],c(.025,0.5,.975), na.rm = T) %>%  round(digits = 3)
  CI_df["r3",] <- quantile(param_df[["r3"]],c(.025,0.5,.975), na.rm = T) %>%  round(digits = 3)
  CI_df["t (y0/r)",] <- quantile(param_df[["t"]],c(.025,0.5,.975)) %>%  round(digits = 1)
  CI_df["t (BD proc)",] <- quantile(param_df[["t_BD"]],c(.025,0.5,.975)) %>%  round(digits = 1)
  CI_df["t1_na_count",] <- sum(is.na(param_df$t1_in_years))
  CI_df["t2_na_count",] <- sum(is.na(param_df$t2_in_years))
  CI_df["t3_na_count",] <- sum(is.na(param_df$t3_in_years))
  
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

    filename_CIs <- paste0("CIs_final/pt",pt,"_CIs.csv")
    filename_params <- paste0("CIs_final/pt",pt,"_params.csv")
    write.csv(out[[1]], file = filename_CIs)
    write.csv(out[[2]], file = filename_params)
    
  }
}
  

