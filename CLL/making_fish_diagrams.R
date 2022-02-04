library(ggplot2)
library(dplyr)
library(tidyr)


# Making the fish diagrams ========================

CCF_mean <- function(timepoint = NA, cluster = NA, CCF_post_df = NA){
  t <- paste0("sample_0",timepoint)
  id <- cluster
  CCF <- filter(CCF_post_df, (Sample_ID == t & Cluster_ID == id))$postDP_ccf_mean
  return(CCF)
}

make_fish_diagram <- function(pt = ""){
  # clone IDs as used for analysis
  clone_ids_list <- list("3" = c(0,1),
                         "6" = c(0,1,2,3),
                         "9" = c(0,1,2,3),
                         "21" = c(0,1))
  clone_ids <- clone_ids_list[[pt]]
  
  # clone IDs conistent with numbering from Gruber et al paper, used in PhylogicNDT
  phylogic_clone_ids_list <- list("3" = c(2),
                                  "6" = c(2,3,4),
                                  "9"=c(3,2,4),
                                  "21" = c(2))
  phylogic_clone_ids <- phylogic_clone_ids_list[[pt]]
  clones <- paste0("cl.",clone_ids)
  
  #set ordering of clone segments for plotting purposes
  clone_plotting_order_list <- list("3" = c("cl.0","cl.1","cl.0-2"),
                                    "6" = c("cl.0","cl.1","cl.3","cl.1-2","cl.0-2","cl.2","cl.0-3"),
                                    "9" = c("cl.0","cl.1","cl.0-2","cl.2","cl.3","cl.2-2","cl.0-3"),
                                    "21" =  c("cl.0","cl.1","cl.0-2"))
  clone_plotting_order <- clone_plotting_order_list[[pt]]
  
  # number of sections to split up clone into for plotting purposes
  num_clone_splits_list <- list("3" = c(2,1),
                                    "6" = c(3,2,1,1),
                                    "9" = c(3,1,2,1),
                                    "21" =  c(2,1))
  num_clone_splits <- num_clone_splits_list[[pt]]
  
  # set colors of all clones, including for split up clone segments
  clone_colors_list <- list("3" = c("dodgerblue3", "firebrick3","dodgerblue3"),
                            "6" = c("dodgerblue3", "firebrick3","gray80","firebrick3","dodgerblue3","goldenrod2","dodgerblue3"),
                            "9" = c("dodgerblue3", "firebrick3","dodgerblue3", "goldenrod2","gray80","goldenrod2","dodgerblue3"),
                            "21" = c("dodgerblue3", "firebrick3","dodgerblue3"))
  clone_colors <- clone_colors_list[[pt]]
  
  # number of sequencing timepoints
  n_list <- list("3" = 3,
                 "6" = 3,
                 "9" = 3,
                 "21" = 2)
  n <- n_list[[pt]]
  
  # create vector of column names for time of events
  t0_vec <- c("MRCA (years)",paste0(paste0("t",1:n)," (years)"))
  
  CCF_filename <- paste0("data/pt",pt,".cluster_ccfs.txt")
  CCF_post_df <- read.table(CCF_filename,
                            header = TRUE, fill = TRUE, na.strings = c("NA"," ",""), sep = "\t")
  
  CI_filename <- paste0("data/CIs_final/pt",pt,"_CIs.csv")
  CIs <- read.csv(CI_filename, header = TRUE, fill = TRUE, na.strings = c("NA"," ",""), row.names = 1)
  inputs_df <- read.csv("data/CLL_inputs.csv", header = TRUE, fill = TRUE, na.strings = c("NA"," ",""))
  
  age <- inputs_df[inputs_df$pt == pt,"age"]
  seq_times <- c()
  M_j <- c()
  for (j in 1:n){
    assign(paste0("M",j),inputs_df[inputs_df$pt == pt,paste0("M",j)])
    assign(paste0("tp",j),inputs_df[inputs_df$pt == pt,paste0("tp",j)])
    seq_times[j] <- age + inputs_df[inputs_df$pt == pt,paste0("tp",j)]
    M_j[j] <- inputs_df[inputs_df$pt == pt,paste0("M",j)]
  }
  
  ccf_stems <- c("alpha","beta","ccf_last")
  ccf_stems_pt <- ccf_stems[1:n]
  
  for (j in 1:n){
    for (cl in 1:length(phylogic_clone_ids)){
      ccf_mean <- CCF_mean(timepoint = j, 
                           cluster = phylogic_clone_ids[cl], 
                           CCF_post_df = CCF_post_df)
      assign(paste0(ccf_stems_pt[j],cl),ccf_mean)
    }
  }
  
  if (pt == "3"){
    WBC_df <- 5*10^9*data.frame(cl.0 = c(M1*(1-alpha1),M2*(1-beta1),M3*(1-ccf_last1)),
                         cl.1 = c(M1*alpha1,M2*beta1,M3*ccf_last1))
  }else if (pt == "6"){
    WBC_df <- 5*10^9*data.frame(cl.0 = c(M1*(1-alpha1 - alpha2),M2*(1-beta1 - beta2),M3*(1-ccf_last1 - ccf_last2)),
                         cl.1 = c(M1*(alpha1-alpha3),M2*(beta1-beta3),M3*(ccf_last1-ccf_last3)),
                         cl.2 = c(M1*(alpha2),M2*(beta2),M3*(ccf_last2)),
                         cl.3 = c(M1*alpha3,M2*beta3,M3*ccf_last3))
  }else if (pt == "9"){
    WBC_df <- 5*10^9*data.frame(cl.0 = c(M1*(1-alpha1 - alpha2),M2*(1-beta1 - beta2),M3*(1-ccf_last1 - ccf_last2)),
                                cl.1 = c(M1*alpha1,M2*beta1,M3*ccf_last1),
                                cl.2 = c(M1*(alpha2-alpha3),M2*(beta2-beta3),M3*(ccf_last2-ccf_last3)),
                                cl.3 = c(M1*alpha3,M2*beta3,M3*ccf_last3))
  }else if (pt == "21"){
    WBC_df <- 5*10^9*data.frame(cl.0 = c(M1*(1-alpha1),M2*(1-beta1)),
                         cl.1 = c(M1*alpha1,M2*beta1))
  }
  
  pre_dx_times <- seq(from = CIs["MRCA (years)","X50"], to = age, by = .1)
  all_times <- c(pre_dx_times,seq_times)
  col.names <- c("time", clones)
  fish_df <- data.frame(matrix(1,nrow = length(all_times), ncol = length(col.names)))
  colnames(fish_df) <- col.names
  fish_df["time"] <- all_times
  
  zeros_df <- data.frame(matrix(1,nrow = length(all_times), ncol = ncol(WBC_df)))
  colnames(zeros_df) <- colnames(WBC_df)

  # initialize
  for (i in 1:length(clones)){
    t0 <- CIs[t0_vec[i],"X50"]
    for (j in 1:nrow(zeros_df)){
      if (all_times[j] < t0){
        zeros_df[j,i] <- 0
      }
    }
  }
  
  # setting how the clones will be nested
  for (i in 1:length(clones)){
    t0 <- CIs[t0_vec[i],"X50"]
    ri <- log(WBC_df[1,clones[i]])/(seq_times[1] - t0)
    
    fish_df[1:length(pre_dx_times),clones[i]] <- log10(exp(ri*(fish_df[1:length(pre_dx_times),"time"]-t0)))/num_clone_splits[i]
    fish_df[(length(pre_dx_times)+1):length(all_times),clones[i]] <- log10(WBC_df[[clones[i]]])/num_clone_splits[i]
    fish_df[,clones[i]] <- fish_df[,clones[i]]*zeros_df[,clones[i]]
    if (num_clone_splits[i] > 1){
      for (k in 2:num_clone_splits[i])
      fish_df[paste0(clones[i],"-",k)] <- fish_df[,clones[i]]
    }
    
  }
  # reshape dataframe for plotting
  fish_df <- fish_df %>% pivot_longer(colnames(fish_df)[2:length(colnames(fish_df))], names_to = "clone", values_to = "WBC")
  fish_df$factor <- factor(fish_df$clone, levels = clone_plotting_order)
  
  lwd <- 1.25
  arrow_length <- .2*max(c(rowSums(log10(WBC_df))))
  
  p <- ggplot(fish_df)+
    geom_area(aes(x = time, y = WBC, fill = factor))+
    scale_fill_manual(values=clone_colors, breaks = clones)+
    geom_linerange(data = data.frame(x1 = age,y1 = 0,
                                     lower = 0,
                                     upper = max(c(rowSums(log10(WBC_df))))),
                   aes(x = x1, y = y1,ymin = lower, ymax = upper, fill = NULL),
                   linetype = "longdash",
                   color = "white",
                   lwd = lwd)+  
    geom_segment(data = data.frame(x1 = seq_times,y1 = c(rowSums(log10(WBC_df)))+arrow_length,
                                   x2 = seq_times,
                                   y2 =  c(rowSums(log10(WBC_df)))),
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 lineend = "round",
                 linejoin = "bevel",
                 size = 1,
                 arrow = arrow(length = unit(0.25, "cm")))+
    scale_y_continuous(limits = c(0, NA), expand = c(0,0))+
    scale_x_continuous(breaks = ceiling(seq(min(fish_df$time), max(fish_df$time), by = 10)), limits =c(2.8,69.78))+
    theme_bw()+
    labs(fill="Clone")+
    xlab("patient age (years)")+
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          panel.border = element_blank(),
          axis.text.x = element_text(family="Arial", size = 20),
          axis.title.x = element_text(family="Arial", size = 25),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length.x =unit(.25, "cm"),
          aspect.ratio = 0.3)
  filename_out <- paste0("pt",pt,"_fish_diagram.png")
  ggsave(filename = filename_out,plot = p, dpi = 500)

}

# generate the figures for each patient
make_fish_diagram("3")
make_fish_diagram("21")
make_fish_diagram("6")
make_fish_diagram("9")


