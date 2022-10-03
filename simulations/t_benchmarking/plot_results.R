# Data processing ======

library(ggplot2)
library(dplyr)
library(scales)
directory <- c("fast_growth","slow_growth","t_benchmarking_no_death")

parameter_sets <- c("fast","slow","no death")
df_summary_list <- list()
df_list <- list()
for (i in 1:3){
  df <- read.csv(paste0(directory[i],"/t_benchmarking_simulation_results.txt"),header = T)
  df$PE <- 100*(df$t_estimate - df$t_true)/df$t_true
  df_list[[i]] <- df
  df_summary <- df %>% 
    group_by(tumor_size) %>% 
    summarise(mean = mean(PE),
              sd = sd(PE))
  df_summary$parameter_set <- parameter_sets[i]
  df_summary_list[[i]] <- df_summary
}
df_combined <- bind_rows(df_summary_list)

# combined plot =========
ggplot(data = df_combined) + 
  geom_hline(yintercept=0, color = "black", size=1)+
  geom_ribbon(aes(ymin = mean - sd, 
                  ymax = mean + sd,
                  x = tumor_size,
                  fill = factor(parameter_set),
                  color = factor(parameter_set)), 
              alpha = 0.3,
              linetype = "solid")+
  geom_line(aes(x = tumor_size, 
                y = mean,
                color = factor(parameter_set)),
            size = 1,
            alpha = .9)+
  geom_point(aes(x = tumor_size, 
                 y = mean,
                 fill = factor(parameter_set)),
             color = 'black',
             shape = 21,
             size = 3,
             alpha = 1.0)+
  scale_y_continuous(name = "percent error of t estimate")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                name = "number of cancer cells")+
  theme_bw()+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        aspect.ratio = .7)+
  scale_color_manual(values=c("tomato", "grey50","dodgerblue2"),
                     name = "parameter regime")+
  scale_fill_manual(values=c("tomato", "grey87", "dodgerblue2"),
                    name = "parameter regime")


# plotting separately ===============
p <- list()
x_labs <- c("","","number of cancer cells")
for (i in 1:3){
  x_lab <- x_labs[i]
  plot_name <- paste0(parameter_sets[i],"_t_benchmarking.svg")
  p[[i]] <- ggplot() + 
    geom_hline(yintercept=0, color = "grey", size=1)+
    geom_jitter(data = df_list[[i]], 
                aes(x = tumor_size, y = PE),
                width = 0.1,
                height = 0,
                size = .6,
                shape = 16,
                alpha = 0.3,
                fill = "grey")+
    geom_ribbon(data = df_summary_list[[i]], aes(ymin = mean - sd, 
                                       ymax = mean + sd,
                                       x = tumor_size), 
                #width=.1,
                #size = 1,
                fill = "firebrick3",
                alpha = 0.3)+
    geom_point(data = df_summary_list[[i]], 
               aes(x = tumor_size, 
                   y = mean),
               fill = "firebrick3",
               color = 'black',
               shape = 21,
               size = 3,
               alpha = .9)+
    geom_line(data = df_summary_list[[i]], 
              aes(x = tumor_size, 
                  y = mean),
              color = "firebrick3",
              size = 1,
              alpha = .9)+
    scale_y_continuous(name = "percent error")+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", scales::math_format(10^.x)),
                  name = x_lab)+
    theme_bw()+
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20),
          aspect.ratio = .7)
    ggsave(plot_name,p[[i]],bg = "transparent")
}
