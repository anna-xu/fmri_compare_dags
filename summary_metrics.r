#==============================================================
# SETUP
#==============================================================

library(dplyr)
library(ggplot2)
library(gridExtra)

graph_metrics = read.csv("final_sims_subjs_graph_metrics.csv")

#==============================================================
# HELPER FUNCTIONS FOR PLOTTING
#==============================================================

plot_metric_all_sim = function(directed, metric, num_sims){
  sim_plots = list()
  for (curr_sim in 1:num_sims){
    sim_plots[[curr_sim]] = plot_metric_in_sim(directed, curr_sim, metric)
  }
  do.call(grid.arrange, sim_plots)
}

plot_metric_in_sim = function(directed, sim_id, metric){
  # =============================================================
  # Function to plot distribution of metric value across methods
  # (for a single simulation)
  # =============================================================
  
  df_subset = graph_metrics %>% subset(is_directed==directed) %>% subset(sim==sim_id)
  
  if (metric=="shd"){
   plot <- df_subset %>%
    ggplot(aes(x=method, y=shd, fill=method)) 
  } else if (metric=="tpr"){
    plot <- df_subset %>%
    ggplot(aes(x=method, y=tpr, fill=method))
  } else if (metric == "fpr"){
    plot <- df_subset %>%
    ggplot(aes(x=method, y=fpr, fill=method))
  } else {
    plot <- df_subset %>%
    ggplot(aes(x=method, y=tdr, fill=method))
  }
  
  plot <- plot + 
    ggtitle(as.character(sim_id)) +
    geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.1) + 
    geom_jitter(shape=16, position=position_jitter(0.15)) + 
    stat_summary(fun.y=median, shape=15, size=1, color="white") +
    theme(plot.title=element_text(size=100, face="bold", margin=margin(t=40,b=-30)),
          legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=60, face="bold"))
  
  return(plot)
}

run_and_save_plot = function(directed, method, num_sims){
  # =============================================================
  # Function to run and save the plots for metric across sims
  # =============================================================
  png(paste0(method, "_", directed, ".png"),
      width = 7000,
      height = 7000,
      units = "px")
  plot_metric_all_sim(directed, method, 28)
  dev.off()
}


# =====================================================================
# EXECUTE FIGURE GENERATION
# =====================================================================
run_and_save_plot("directed", "shd", 28)
run_and_save_plot("directed", "tpr", 28)
run_and_save_plot("directed", "tdr", 28)
run_and_save_plot("directed", "fpr", 28)

# Options for undirected graphs (did not run for project)
#run_and_save_plot("undirected", "tpr", 28)
#run_and_save_plot("undirected", "tdr", 28)
#run_and_save_plot("undirected", "fpr", 28)
#run_and_save_plot("undirected", "shd", 28)

# =====================================================================
# CSV FILES OF MEDIAN FOR BASELINE SIMULATIONS FOR EACH METHOD
# =====================================================================

#==============================================================
# TPR median and SE for simulations 1-4 and 17
#==============================================================

summary_graph = as.data.frame(matrix(data=NA, nrow=1, ncol=4))

for (curr_sim in 1:4){
  curr_graph = graph_metrics %>%
    subset(is_directed=="directed") %>%
    subset(sim==curr_sim) %>%
    group_by(method) %>%
    summarise(median(tpr), se=sd(tpr)/length(tpr)) 
  curr_graph$sim = curr_sim
  names(summary_graph) = names(curr_graph)
  summary_graph = rbind(summary_graph, curr_graph)
}

# for simulation 17 only
curr_graph = graph_metrics %>%
  subset(is_directed=="directed") %>%
  subset(sim==17) %>%
  group_by(method) %>%
  summarise(median(tpr), se=sd(tpr)/length(tpr)) 
curr_graph$sim = 17
names(summary_graph) = names(curr_graph)
summary_graph = rbind(summary_graph, curr_graph)

summary_graph = summary_graph[-1,]
write.csv(summary_graph, "summary_metrics_tpr.csv")

#==============================================================
# TPR median and SE for simulations 1-4 and 17
#==============================================================

summary_graph = as.data.frame(matrix(data=NA, nrow=1, ncol=4))

for (curr_sim in 1:4){
  curr_graph = graph_metrics %>%
    subset(is_directed=="directed") %>%
    subset(sim==1) %>%
    group_by(method) %>%
    summarise(median(tpr), se=sd(tpr)/length(tpr)) 
  curr_graph$sim = curr_sim
  names(summary_graph) = names(curr_graph)
  summary_graph = rbind(summary_graph, curr_graph)
}

# for simulation 17 only
curr_graph = graph_metrics %>%
  subset(is_directed=="directed") %>%
  subset(sim==17) %>%
  group_by(method) %>%
  summarise(median(tpr), se=sd(tpr)/length(tpr)) 
curr_graph$sim = 17
names(summary_graph) = names(curr_graph)
summary_graph = rbind(summary_graph, curr_graph)

summary_graph = summary_graph[-1,]
write.csv(summary_graph, "summary_metrics_tpr.csv")

#==============================================================
# TDR median and SE for simulations 1-4 and 17
#==============================================================

summary_graph = as.data.frame(matrix(data=NA, nrow=1, ncol=4))

for (curr_sim in 1:4){
  curr_graph = graph_metrics %>%
    subset(is_directed=="directed") %>%
    subset(sim==curr_sim) %>%
    group_by(method) %>%
    summarise(median(tdr), se=sd(tdr)/length(tdr)) 
  curr_graph$sim = curr_sim
  names(summary_graph) = names(curr_graph)
  summary_graph = rbind(summary_graph, curr_graph)
}

# for simulation 17 only
curr_graph = graph_metrics %>%
  subset(is_directed=="directed") %>%
  subset(sim==17) %>%
  group_by(method) %>%
  summarise(median(tdr), se=sd(tdr)/length(tdr)) 
curr_graph$sim = 17
names(summary_graph) = names(curr_graph)
summary_graph = rbind(summary_graph, curr_graph)

summary_graph_tdr = summary_graph[-1,]
write.csv(summary_graph_tdr, "summary_metrics_tdr.csv")

#==============================================================
# FPR median and SE for simulations 1-4 and 17
#==============================================================

summary_graph = as.data.frame(matrix(data=NA, nrow=1, ncol=4))

for (curr_sim in 1:4){
  curr_graph = graph_metrics %>%
    subset(is_directed=="directed") %>%
    subset(sim==curr_sim) %>%
    group_by(method) %>%
    summarise(median(fpr), se=sd(fpr)/length(fpr)) 
  curr_graph$sim = curr_sim
  names(summary_graph) = names(curr_graph)
  summary_graph = rbind(summary_graph, curr_graph)
}

# for simulation 17 only
curr_graph = graph_metrics %>%
  subset(is_directed=="directed") %>%
  subset(sim==17) %>%
  group_by(method) %>%
  summarise(median(fpr), se=sd(fpr)/length(fpr)) 
curr_graph$sim = 17
names(summary_graph) = names(curr_graph)
summary_graph = rbind(summary_graph, curr_graph)

summary_graph = summary_graph[-1,]
write.csv(summary_graph, "summary_metrics_fpr.csv")

#==============================================================
# SHD median and SE for simulations 1-4 and 17
#==============================================================

summary_graph = as.data.frame(matrix(data=NA, nrow=1, ncol=4))

for (curr_sim in 1:4){
  curr_graph = graph_metrics %>%
    subset(is_directed=="directed") %>%
    subset(sim==curr_sim) %>%
    group_by(method) %>%
    summarise(median(shd), se=sd(shd)/length(shd)) 
  curr_graph$sim = curr_sim
  names(summary_graph) = names(curr_graph)
  summary_graph = rbind(summary_graph, curr_graph)
}

# for simulation 17 only
curr_graph = graph_metrics %>%
  subset(is_directed=="directed") %>%
  subset(sim==17) %>%
  group_by(method) %>%
  summarise(median(shd), se=sd(shd)/length(shd)) 
curr_graph$sim = 17
names(summary_graph) = names(curr_graph)
summary_graph = rbind(summary_graph, curr_graph)

summary_graph = summary_graph[-1,]
write.csv(summary_graph, "summary_metrics_shd.csv")