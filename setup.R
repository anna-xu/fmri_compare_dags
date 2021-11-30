# ========================================================================
# Setup libraries
# ========================================================================
library(CompareCausalNetworks)
library(backShift)
library(dplyr)
library(pcalg)
library(graph)


# ========================================================================
# Helper functions
# ========================================================================

load_network = function(network_csv, num_nodes){
  # ========================================================================
  # Function that takes in a csv of network data and loads it
  # CSV can be nxp matrix where n is num_observations and p is num_nodes
  # Or CSV can be a single adjacency matrix
  # ========================================================================
  network = read.csv(network_csv)
  network = subset(network, select = names(network)[2:(num_nodes+1)])
  return(network)
}

convert_to_matrix = function(network_dataframe){
  # ========================================================================
  # converts dataframe of adjacency matrix to just matrix form
  # useful for graphs
  # ========================================================================
  network_matrix = as.matrix(unname(network_dataframe))
  return(network_matrix)
}

preproc_ground_truth = function(ground_truth_columns, ignore_cycles){
  # ========================================================================
  # takes in columns of the ground truth data frame and 
  # binarizes them (1 if there is an edge, 0 else, ignores sign)
  # if ignore_cycles==1, sets cycles to 0 by setting diagonals of matrix = 0
  # (useful if comparing with DAG estimates where there are no cycles)
  # ========================================================================
  bin_ground_truth = ifelse(ground_truth_columns != 0, 1, 0)
  bin_ground_truth_mat = convert_to_matrix(bin_ground_truth)
  if (ignore_cycles) {
    diag(bin_ground_truth_mat) = 0
  }
  return(bin_ground_truth_mat)
}

graph_networks = function(graph_matrix, title){
  # ========================================================================
  # function to plot graphs and save as output
  # ========================================================================
  network_graph = as(as.matrix(unname(graph_matrix)), "graphNEL")
  png(filename=paste0(title,".png"))
  plot(network_graph, main = title)
  dev.off()
}

compute_metrics = function(ground_truth, estimate_graph){
  # ========================================================================
  # takes in ground_truth, estimate_graph as adjacency matrices
  # outputs dataframe with tpr, fpr, tdr, and ppv 
  # (positive predictive value as measure of precision)
  # ========================================================================
  estimates_where_positive_ground_truth = estimate_graph[ground_truth == 1]
  true_positives = sum(estimates_where_positive_ground_truth)
  positives = length(estimates_where_positive_ground_truth)
  
  estimates_where_gap_ground_truth = estimate_graph[ground_truth == 0]
  false_positives = sum(estimates_where_gap_ground_truth)
  negatives = length(estimates_where_gap_ground_truth)
  
  predicted_positives = sum(estimate_graph==1)
  
  # metrics
  tpr = true_positives / positives
  fpr = false_positives / negatives
  tdr = true_positives / predicted_positives
  ppv = true_positives / (true_positives + false_positives)
  metrics = data.frame(t(matrix(c(tpr, fpr, tdr, ppv))))
  names(metrics) = c("tpr", "fpr", "tdr", "ppv")
  
  return(metrics)
}

run_single_subject_workflow = function(network_csv, num_nodes, ground_truth_csv, methods, sim_ind, save){
  # ========================================================================
  # Function takes in a network_csv for a single subject, the number of nodes,
  # Columns from an original csv containing the ground truth
  # A vector of the methods to be used (as strings)
  # A string that specifies which simulation this is
  # ========================================================================
  # load network data nxp matrix
  network_data = load_network(network_csv, num_nodes)
  # load ground truth
  ground_truth = read.csv(ground_truth_csv)
  ground_truth = subset(ground_truth, select = names(ground_truth)[2:(num_nodes+1)])
  ground_truth_columns = ground_truth[,1:num_nodes]
  ground_truth_mat = preproc_ground_truth(ground_truth_columns, ignore_cycles = FALSE)
  # load ground truth without cycles
  ground_truth_acyclic_mat = preproc_ground_truth(ground_truth_columns, ignore_cycles = TRUE)
  
  # initiate graph metrics dataframe to store summary metrics
  graph_metrics = as.data.frame(matrix(data=NA, nrow=1, ncol=5))
  names(graph_metrics) = c("tpr", "fpr", "tdr", "ppv", "method")
  
  # for each method, 
  for (method_ind in 1:length(methods)){
    # get a graph estimate and convert it to matrix format
    graph_estimate = getParents(network_data, method=methods[method_ind], alpha=0.1)
    graph_estimate_mat = convert_to_matrix(graph_estimate)
    
    # plot and save the graphs
    graph_networks(graph_estimate_mat, paste0("simulation_", sim_ind, "_", methods[method_ind], "_estimate"))
    graph_networks(ground_truth_mat, paste0("simulation_",sim_ind,"_ground_truth"))
    graph_networks(ground_truth_acyclic_mat, paste0("simulation_", sim_ind, "_ground_truth_acyclic"))
    
    # compute summary metrics and aggregate
    curr_graph_metrics = compute_metrics(ground_truth_mat, graph_estimate_mat)
    curr_graph_metrics$method = methods[method_ind]
    graph_metrics = rbind(graph_metrics, curr_graph_metrics)
    
    # output graph estimates as a csv of the adjacency matrix
    # this will be useful for Python toolbox for other summary metrics
    write.csv(graph_estimate_mat, paste0("simulation_", sim_ind, "_", methods[method_ind], "_estimate.csv"))
  }
  
  graph_metrics = graph_metrics[-1,]
  
  if (save){
    write.csv(graph_metrics, paste0("simulation_", sim_ind, "_graph_metrics.csv")) 
  }
}