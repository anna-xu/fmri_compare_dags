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

compute_metrics = function(ground_truth, estimate_graph, symm){
  # ========================================================================
  # takes in ground_truth, estimate_graph as adjacency matrices
  # outputs dataframe with tpr, fpr, tdr
  # ========================================================================
  if (symm){
    ground_truth = (ground_truth + t(ground_truth)) != 0
    estimate_graph = (estimate_graph + t(estimate_graph)) != 0 
  }
  
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
  metrics = data.frame(t(matrix(c(tpr, fpr, tdr))))
  names(metrics) = c("tpr", "fpr", "tdr")
  
  return(metrics)
}

compute_shd_with_matrices = function(ground_truth, estimate_graph, symm){
  # ========================================================================
  # `shd` function copied from the `Causal Discovery Toolbox` in Python
  # (adapted for use in R); sets double_for_anticausal as False
  # This means badly oriented edges are only counted as one mistake
  # This function is a simple implementation of the Structural Hamming Distance
  # With the option to only look at the skeleton if symm=True
  # ========================================================================
  if (symm){
    ground_truth = (ground_truth + t(ground_truth)) != 0
    estimate_graph = (estimate_graph + t(estimate_graph)) != 0 
  }
  
  diff = abs(ground_truth - estimate_graph)
  diff = diff + t(diff)
  diff[diff > 1] = 1
  return(sum(diff)/2)
}


run_single_subject_workflow = function(network_csv, num_nodes, ground_truth_csv, methods, sim_ind, subj_id, save){
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
  graph_metrics = as.data.frame(matrix(data=NA, nrow=1, ncol=6))
  names(graph_metrics) = c("tpr", "fpr", "tdr", "shd", "method", "is_directed")
  
  # for each method, 
  for (method_ind in 1:length(methods)){
    # get a graph estimate and convert it to matrix format
    graph_estimate = getParents(network_data, method=methods[method_ind], alpha=0.1)
    graph_estimate_mat = convert_to_matrix(graph_estimate)
    
    # plot and save the graphs
    graph_networks(graph_estimate_mat, paste0("simulation_", sim_ind, "_", subj_id, "_", methods[method_ind], "_estimate"))
    graph_networks(ground_truth_mat, paste0("simulation_",sim_ind, "_", subj_id, "_ground_truth"))
    graph_networks(ground_truth_acyclic_mat, paste0("simulation_", sim_ind, "_", subj_id, "_ground_truth_acyclic"))
    
    # compute summary metrics and aggregate
    # uses acyclic version of ground truth
    curr_graph_metrics = compute_metrics(ground_truth_acyclic_mat, graph_estimate_mat, symm=T)
    curr_graph_metrics$shd = compute_shd_with_matrices(ground_truth_acyclic_mat, graph_estimate_mat, symm=T)
    curr_graph_metrics$method = methods[method_ind]
    curr_graph_metrics$is_directed = "undirected"
    graph_metrics = rbind(graph_metrics, curr_graph_metrics)
    curr_graph_metrics = compute_metrics(ground_truth_acyclic_mat, graph_estimate_mat, symm=F)
    curr_graph_metrics$shd = compute_shd_with_matrices(ground_truth_acyclic_mat, graph_estimate_mat, symm=F)
    curr_graph_metrics$method = methods[method_ind]
    curr_graph_metrics$is_directed = "directed"
    graph_metrics = rbind(graph_metrics, curr_graph_metrics)
    
    # output graph estimates as a csv of the adjacency matrix
    # this will be useful for Python toolbox for other summary metrics
    write.csv(graph_estimate_mat, paste0("simulation_", sim_ind, "_", subj_id, "_", methods[method_ind], "_estimate.csv"))
  }
  
  graph_metrics = graph_metrics[-1,]
  
  if (save){
    write.csv(graph_metrics, paste0("simulation_", sim_ind, "_", subj_id, "_graph_metrics.csv")) 
  }
  
  return(graph_metrics)
}

run_single_simulation_workflow = function(sim_ind, save){
  # ========================================================================
  # Function that runs the single subject workflow for all subjects in sim
  # Takes in an integer as a string that indicates which simulation
  # Option to save output in CSV
  # Returns a dataframe with graph metrics
  # ========================================================================
  
  # get total number of subjects per simulation (assume will always be <100)
  file_paths_single_digit_subjs = Sys.glob(file.path(".", paste0("simulation_", sim_ind, "_[0-9].csv")))
  file_paths_double_digit_subjs = Sys.glob(file.path(".", paste0("simulation_", sim_ind, "_[0-9][0-9].csv")))
  num_subjects = length(file_paths_single_digit_subjs) + length(file_paths_double_digit_subjs)
  
  # initiate graph metrics dataframe
  all_graph_metrics = as.data.frame(matrix(data=NA, nrow=1, ncol=7))
  names(all_graph_metrics) = c("tpr", "fpr", "tdr", "shd", "method", "is_directed", "subj")
  # loop through each subject
  for (subj in 1:num_subjects){
    sim_data_csv_path = paste0("simulation_", sim_ind, "_", as.character(subj), ".csv")
    ground_truth_csv_path = paste0("ground_truth_", sim_ind, "_", as.character(subj), ".csv")
    num_nodes = ncol(read.csv(ground_truth_csv_path))-1 # have to do - 1 because 1st column is indexing
    curr_subj_graph_metrics = run_single_subject_workflow(sim_data_csv_path, num_nodes, ground_truth_csv_path, methods = c("pc","ges", "arges"), sim_ind, as.character(subj), save=FALSE)
    curr_subj_graph_metrics$subj = subj
    # add to all_graph_metrics
    all_graph_metrics = rbind(all_graph_metrics, curr_subj_graph_metrics)
  }
  
  all_graph_metrics = all_graph_metrics[-1,]
  
  if (save){
    write.csv(all_graph_metrics, paste0("simulation_", sim_ind, "_graph_metrics.csv"))
  }
  
  return(all_graph_metrics)
}

run_all_simulations = function(num_simulations, save){
  # ========================================================================
  # Function that loops through all subjects for all simulations 
  # Takes in an integer that specifies the number of simulations
  # Option to save csv of output
  # Returns a dataframe with graph metrics for all subjects, all simulations
  # ========================================================================
  all_graph_metrics = as.data.frame(matrix(data=NA, nrow=1, ncol=8))
  names(all_graph_metrics) = c("tpr", "fpr", "tdr", "shd", "method", "is_directed", "subj", "sim")
  
  for (curr_sim in 1:num_simulations){
    curr_sim_metrics = run_single_simulation_workflow(as.character(curr_sim), save=FALSE)
    curr_sim_metrics$sim = curr_sim
    all_graph_metrics = rbind(all_graph_metrics, curr_sim_metrics)
    
  }
  
  all_graph_metrics = all_graph_metrics[-1,]  
  
  if (save){
    write.csv(all_graph_metrics, "all_sims_subjs_graph_metrics.csv")
  }
  
  return(all_graph_metrics)
}