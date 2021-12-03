# fmri_compare_dags

This class project for STATS 209 compares the performance of PC, GES, ARGES, and notears for estimating DAGs in simulated brain network data from Smith et al. (2011). As of 12/02/2021, 11:59 PM, this README reflects the latest version of the project (submitted). 

## Workflow Description
1. Load simulated data and save in workable format for R and Python using `load_and_save_network_data.ipynb`.

2. Run PC, GES, and ARGES based on the "CompareCausalNetworks" package. The R script `run_pc_ges_arges.r` goes through this process and outputs graphs of the ground truth network, ground truth network without cycles, and estimated DAGs. It also outputs a performance metrics CSV for each simulation for each subject that contains the metrics `true positive rate`, `false positive rate`, `true discovery rate` both for obtaining the skeleton of the graph (undirected) and the directed graph. Because the ground truth network sometimes contained cycles, this R notebook only compared the estimated network with an acyclic version of ground truth network (i.e., removing the cycles in the ground trouth network). However, to get other relevant summary metrics such as the `structural hamming distance`, I used the Python toolbox `Causal Discovery Toolbox` and adapted this to R. 

3. Since `notears` is available in Python, the implementation of this is in `run_notears.py`. This script outputs performance metrics and appends them to the performance metrics in step 2. It then saves plots of networks estimated from this algorithm.

4. Finally, the R script `summary_metrics.r` finalizes this project by outputting median scores, SE, and plots of the distribution of each performance metric across simulations. 

Running the scripts in sequential order will result in `final_sims_subjs_graph_metrics.csv`, `summary_metrics_{name of performance metric}.csv`,  and `{performance metric}_directed.png`. These files were used for the final project; this GitHub repo retains a copy of the `final_sims_subjs_graphs_metrics.csv` just in case (so that I don't need to run this again for whatever reason). There is an additional option built into these functions to also evaluate performance of undirected graphs (this was not reported in the final project but can be used for extensions in the future). 