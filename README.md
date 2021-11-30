# fmri_compare_dags

This class project compares the performance of PC, GES, ARGES, and notears for estimating DAGs in simulated brain network data from Smith et al. (2011). 

## Workflow Description
1. Load simulated data and save in workable format for R and Python using `load_and_save_network_data.ipynb`.

2. Run PC, GES, and ARGES based on the "CompareCausalNetworks" package. The R notebook `compare_causal_networks.Rmd` goes through this process and outputs graphs of the ground truth network, ground truth network without cycles, and estimated DAGs. It also outputs a performance metrics CSV for each simulation for each subject that contains the metrics `true positive rate`, `false positive rate`, `true discovery rate` both for obtaining the skeleton of the graph (undirected) and the directed graph. Because the ground truth network sometimes contained cycles, this R notebook only compared the estimated network with an acyclic version of ground truth network (i.e., removing the cycles in the ground trouth network). However, to get other relevant summary metrics such as the `structural hamming distance`, I used the Python toolbox `Causal Discovery Toolbox`; this required outputting the adjacency matrices into a format to read into Python. 

3. The Python notebook `additional_metrics_and_clean_up.ipynb` will now add the `structural hamming distance` metric to the performance metrics csv for each simulation. Since step 2 saved a lot of graphs and adjacency matrices for each subject for each simulation, this Python notebook also does some clean up by creating separate folders for each simulation and each subject, and then it moves the files appropriately. 

4. Since `notears` is available in Python, the implementation of this is in `run_notears.ipynb`. This notebook then adds the performance metrics into the performance metrics csv and saves the appropriate figures into each folder (from the cleanup process above).

5. Finally, the R notebook `summarize_dag_estimates.Rmd` generates summary statistics and creates summary figures showing the performance for each method across each simulation for each subject. These figures and summary statistics are reported in the project writeup.