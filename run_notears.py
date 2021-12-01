import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import networkx as nx
from notears.linear import notears_linear

#plt.figure()
#plt.plot(sims['ts'][2000:2200, :])
#plt.figure()
#plt.plot(ts[10])

def load_sims(i):
    path = f'sims/sim{i}.mat'
    sims = scipy.io.loadmat(path)
    ts = sims['ts'].reshape(
        int(sims['Nsubjects']),
        int(sims['Ntimepoints']),
        int(sims['Nnodes']),
    )
    nets = sims['net']
    return ts, nets

def tabulate_accuracy(ground_truth, estimate_graph, symm=False):
    if symm:
        ground_truth = (ground_truth + ground_truth.T) != 0
        estimate_graph = (estimate_graph + estimate_graph.T) != 0 
    
    estimates_where_positive_ground_truth = estimate_graph[ground_truth == 1]
    true_positives = np.sum(estimates_where_positive_ground_truth)
    positives = len(estimates_where_positive_ground_truth)

    estimates_where_gap_ground_truth = estimate_graph[ground_truth == 0]
    false_positives = np.sum(estimates_where_gap_ground_truth)
    negatives = len(estimates_where_gap_ground_truth)

    predicted_positives = np.sum(estimate_graph==1)

    diff = np.abs(ground_truth.astype(int) - estimate_graph.astype(int))
    diff = diff + diff.T
    diff[diff > 1] = 1
    shd = np.sum(diff)/2
    
    # metrics
    tpr = true_positives / positives
    fpr = false_positives / negatives
    tdr = true_positives / predicted_positives
    
    return {
        'tpr': tpr,
        'fpr': fpr,
        'tdr': tdr,
        'shd': shd
    }

def tune_params(ts, net):

    lambdas = [0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    losses = ['l2', 'poisson']

    true = (net - np.diag(np.diag(net)) != 0)
    acc = {}
    for lm in lambdas:
        for loss in losses:
            print(f'Running: lambda={lm}, loss func={loss}')
            out = notears_linear(
                X=ts,
                lambda1=lm,
                loss_type=loss,
                max_iter=1000
            )
            acc[(lm, loss)] = tabulate_accuracy(true, (np.abs(out) > 0.01))
    return acc

n_sims = 28

dfs = []
for sim in range(n_sims):
    print(f'\nSIMULATION {sim + 1}\n=======================')
    ts, nets = load_sims(sim + 1)
    n_sub = ts.shape[0]
    selected = np.random.permutation(n_sub)[:(n_sub // 25)]
    for s in selected:
        df_new = pd.DataFrame(tune_params(ts[s], nets[s]))
        df_new.index = pd.Index([(i, sim + 1, s) for i in df_new.index])
        dfs += [df_new]
df = pd.concat(dfs)
df.index = df.index.set_names(['Metric', 'Simulation', 'Subject'])
df = df.iloc[df.index.get_level_values('Metric') == 'shd']

df.mean(level='Simulation').idxmin(axis=1)

fig, ax = plt.subplots(figsize=(12, 9))
ax.plot(df.mean(level='Simulation').values.T)
ax.legend(range(1, 29), bbox_to_anchor=(0.2, 1.05),
           ncol=7, fancybox=True, shadow=True)
ax.set_xticks(range(12))
ax.set_xticklabels([f'lambda={la}, loss={lo}' for la, lo in df.columns], rotation=45, horizontalalignment='right')

def run_notears_linear(ts, net, ids):
    sim_id, sub_id = ids
    true = (net - np.diag(np.diag(net)) != 0)
    acc = {}
    out = notears_linear(
        X=ts,
        lambda1=0,
        loss_type='l2',
        max_iter=1000
    )
    #print(out)
    acc = tabulate_accuracy(true, (np.abs(out) > 0.01))
    acc_symm = tabulate_accuracy(true, (np.abs(out) > 0.01), symm=True)
    G = nx.DiGraph(out)
    plt.figure()
    nx.draw(G, with_labels=True, node_color='#DDDDDD', node_size=999)
    plt.savefig(f'simulation_{sim_id}_{sub_id}_notears_linear.png')
    return acc, acc_symm

n_sims = 28

accs = {
    'sim': [],
    'subj': [],
    'is_directed': [],
    'tpr': [],
    'fpr': [],
    'tdr': [],
    'shd': []
}
for sim in range(n_sims):
    print(f'\nSIMULATION {sim + 1}\n=======================')
    ts, nets = load_sims(sim + 1)
    n_sub = ts.shape[0]
    for s in range(n_sub):
        acc, acc_symm = run_notears_linear(ts[s], nets[s], (sim, s))
        accs['sim'] += [sim + 1, sim + 1]
        accs['subj'] += [s + 1, s + 1]
        accs['is_directed'] += ['directed', 'undirected']
        for k in acc.keys():
            accs[k] += [acc[k], acc_symm[k]]
df = pd.DataFrame(accs)

