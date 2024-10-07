#!/usr/bin/env python3

import os
import re
import numpy as np
import bottleneck as bn
import pandas as pd
from scipy.special import gamma, binom
from scipy.stats import chi2
from scipy.cluster.hierarchy import linkage, cut_tree, leaves_list
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics.cluster import v_measure_score
from sklearn.cluster import AgglomerativeClustering

EPSILON = np.finfo(np.float64).resolution
log_EPSILON = np.log(EPSILON)

DOT_HEADER = 'digraph G {\n' \
    'node [width=0.75 fillcolor="#a6cee3", style=filled, fontcolor=black, ' \
    'shape=circle, fontsize=20, fontname="arial", fixedsize=True];\n' \

DOT_CELLS = 'node [width=0.5, fillcolor="#e8bdc9", fontcolor=black, ' \
    'style=filled, shape=square, fontsize=8, fontname="arial", fixedsize=True];\n'


# ------------------------------------------------------------------------------
# Mathematical functions
# ------------------------------------------------------------------------------

def check_beta_params(mean, var):
    ''' Check if parameters can be used for a beta function

    Arguments:
        mean (float): Beta function mean
        var (float): Beta function variance

    Returns:
        bool: True if parameters can be used, False otherwise

    '''
    return mean > .5 * (1 - (1 - 4 * var) ** .5)


# ------------------------------------------------------------------------------
# Evaluation functions
# ------------------------------------------------------------------------------

def get_v_measure(pred_clusters, true_clusters, out_file=''):
    score = v_measure_score(true_clusters, pred_clusters)
    if out_file:
        _write_to_file(out_file, score)
    return score


def get_ARI(pred_clusters, true_clusters, out_file=''):
    score = adjusted_rand_score(true_clusters,pred_clusters)
    if out_file:
        _write_to_file(out_file, score)
    return score


def get_hamming_dist(df_pred, df_true):
    if df_true.shape != df_pred.shape:
        score = np.count_nonzero(df_pred.round() != df_true.T)
    else:
        score = np.count_nonzero(df_pred.round() != df_true)
        # To catch caes where NxN dataframes got mixed up
        score_t = np.count_nonzero(df_pred.round() != df_true.T)
        if score_t < score:
            score = score_t
    return score


def _get_genotype_all(df_in, assign):
    df_out = pd.DataFrame(
        index=np.arange(df_in.shape[0]), columns=np.arange(len(assign))
    )
    if df_in.shape == df_out.shape:
        return df_in

    for cell_id, cl_id in enumerate(assign):
        if isinstance(cl_id, tuple):
            df_out[cell_id] = df_in[list(cl_id)].max(axis=1)
        else:
            df_out[cell_id] = df_in[cl_id]
    return df_out


def get_dist(assignments):
    steps, cells = assignments.shape
    dist = np.zeros(np.arange(cells).sum(), dtype=np.int32)
    # Sum up Hamming distance between cells for each posterior sample
    for assign in assignments:
        dist += pdist(np.stack([assign, assign]).T, 'hamming').astype(np.int32)
    # Return mean posterior cellwise hamming distance
    return dist / steps


def _get_MPEAR(assignments):
    dist = get_dist(assignments)
    sim = 1 - dist
    # dist = squareform(dist)
    Z = linkage(dist, method='ward')
    
    cl_no = []
    for assignment in assignments:
        cl_no.append(
            len([i for i in zip(*np.unique(assignment, return_counts=True)) \
                if i[1] > 2]))
    avg_cl_no = np.mean(cl_no)
    
    n_range = np.arange(max(2, avg_cl_no * 0.2),
        min(avg_cl_no * 2.5, assignments.shape[1]), dtype=int)

    best_MPEAR = -np.inf
    best_assignment = None

    for n in n_range:
        # model = AgglomerativeClustering(
        #     metric='precomputed', n_clusters=n, linkage='complete'
        # ).fit(dist)
        clusters = cut_tree(Z, n_clusters=n).flatten()
        score = _calc_MPEAR(sim, clusters)
        if score > best_MPEAR:
            best_assignment = clusters
            best_MPEAR = score

    return best_assignment


def _calc_MPEAR(pi, c):
    # Fritsch, A., Ickstadt, K. (2009) - Eq. 13
    I = 1 - pdist(np.stack([c, c]).T, 'hamming')

    I_sum = I.sum()
    pi_sum = pi.sum()
    index = (I * pi).sum()

    expected_index = (I_sum * pi_sum) / binom(c.size, 2)
    max_index = .5 * (I_sum + pi_sum)

    return (index - expected_index) / (max_index - expected_index)


def get_mean_hierarchy_assignment(assignments, params_full):
    steps = assignments.shape[0]
    assign = _get_MPEAR(assignments)
    clusters = np.unique(assign)

    params = np.zeros((clusters.size, params_full.shape[2]))
    for i, cluster in enumerate(clusters):
        cells_cl_idx = assign == cluster
        cells = np.nonzero(cells_cl_idx)[0]
        other = np.nonzero(~cells_cl_idx)[0]
        # Paper - section 2.3: first criteria
        if cells.size == 1:
            same_cluster = np.ones(steps).astype(bool)
        else:
            same_cluster = 0 == bn.nansum(
                bn.move_std(assignments[:, cells], 2, axis=1), axis=1
            )
        # Paper - section 2.3: second criteria
        cl_ids = np.array(
            [np.argmax(np.bincount(i)) for i in assignments[:,cells]])
        other_cl_id = assignments[:,other]
        no_others = [cl_ids[j] not in other_cl_id[j] for j in range(steps)]

        # At least criteria 1 fullfilled
        if any(same_cluster):
            # Both criteria fullfilled in at least 1 posterior sample
            if any(same_cluster & no_others):
                step_idx = np.argwhere(same_cluster & no_others).flatten()
            else:
                step_idx = np.argwhere(same_cluster).flatten()

            for step in step_idx:
                all_cl_ids = np.append(np.unique(other_cl_id[step]), cl_ids[step])
                rel_cl_id = np.argwhere(np.sort(all_cl_ids) == cl_ids[step])[0][0]
                params[i] += params_full[step][rel_cl_id]
            params[i] /= step_idx.size
        # If not, take parameters from all posterior samples
        else:
            for step, step_assign in enumerate(assignments):
                cl_id_all = np.unique(step_assign)
                cl_id, cnt = np.unique(step_assign[cells], return_counts=True)
                cl_id_new = np.argwhere(np.in1d(cl_id_all, cl_id)).flatten()
                params[i] += np.dot(cnt, params_full[step][cl_id_new])
            params[i] /= steps * cells.size

    params_df = pd.DataFrame(params).T[assign]
    return assign, params_df


def get_latents_posterior(results, data, single_chains=False):
    latents = []
    if single_chains:
        for result in results:
            latents.append(_get_latents_posterior_chain(result, data))
    else:
        result = _concat_chain_results(results)
        latents.append(_get_latents_posterior_chain(result, data))
    return latents


def _concat_chain_results(results):
    assign = np.concatenate([i['assignments'][i['burn_in']:] for i in results])
    a = np.concatenate([i['DP_alpha'][i['burn_in']:] for i in results])
    ML = np.concatenate([i['ML'][i['burn_in']:] for i in results])
    MAP = np.concatenate([i['MAP'][i['burn_in']:] for i in results])
    FN = np.concatenate([i['FN'][i['burn_in']:] for i in results])
    FP = np.concatenate([i['FP'][i['burn_in']:] for i in results])

    # Fill clusters not used by all chains with zeros
    params = [i['params'] for i in results]
    cl_max = np.max([i.shape[1] for i in params])
    for i, par_chain in enumerate(params):
        cl_diff = cl_max - par_chain.shape[1]
        params[i] = np.pad(par_chain, [(0, 0), (0, cl_diff), (0, 0)])
    par = np.concatenate(params)

    return {'assignments': assign, 'params': par, 'DP_alpha': a, 'FN': FN,
        'FP': FP, 'burn_in': 0, 'ML': ML, 'MAP': MAP}


def _get_latents_posterior_chain(result, data):
    burn_in = result['burn_in']
    assign, geno = get_mean_hierarchy_assignment(
        result['assignments'][burn_in:], result['params'][burn_in:]
    )
    a = _get_posterior_avg(result['DP_alpha'][burn_in:])
    FN = _get_posterior_avg(result['FN'][burn_in:])
    FP = _get_posterior_avg(result['FP'][burn_in:])

    FN_geno = (((geno.T.values.round() == 1) & (data == 0)).sum() + EPSILON)\
        / (geno.values.round().sum() + EPSILON)
    FP_geno = (((geno.T.values.round() == 0) & (data == 1)).sum() + EPSILON)\
        / ((1 - geno.values.round()).sum() + EPSILON)

    return {'a': a, 'assignment': assign, 'genotypes': geno, 'FN': FN, 'FP': FP,
        'FN_geno': FN_geno, 'FP_geno': FP_geno}


def _get_posterior_avg(data):
    return np.mean(data), np.std(data)


def get_latents_point(results, est, data, single_chains=False):
    latents = []
    if single_chains:
        for result in results:
            latents.append(_get_latents_point_chain(result, est, data))
    else:
        scores = [np.max(i[est][i['burn_in']:]) for i in results]
        best_chain = results[np.argmax(scores)]
        latents.append(_get_latents_point_chain(best_chain, est, data))

    return latents


def _get_latents_point_chain(result, est, data):
    step_no_bi = np.argmax(result[est][result['burn_in']:])
    step = step_no_bi + result['burn_in']

    # DPMM conc. parameter
    a = result['DP_alpha'][step]
    FP = result['FP'][step]
    FN = result['FN'][step]
    assignment = result['assignments'][step].tolist()

    cl_names = np.unique(assignment)

    geno_all = result['params'][step_no_bi][np.arange(cl_names.size)] # step_no_bi + 1
    geno = pd.DataFrame(geno_all, index=cl_names).T[assignment]

    FN_geno = (((geno.T.values.round() == 1) & (data == 0)).sum() + EPSILON) \
        / (geno.values.round().sum() + EPSILON)
    FP_geno = (((geno.T.values.round() == 0) & (data == 1)).sum() + EPSILON) \
        / ((1 - geno.values.round()).sum() + EPSILON)

    return {'step': step, 'a': a, 'assignment': assignment, 'genotypes': geno,
        'FN': FN, 'FP': FP, 'FN_geno': FN_geno, 'FP_geno': FP_geno}


def _write_to_file(file, content, attach=False):
    if attach and os.path.exists(file):
        open_flag = 'a'
    else:
        open_flag = 'w'

    with open(file, open_flag) as f:
        f.write(str(content))


def newick_to_gv(in_file, out_file=''):
    with open(in_file, 'r') as f:
        tree = f.read().strip().rstrip(';')

    edges, cells = get_edges_from_newick(tree)
    gv_tree = edges_to_gv(edges, cells)

    if out_file:
        _write_to_file(out_file, gv_tree)
    else:
        return gv_tree


def get_edges_from_newick(data):
    cells = sorted(re.findall(r'\w+cell\d*', data))
    for i, cell in enumerate(cells):
        data = data.replace(cell, f'C{i}')

    edges = []
    node_no = len(cells)

    while True:
        pairs = re.findall(r'\((C\d+):(0.\d+),(C\d+):(0.\d+)\)', data)
        if not pairs:
            break
        for i, pair in enumerate(pairs):
            n1, d1, n2, d2 = pair
            edges.append((node_no, int(n1.lstrip('C')), float(d1)))
            edges.append((node_no, int(n2.lstrip('C')), float(d2)))

            data = data.replace('({}:{},{}:{})'.format(*pair), f'C{node_no}')
            node_no += 1

    return edges, cells


def get_edges_from_gz(data):
        mut_edges = []
        muts = set([])
        cell_edges = []
        cells = []

        for line in data.split(';\n')[1:-1]:
            edge_nodes = re.search(r'(\d+)\s+->\s+(\d+)', line)
            attachment_nodes = re.search(r'(\d+)\s+->\s+(s\d+)', line)
            single_node = re.search(r'(s?\d+)$', line)

            if edge_nodes:
                n_from = int(edge_nodes.group(1))
                n_to = int(edge_nodes.group(2))
                n_from -= 1
                n_to -= 1

                if n_from != -1 and n_to != -1:
                    mut_edges.append((n_from, n_to))
                muts.update([n_from, n_to])
            if attachment_nodes:
                n_from = int(attachment_nodes.group(1))
                n_to = attachment_nodes.group(2)
                n_from -= 1
                cell_edges.append((n_from, n_to))
                cells.append(n_to)

            elif single_node:
                node = single_node.group(1)
                if node.startswith('s'):
                    cells.append(cells)
                else:
                    muts.add(int(node) - 1)

        return mut_edges, muts, cell_edges, cells


def edges_to_gv(edges, cells):
    # GraphViy Header: Node style
    out_str = DOT_HEADER

    e_length = [i[2] for i in edges]
    e_scaled = np.ceil(e_length / np.max(e_length) * (100)).astype(int)

    for i, edge in enumerate(edges):
        try:
            n_to = cells[edge[1]]
        except IndexError:
            n_to = edge[1]
        out_str += '{} -> {} [label="{}"];\n' \
            .format(edge[0], n_to, ' ' * e_scaled[i])

    out_str += '}'
    return out_str


def collapse_cells_on_tree(data_folder, out_file=''):
    tree_file = os.path.join(data_folder, 'tree.gv')

    with open(tree_file, 'r') as f:
        tree_str = f.read()
    mut_edges, muts, cell_edges, cells = get_edges_from_gz(tree_str)

    cell_edges_collapse = {}
    for mut_from, cell_to in cell_edges:
        try:
            cell_edges_collapse[mut_from].append(cell_to)
        except KeyError:
            cell_edges_collapse[mut_from] = [cell_to]

    # GraphViy Header: Node style
    out_str = DOT_HEADER
    for mut_edge in mut_edges:
        out_str += '{} -> {};\n'.format(*mut_edge)

    out_str += DOT_CELLS
    i = 0
    for mut_from, cells_to in cell_edges_collapse.items():
        size = 0.5 + len(cells_to) * 1
        out_str += '{f} -> s{t} [label="{s}", size={s}];\n' \
            .format(f=mut_from, t=i, s=size)
        i += 1
    out_str += '}'

    if not out_file:
        out_file = os.path.join(data_folder, 'tree_collapsed.gv')

    _write_to_file(out_file, out_str)

    try:
        from graphviz import render
        render('dot', 'png', out_file)
    except ImportError:
        pass


def get_lugsail_batch_means_est(data_in, steps=None):
    m = len(data_in)
    T_iL = []
    s_i = []
    n_i = []

    for data_chain, burnin_chain in data_in:
        data = data_chain[burnin_chain:steps]
        if data.size < 9: # otherwise b // 3 == 0
            return np.inf
        # [chapter 2.2 in Vats and Knudson, 2018]
        n_ii = data.size
        b = int(n_ii ** (1/2)) # Batch size. Alternative: n ** (1/3)
        n_i.append(n_ii)

        chain_mean = bn.nanmean(data)
        T_iL.append(
            2 * get_tau_lugsail(b, data, chain_mean) \
            - get_tau_lugsail(b // 3, data, chain_mean)
        )
        s_i.append(bn.nanvar(data, ddof=1))

    T_L = np.mean(T_iL)
    s = np.mean(s_i)
    n = np.round(np.mean(n_i))

    sigma_L = ((n - 1) * s + T_L) / n

    # [eq. 5 in Vats and Knudson, 2018]
    try:
        R_L = np.sqrt(sigma_L / s)  
    except FloatingPointError:
        return np.inf
    else:
        return R_L


def get_tau_lugsail(b, data, chain_mean):
    a = data.size // b # Number of batches
    batch_mean = bn.nanmean(np.reshape(data[:a * b], (a, b)), axis=1)
    return (b / (a - 1)) * bn.nansum(np.square(batch_mean - chain_mean))


def get_cutoff_lugsail(e, a=0.05):
    M = (4 * np.pi * chi2.ppf(1 - a, 1)) / (gamma(1/2)**2 * e**2)
    return np.sqrt(1 + 1 / M)


if __name__ == '__main__':
    print('Here be dragons...')