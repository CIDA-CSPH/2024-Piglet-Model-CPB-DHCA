import numpy as np
import pandas as pd
from itertools import permutations
import glob
import os
import pickle
import time
import json
import networkx as nx
from sknetwork.clustering import Louvain, get_modularity
from sknetwork.topology import get_connected_components

def generate_clusters(mytiss,myres,mygenes,mygroups):    
    # set some inputs
    param_dict = {}
    info_file = []
    # from inputs
    param_dict["tissue"] = mytiss
    param_dict["sig_name"] = mygenes
    param_dict["res"] = myres
    param_dict["groups"] = mygroups
    mygroup1 = mygroups.split("_")[0]
    mygroup2 = mygroups.split("_")[1]
    # hard coded
    param_dict["n_runs"] = 10000
    param_dict["edge_thresh"] = 0.7
    param_dict["min_genes"] = 10
    param_dict["net_type"] = 'STRING.edg'
    param_dict["prob_thresh"] = 0.5
    param_dict["robust_thresh"] = 0.9
    param_dict["robust_min_genes"] = 10

    # set some filepaths
    fp_deg = "../../results/Make_DEG_Lists/"
    fp_net = "../../data/networks/"
    fp_results = "../../results/Make_Clusters/"
    save_time = int(time.time())
    fn_end = f"{save_time}__{param_dict['groups']}__{param_dict['tissue']}__{param_dict['sig_name']}__{param_dict['res']}"
    os.mkdir(fp_results+fn_end)
    save_FN = fp_results + f"{fn_end}/"

    # read in the up and down regulated genes
    df = pd.read_csv(f"{fp_deg}{mygroup1}_{mygroup2}_{param_dict['tissue']}_up-regulated.csv")
    up_genes =  df[df['Human gene stable ID'].notna()]["Human gene stable ID"].tolist()
    df = pd.read_csv(f"{fp_deg}{mygroup1}_{mygroup2}_{param_dict['tissue']}_down-regulated.csv")
    down_genes =  df[df['Human gene stable ID'].notna()]["Human gene stable ID"].tolist()
    print("The number of up_genes is ",len(up_genes))
    print("The number of down_genes is ",len(down_genes))

    # set the genes from params
    if param_dict["sig_name"] == "UpGenes":
        sig_genes = up_genes
    elif param_dict["sig_name"] == "DownGenes":
        sig_genes = down_genes
    elif param_dict["sig_name"] == "UpDownGenes":
        sig_genes = up_genes + down_genes

    # generate the clusters
    tic = time.time()
    # read in the STRING network
    df_edge = pd.read_csv(fp_net + param_dict["net_type"], sep="\t")
    info_file.append(f"The total number of edges in the network is {df_edge.shape[0]}")

    df_subset = df_edge[df_edge["Weight"]>param_dict["edge_thresh"]]
    info_file.append(f"The number of edges that pass the weight threshold is {df_subset.shape[0]}")
    df_subset = df_subset[(df_subset["Gene1"].isin(sig_genes))&(df_subset["Gene2"].isin(sig_genes))]
    info_file.append(f"The number of edges in genelist network is {df_subset.shape[0]}")
    info_file.append(f"The number of genes in genelist network is {len(np.unique(df_subset['Gene1'].tolist()+df_subset['Gene2'].tolist()))}")

    # Turn the dataframe into a networkx graph object and scipy sparse matrix
    edgelist = list(df_subset.itertuples(index=False))
    edgelist = [f"{item[0]} {item[1]} {item[2]}" for item in edgelist]
    G = nx.parse_edgelist(edgelist, data=(('weight',float),))
    adjmat = nx.to_scipy_sparse_matrix(G)

    # get n_run numbers for random generator
    run_seeds = np.random.choice(a=1000000, size=param_dict["n_runs"], replace=True)
    # Run Lovain n_runs amount of times
    mods = []
    robust_net = np.zeros(adjmat.shape)
    for run_idx, arun in enumerate(range(param_dict["n_runs"])):
        # Find a partition
        louvain = Louvain(shuffle_nodes=True,resolution=param_dict["res"],
                          modularity="newman",sort_clusters=True,
                          random_state=int(run_seeds[run_idx]))
        louvain.fit(adjmat)
        labels1 = louvain.predict()
        labels2 = louvain.predict_proba()
        # save the modularity score
        mods.append(get_modularity(adjmat, labels1))
        # add move genes with a maximum probability below threshold to null cluster (named 400)
        labels1_thresh = [item for item in labels1]
        for idx in range(len(labels1)):
            if np.max(labels2[idx,:]) < param_dict["prob_thresh"]:
                labels1_thresh[idx] = 400
        labels1_thresh = np.array(labels1_thresh)
        # find number of genes in each cluster and the indices of those genes in the adjmat object
        labels_unique1, counts1 = np.unique(labels1_thresh, return_counts=True)
        ind400 = np.where(labels_unique1==400)
        # delete those from the label and count objects
        labels_unique1 = np.delete(labels_unique1,ind400)
        counts1 = np.delete(counts1,ind400)
        # Subset labels and adjmat for only clusters above size of min_genes
        mininds = np.where(counts1>=param_dict["min_genes"])
        good_mat_inds = []
        good_labels = []
        good_label_dict = {}
        for idx, alabel in enumerate(labels1_thresh):
            if alabel in mininds[0]:
                good_labels.append(alabel)
                good_mat_inds.append(idx)
                if alabel not in good_label_dict:
                    good_label_dict[alabel] = [idx]
                else:
                    good_label_dict[alabel] = good_label_dict[alabel] + [idx]
        # subset the network to these larger clusters
        good_labels = np.array(good_labels)
        good_adjmat = adjmat[good_mat_inds,:]
        good_adjmat = good_adjmat[:,good_mat_inds]
        # save the modularity of this new subset of clustering
        mods.append(get_modularity(good_adjmat, good_labels))
        # add edges robust_network
        for akey in good_label_dict:
            edge_pairs = permutations(good_label_dict[akey],2)
            for apair in edge_pairs:
                robust_net[apair] = robust_net[apair] + 1
    # threshold robust_net and find connected components
    robust_net2 = robust_net / param_dict["n_runs"]
    robust_net3 = np.where(robust_net2 > param_dict["robust_thresh"],1,0)
    cc_labels = get_connected_components(robust_net3)
    cc_unique1, cc_counts = np.unique(cc_labels, return_counts=True)
    # subset robust_net to only include nodes in CCs bigger than robust_min_genes
    cc_mininds = np.where(cc_counts>=param_dict["robust_min_genes"])
    robust_inds = []
    robust_labels = []
    for idx, alabel in enumerate(cc_labels):
        if alabel in cc_mininds[0]:
            robust_labels.append(alabel)
            robust_inds.append(idx)
    robust_labels = np.array(robust_labels)
    robust_adjmat = adjmat[robust_inds,:]
    robust_adjmat = robust_adjmat[:,robust_inds]
    # reorder CC loabels so they start at zero and increase by one
    for idx, alabel in enumerate(np.unique(robust_labels)):
        robust_labels[robust_labels == alabel] = idx
    mods.append(get_modularity(robust_adjmat, robust_labels))
    # find up, down, added labels for the cluster nodes and which genes are in each cluster
    node_names = list(G.nodes())
    labels_diffexp = []
    genes_per_cluster = {}
    for idx, anind in enumerate(robust_inds):
        # get differential expression node labels
        node_tmp = node_names[anind]
        if node_tmp in down_genes:
            labels_diffexp.append(0)
        elif node_tmp in up_genes:
            labels_diffexp.append(1)
        else:
            labels_diffexp.append(2)
        # get genes per cluster
        if robust_labels[idx] not in genes_per_cluster:
            genes_per_cluster[robust_labels[idx]] = [node_tmp]
        else:
            genes_per_cluster[robust_labels[idx]] = genes_per_cluster[robust_labels[idx]] + [node_tmp]
    df_clus = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in genes_per_cluster.items() ]))
    info_file.append(f"The number of clusters is {df_clus.shape[1]}")
    info_file.append(f"The total number of genes is {np.sum([len(genes_per_cluster[akey]) for akey in genes_per_cluster])}")
    ## save some things
    # save cluster memberships from robust network
    df_clus.to_csv(save_FN+"cluster_membership.tsv",sep="\t",header=True,index=False)
    np.save(save_FN+"robust_adjmat.npy",robust_adjmat)
    np.save(save_FN+"robust_labels.npy",robust_labels)
    np.save(save_FN+"run_seeds.npy",run_seeds)
    np.save(save_FN+"modularities.npy",mods)
    with open(save_FN+"params.json", 'w') as handle:
        json.dump(param_dict, handle)
    info_file.append(f"The time it took in minutes is {int((time.time()-tic)/60)}")
    with open(save_FN+"info.json", 'w') as thefile:
        for item in info_file:
            thefile.write("%s\n" % item)
            
for atiss in ["Kidney","Liver","Lung","Heart","Ileum"]:
    for agenes in ["UpDownGenes","DownGenes","UpGenes"]: 
        for ares in [1]:
            for agroup in ["CPBnodrug_Controlnodrug","CPBLind_CPBnodrug","CPBLind_Controlnodrug"]:
                try:
                    generate_clusters(atiss,ares,agenes,agroup)
                except ValueError:
                    pass
                except FileNotFoundError: # this was added as not all groups had DEGs
                    pass
                except nx.NetworkXError: # this is beacuse sometimes no IDs map into the network
                    pass