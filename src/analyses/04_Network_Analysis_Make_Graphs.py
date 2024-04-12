import numpy as np
import pandas as pd
from itertools import permutations
import glob
import os
import pickle
import time
import json
import networkx as nx
from sknetwork.clustering import aggregate_graph
from sknetwork.linalg import normalize
from sknetwork.utils import get_membership
from sknetwork.visualization import svg_graph
from sknetwork.embedding import Spring

fp_deg= "../../results/Make_DEG_Lists/"
fp_net = "../../data/networks/"
fp_results = "../../results/Make_Clusters/"
mydirs = [f.path for f in os.scandir(fp_results) if f.is_dir()]

for adir in mydirs:
    fp_save = adir
    try:
        robust_adjmat = np.load(fp_save+"/robust_adjmat.npy",allow_pickle=True)[()]
    except FileNotFoundError:
        continue
    # load needed data
    robust_adjmat = np.load(fp_save+"/robust_adjmat.npy",allow_pickle=True)[()]
    robust_labels = np.load(fp_save+"/robust_labels.npy")
    df_clus = pd.read_csv(fp_save+"/cluster_membership.tsv",sep="\t")
    with open(fp_save+"/params.json","r") as handle:
        param_dict = json.load(handle)
    # make df_subset
    df_edge = pd.read_csv(fp_net + param_dict["net_type"], sep="\t")
    df_subset = df_edge[df_edge["Weight"]>param_dict["edge_thresh"]]
    # read in the up and down regulated genes
    mygroup1 = param_dict["groups"].split("_")[0]
    mygroup2 = param_dict["groups"].split("_")[1]
    df = pd.read_csv(f"{fp_deg}{mygroup1}_{mygroup2}_{param_dict['tissue']}_up-regulated.csv")
    up_genes =  df[df['Human gene stable ID'].notna()]["Human gene stable ID"].tolist()
    df = pd.read_csv(f"{fp_deg}{mygroup1}_{mygroup2}_{param_dict['tissue']}_down-regulated.csv")
    down_genes =  df[df['Human gene stable ID'].notna()]["Human gene stable ID"].tolist()
    # make the layouts for agg of robust network
    spring = Spring(2,position_init="spectral")
    embedding = spring.fit_transform(robust_adjmat)
    mycolors = ["blue","red","green","orange","purple","yellow","aqua","greenyellow","deeppink"]
    average = normalize(get_membership(robust_labels).T)
    position_aggregate = average.dot(embedding)
    labels_unique, counts = np.unique(robust_labels, return_counts=True)
    adjacency_aggregate = aggregate_graph(robust_adjmat,labels=robust_labels)
    image = svg_graph(adjacency_aggregate, position_aggregate, counts, labels=labels_unique,
                      display_node_weight=True, node_weights=counts,
                      height=500,width=500,filename=fp_save+"/agg_network_plot.svg")
    for acol in list(df_clus):
        clus_genes = df_clus[df_clus[acol].notna()][acol].tolist()
        df_tmp = df_subset[(df_subset["Gene1"].isin(clus_genes))&(df_subset["Gene2"].isin(clus_genes))]
        # Turn the dataframe into a networkx graph object and scipy sparse matrix
        edgelist_tmp = list(df_tmp.itertuples(index=False))
        edgelist_tmp = [f"{item[0]} {item[1]} {item[2]}" for item in edgelist_tmp]
        G_tmp = nx.parse_edgelist(edgelist_tmp, data=(('weight',float),))
        adjmat_tmp = nx.to_scipy_sparse_matrix(G_tmp)
        # get up_down label
        diffexp_tmp = []
        for anode in G_tmp.nodes():
            if anode in down_genes:
                diffexp_tmp.append(0)
            elif anode in up_genes:
                diffexp_tmp.append(1)
        spring_tmp = Spring(2,position_init="spectral")
        embedding_tmp = spring_tmp.fit_transform(adjmat_tmp)
        image = svg_graph(adjmat_tmp, embedding_tmp, labels=diffexp_tmp,
                          label_colors=mycolors[0:2],
                          height=500,width=500,display_edge_weight=False,
                          filename=fp_save+f"/{acol}_network_plot.svg")