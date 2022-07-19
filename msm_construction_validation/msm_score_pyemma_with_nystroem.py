import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pyemma
import glob
import pickle

def load_data(feat_all, column_ind):

    reduced_set = []

    for traj in feat_all:
        reduced_set.append(traj[:,column_ind])

    return feat_all

def single_run(feat_all, n_tica_components=10, n_clusters=100):

    tica_object = pyemma.coordinates.tica(data=feat_all, lag=4, dim=n_tica_components)
    tica_trajs = tica_object.get_output()

    cluster_object = pyemma.coordinates.cluster_mini_batch_kmeans(tica_trajs, max_iter=300, k=n_clusters)
    dtrajs = cluster_object.assign(tica_trajs)

    msm_object = pyemma.msm.estimate_markov_model(dtrajs, 30, dt_traj='1000 ps', score_method='VAMP1', score_k=10)

    score = msm_object.score_cv(dtrajs, score_method='VAMP1', score_k=5)

    return score

def parameter_search(feat_all, tica_list, clusters_list):

    score_results = []
    for n_tica in tica_list:
        for n_cls in clusters_list:
            try:
                score = single_run(feat_all, n_tica_components=n_tica, n_clusters=n_cls)
                score_results.append((n_tica, n_cls, score))
                pickle.dump(score_results, open("score_results.pkl",'wb'))
            except:
                print("WARNING: Could not compute score for paramter set %d TICs and %d clusters"%(n_tica, n_cls))

    return score_results

if __name__=="__main__":

    #tn_obj = pickle.load(open("t_nystroem_obj_80.pkl",'rb'))
    tn_obj = pickle.load(open("t_nystroem_obj_100.pkl",'rb'))
    feat_all = pickle.load(open("dist_all.pkl",'rb'))

    feat_nystroem = load_data(feat_all, tn_obj.column_indices)
    #score_results = parameter_search(feat_all, [1,2,3,4,5,6,7,8,9,10], [100, 150, 200, 250, 300, 350, 400, 450])
    score_results = parameter_search(feat_nystroem, range(2,11,2), range(100,601,25))
    #score_results = parameter_search(feat_all, [3,4], [100,200])
    pickle.dump(score_results, open("score_results_final.pkl",'wb'))

