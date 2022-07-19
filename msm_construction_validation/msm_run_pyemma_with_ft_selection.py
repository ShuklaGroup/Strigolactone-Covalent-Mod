import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import pyemma

import glob
import pickle
import time
import os

def calc_full_oasis_features(contact_list, column_ind, traj_list, top):

    pairs_use = contact_list[column_ind,:]

    dist_all = []

    for traj_file in traj_list:
        traj = md.load(traj_file, top=top)
        distances, contacts_use = md.compute_contacts(traj, contacts=pairs_use, scheme='ca', ignore_nonprotein=False)
        dist_all.append(distances)

    return contacts_use, dist_all

def calc_oasis(feat_all, n_components=50):

    if not os.path.exists("t_nystroem_obj_%d.pkl"%n_components):

        tn_object = pyemma.coordinates.tica_nystroem(data=feat_all, lag=4, max_columns=n_components)
        pickle.dump(tn_object, open("t_nystroem_obj_%d.pkl"%n_components,'wb'))

    else:
        tn_object = pickle.load(open("t_nystroem_obj_%d.pkl"%n_components,'rb'))

    column_ind = tn_object.column_indices

    return column_ind

def get_reduced_set(feat_all, column_ind):

    reduced_set = []
    
    for traj in feat_all:
        reduced_set.append(traj[:,column_ind])

    return reduced_set

def run_pipeline(feat_all, lag_time, tica_components, n_clusters):

    params_used = open("Parameters.txt",'w')
    params_used.write("Lag time: %d \n"%lag_time)
    params_used.write("TICA components: %d \n"%tica_components)
    params_used.write("Clusters: %d"%n_clusters)
    params_used.close()

    localtime = time.asctime(time.localtime(time.time()))
    print("Began TICA at: ",localtime)

    tica_object = pyemma.coordinates.tica(data=feat_all, lag=4, dim=tica_components)
    tica_trajs = tica_object.get_output()
    pickle.dump(tica_object, open("tica_object.pkl",'wb'))
    pickle.dump(tica_trajs, open("tica_trajs.pkl",'wb'))

    localtime = time.asctime(time.localtime(time.time()))
    print("Began clustering at: ",localtime)

    cluster_object = pyemma.coordinates.cluster_mini_batch_kmeans(tica_trajs, max_iter=200, k=n_clusters)
    dtrajs = cluster_object.assign(tica_trajs)
    pickle.dump(dtrajs, open("dtrajs.pkl",'wb'))

    #cluster_object = pyemma.coordinates.cluster_mini_batch_kmeans(feat_all, k=n_clusters, max_iter=500)
    #dtrajs = cluster_object.assign(feat_all)
    #pickle.dump(dtrajs, open("dtrajs.pkl",'wb'))

    its_object = pyemma.msm.its(dtrajs, errors='bayes', lags=80, nits=10)
    pyemma.plots.plot_implied_timescales(its_object, outfile="its_plot.png", units='ns', dt=0.1)
    #pyemma.plots.plot_implied_timescales(its_object, outfile="its_plot.png", units='ns', dt=1.0)

    localtime = time.asctime(time.localtime(time.time()))
    print("Began estimating MSM at: ",localtime)

    msm_object = pyemma.msm.estimate_markov_model(dtrajs, lag_time, dt_traj='100 ps', score_method='VAMP1', score_k=10)
    #msm_object = pyemma.msm.estimate_markov_model(dtrajs, lag_time, dt_traj='1000 ps', score_method='VAMP1', score_k=10)
    pickle.dump(msm_object, open("msm_object.pkl",'wb'))
    #fig, pos = pyemma.plots.plot_markov_model(msm_object)
    #plt.savefig("network.png")

    localtime = time.asctime(time.localtime(time.time()))
    print("Data usage: ", msm_object.active_count_fraction)

    print("Began scoring at: ",localtime)
    score = msm_object.score_cv(dtrajs, score_method='VAMP1', score_k=5)
    print(score)

    localtime = time.asctime(time.localtime(time.time()))
    print("Finished scoring at: ",localtime)

    weights = msm_object.trajectory_weights()
    pickle.dump(weights, open("msm_weights.pkl",'wb'))

    #eigs = msm_object.eigenvalues()
    #plt.figure()
    #plt.plot(list(range(len(eigs))), eigs, 'o')
    #plt.savefig("eigenvalues.png")

    #cktest = msm_object.cktest(5, mlags=range(4))
    #fig, ax = pyemma.plots.plot_cktest(cktest)
    #fig.savefig("cktest.png")

    return score

def select_ftr_size(ft_all, lag_time, tica_components, n_clusters):

    scores_all = []

    for i in range(10,151,10):

        column_ind_all = calc_oasis(ft_all, n_components=i)
        print(len(column_ind_all))
        reduced_ft = get_reduced_set(ft_all, column_ind_all)
        score = run_pipeline(reduced_ft, lag_time, tica_components, n_clusters)
        scores_all.append(score)

    pickle.dump(scores_all, open("scores_all.pkl",'wb'))

    return scores_all

if __name__=="__main__":

    localtime = time.asctime(time.localtime(time.time()))
    print("Began run at: ",localtime)

    ##MSM Parameters
    lag_time = 300 #Trajectory steps
    #lag_time = 30 #Trajectory steps
    #tica_components = 8
    #tica_components = 6
    tica_components = 4
    n_clusters = 375
    #n_clusters = 400
    #n_clusters = 300
    #n_clusters = 325

    ft_all = pickle.load(open("dist_all.pkl",'rb'))
    #ft_all = pickle.load(open("../dist_all.pkl",'rb'))
    #column_ind = calc_oasis(ft_all, n_components=120)
    column_ind = calc_oasis(ft_all, n_components=100)

    #scores_all = select_ftr_size(ft_all, lag_time, tica_components, n_clusters)
    
    try:
        ft_oasis_selected = pickle.load(open("oasis_ftr_full.pkl",'rb'))
    except:
        contact_list = np.load("contact_list.npy")
        traj_list = sorted(glob.glob("../stripped/*"))
        #traj_list = sorted(glob.glob("tmd_trajs/*xtc"))
        top="../stripped.AtD14_DOH-covH.prmtop"
        #top="../stripped.ShHTL7_DOH-covH.prmtop"
        #print(traj_list)
        ft_oasis_selected = calc_full_oasis_features(contact_list, column_ind, traj_list, top)
        pickle.dump(ft_oasis_selected, open("oasis_ftr_full.pkl",'wb'))

    score = run_pipeline(ft_oasis_selected[1], lag_time, tica_components, n_clusters)
