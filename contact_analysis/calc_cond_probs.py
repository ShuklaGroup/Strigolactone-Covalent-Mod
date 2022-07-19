import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import glob
import pickle
import itertools

plt.rc('savefig', dpi=300)
matplotlib.rc('font',family='Helvetica-Normal',size=16)

def load_contact_distances(identifiers):
  
    file_lists = []

    for ft in feat_identifiers:
        file_lists.append(sorted(glob.glob("../analysis/%s"%ft)))

    feat_all = []
    for i in range(len(file_lists[0])): #Iterate through features
        feat = [] #Initialize feature for j-th trajectory
        for j in range(len(file_lists)):
            feat.append(np.load(file_lists[j][i])) #Load i-th feature for j-th trajectory

        feat_all.append(np.hstack(feat))

    feat_all = np.vstack(feat_all)

    return feat_all

def compute_probs(feat_all, cutoffs=[1.0, 1.0, 1.0, 1.0], conditions=[0, 0, 0, 0], weights=None):

    #Identify contacts
    contacts = np.zeros(np.shape(feat_all)) #Initialize array
    #print(contacts)

    for i in range(len(cutoffs)):
        if conditions[i] == 0:
            contacts[feat_all[:,i] < cutoffs[i],i] = 1
        else: 
            contacts[feat_all[:,i] > cutoffs[i],i] = 1

    #print(contacts)
    #Compute ligand contact probability by residue
    if weights is None:
        contact_probs = np.sum(contacts, axis=0)/np.shape(contacts)[0]
    else:
        contact_probs = np.matmul(contacts.T, weights)

    return contacts, contact_probs

def compute_conditional_probs(contacts, prob_index=2, condition_indices=[0,1], weights=None):

    switches = np.zeros((np.shape(contacts)[0],1))

    all_conditions = np.zeros((np.shape(contacts)[0],1))

    for i in range(np.shape(switches)[0]):
        conditions = []
        for j in range(len(condition_indices)):
            conditions.append(contacts[i,condition_indices[j]])

        if all(conditions):
            all_conditions[i] = 1 #Indicate if all conditions are true

        if all(conditions) and contacts[i, prob_index]:
            switches[i] = 1 

    if weights is None:
        switch_cond_probs = np.sum(switches, axis=0)/np.sum(all_conditions, axis=0)
    else:
        cond_probs = np.matmul(all_conditions.T, weights)
        switch_cond_probs = np.matmul(switches.T, weights)/cond_probs

    return switch_cond_probs

if __name__=="__main__":

    feat_identifiers = ["*helix-T1-T3-distance.npy","*helical-content-unfolded*","*T1-T2-hinge-contact-1*","*loop-distance-3*"]
    #feat_identifiers = ["*helix-T1-T3-distance-1.npy","*helical-content-unfolded*","*T1-T2-hinge-contact-1*","*loop-distance-3*"]
    feat_all = load_contact_distances(feat_identifiers)

    #weights_all = pickle.load(open("../msm_use/msm_weights.pkl", 'rb'))
    #weights_all = pickle.load(open("../oasis2/msm_weights.pkl", 'rb'))
    weights_all = pickle.load(open("../oasis4/clusters_350/msm_weights.pkl", 'rb'))

    weights = []
    for traj in weights_all:
        weights.extend(traj)

    weights = np.array(weights)

    #contacts, contact_probs = compute_probs(feat_all, cutoffs=[1.3, 0.6, 0.75, 1.0], conditions=[0, 0, 0, 1], weights=None)
    contacts, contact_probs = compute_probs(feat_all, cutoffs=[1.3, 0.6, 0.75, 1.0], conditions=[0, 0, 0, 1], weights=weights)

    print(contact_probs)

    del feat_all

    #Single conditions
    for i in range(len(feat_identifiers)):
        for j in range(len(feat_identifiers)):
            if j != i:
                #switch_cond_probs = compute_conditional_probs(contacts, prob_index=i, condition_indices=[j], weights=weights)
                switch_cond_probs = compute_conditional_probs(contacts, prob_index=i, condition_indices=[j], weights=None)
                print("Probability of %s given %s"%(feat_identifiers[i], feat_identifiers[j]))
                print(switch_cond_probs[0])

    #Double conditions
    for i in range(len(feat_identifiers)):
        for j in itertools.combinations(range(len(feat_identifiers)),2):
            if i not in j:
                #switch_cond_probs = compute_conditional_probs(contacts, prob_index=i, condition_indices=j, weights=weights)
                switch_cond_probs = compute_conditional_probs(contacts, prob_index=i, condition_indices=j, weights=None)
                print("Probability of %s given %s, %s"%(feat_identifiers[i], feat_identifiers[j[0]], feat_identifiers[j[1]]))
                print(switch_cond_probs[0])

    #Triple conditions
    for i in range(len(feat_identifiers)):
        for j in itertools.combinations(range(len(feat_identifiers)),3):
            if i not in j:
                #switch_cond_probs = compute_conditional_probs(contacts, prob_index=i, condition_indices=j, weights=weights)
                switch_cond_probs = compute_conditional_probs(contacts, prob_index=i, condition_indices=j, weights=None)
                print("Probability of %s given %s, %s, %s"%(feat_identifiers[i], feat_identifiers[j[0]], feat_identifiers[j[1]], feat_identifiers[j[2]]))
                print(switch_cond_probs[0])



