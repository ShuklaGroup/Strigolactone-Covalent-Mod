import numpy as np
import mdtraj as md
import pyemma
import itertools
import pickle
import glob
import os

def get_pairs(res_list):

    pairs = []
    p_all = itertools.combinations(res_list, 2)

    for pair in p_all:
        if np.abs(pair[1] - pair[0]) > 3:
            pairs.append(pair)

    pairs = np.array(pairs)

    return pairs

def compute_contact_distances(traj_list, pairs, top):

    dist_all = []
    for traj_file in traj_list:
        traj = md.load(traj_file, top=top)[::10] #Subsample because memory :(
        distances, contact_list = md.compute_contacts(traj, contacts=pairs, scheme='ca', ignore_nonprotein=False)
        dist_all.append(distances)

    return contact_list, dist_all

def calc_oasis(feat_all, n_components=50):

    if not os.path.exists("t_nystroem_obj.pkl"):

        tn_object = pyemma.coordinates.tica_nystroem(data=feat_all, lag=4, max_columns=n_components)
        pickle.dump(tn_object, open("t_nystroem_obj.pkl",'wb'))

    else:
        tn_object = pickle.load(open("t_nystroem_obj.pkl",'rb'))

    column_ind = tn_object.column_indices

    return column_ind


if __name__=="__main__":

    if os.path.exists("dist_all.pkl"):

        contact_list = np.load("contact_list.npy")
        dist_all = pickle.load(open("dist_all.pkl",'rb'))

    else:
        res_list = list(range(121,178)) + list(range(212,220)) + list(range(240,245)) #AtD14
        #res_list = list(range(123,181)) + list(range(215,223)) + list(range(243,248)) #ShHTL7
        pairs = get_pairs(res_list)

        traj_list = sorted(glob.glob("../stripped/*"))

        #contact_list, dist_all = compute_contact_distances(traj_list, pairs, "../stripped.4ih4_modelled_active.prmtop")
        #contact_list, dist_all = compute_contact_distances(traj_list, pairs, "../stripped.AtD14_D218A_apo.prmtop")
        #contact_list, dist_all = compute_contact_distances(traj_list, pairs, "../stripped.ShHTL7_apo_active.prmtop")
        contact_list, dist_all = compute_contact_distances(traj_list, pairs, "../stripped.AtD14_DOH-covH.prmtop")
        #contact_list, dist_all = compute_contact_distances(traj_list, pairs, "../stripped.ShHTL7_DOH-covH.prmtop")

        np.save("contact_list.npy", contact_list)
        pickle.dump(dist_all, open("dist_all.pkl",'wb'))

    #column_ind = calc_oasis(dist_all)
    #print(column_ind)
    #print(len(column_ind))

    #for ind in column_ind:
    #    print(contact_list[ind])
