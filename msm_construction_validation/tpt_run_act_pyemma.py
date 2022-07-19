import numpy as np
import mdtraj as md
import pyemma
import pickle
import glob

def get_samples(msm_obj, topfile, n=50):
    """Draw samples from each MSM state"""

    samples = msm_obj.sample_by_state(n)
    trajs = sorted(glob.glob("../stripped/*"))
    #Save xtc files
    for i in range(len(samples)):
        state_samples = samples[i]
        state_traj = []
        for j in range(len(state_samples)):
            state_traj.append(md.load_frame(trajs[state_samples[j,0]], state_samples[j,1], top=topfile))
        md.join(state_traj).save_xtc("state_%d.xtc"%i)

    return samples

def id_active_inactive_states(samples):
    """Identify ligand bound, unbound, and inverse bound states"""

    T1_T3_dist = sorted(glob.glob("../analysis/*T1-T3-distance*")) #Ligand pocket features
    Helical_Content = sorted(glob.glob("../analysis/*helical-content-unfolded*")) #Ligand pocket features
    T1_ext = sorted(glob.glob("../analysis/*T1-T2-hinge-contact-1*")) #Ligand pocket features
    D_loop = sorted(glob.glob("../analysis/*loop-distance-3*")) #Ligand pocket features
    
    active_states = []
    inactive_states = []
    
    for i in range(len(samples)):
        state_samples = samples[i]
        t1t3_dist = []
        helical = []
        t1_ext = []
        dloop_dist = []
        for j in range(len(state_samples)):
            t1t3_dist.append(np.load(T1_T3_dist[state_samples[j,0]])[state_samples[j,1]])
            helical.append(np.load(Helical_Content[state_samples[j,0]])[state_samples[j,1]])
            t1_ext.append(np.load(T1_ext[state_samples[j,0]])[state_samples[j,1]])
            dloop_dist.append(np.load(D_loop[state_samples[j,0]])[state_samples[j,1]])
        if all((np.mean(t1t3_dist) < 1.3, np.mean(helical) < 0.6, np.mean(t1_ext) < 0.75, np.mean(dloop_dist) > 1.0)):
            active_states.append(i)
        elif all((np.mean(t1t3_dist) > 2.0, np.mean(helical) > 0.7, np.mean(t1_ext) > 0.9, np.mean(dloop_dist) < 0.8)):
            inactive_states.append(i)

    return active_states, inactive_states

def coarse_grain(msm):

    pcca_object = msm.pcca(8)

    return pcca_object

def calc_tpt(msm, source, sink):

    tpt_object = pyemma.msm.tpt(msm, source, sink)
    mfpt = tpt_object.mfpt
    intermediates = tpt_object.I
    pathways = tpt_object.pathways

    return tpt_object, mfpt, intermediates, pathways

def pathway_fluxes(pathways, capacities, anchored):

    sticky_paths = []
    direct_paths = []
    sticky_flux = []
    direct_flux = []

    for i in range(len(pathways)):
        if any(j in anchored for j in pathways[i]):
            sticky_paths.append(pathways[i])
            sticky_flux.append(capacities[i])
        else:
            direct_paths.append(pathways[i])
            direct_flux.append(capacities[i])
    
    return sticky_paths, direct_paths, sticky_flux, direct_flux

if __name__=="__main__":

    msm = pickle.load(open("../oasis4/clusters_350/msm_object.pkl",'rb'))
    #msm = pickle.load(open("../oasis2/msm_object.pkl",'rb'))

    try:
        samples = pickle.load(open("state_samples.pkl",'rb'))
    except:
        #samples = get_samples(msm, topfile="../stripped.AtD14_DOH-covH.prmtop")
        samples = get_samples(msm, topfile="../stripped.ShHTL7_DOH-covH.prmtop")
        pickle.dump(samples, open("state_samples.pkl",'wb'))

    active, inactive = id_active_inactive_states(samples)

    print("State assignments")
    print(active)
    print(inactive)

    pcca_object = msm.pcca(5)
    print(msm.metastable_assignments)

    state_sets = []
    for state in range(5):
        print(state)
        microstates = [i for i, x in enumerate(msm.metastable_assignments) if x == state]
        print(microstates)
        state_sets.append(microstates)
    
    print(active)
    for i in active:
        print(msm.metastable_assignments[i])
    print(inactive)
    for i in inactive:
        print(msm.metastable_assignments[i])

    print("Productive binding")
    tpt1, mfpt1, int1, path1 = calc_tpt(msm, inactive, active)
    pickle.dump(tpt1, open("tpt_bound.pkl",'wb'))
    print("MFPT:")
    print(mfpt1)
    print(int1)
    print(path1)

    print(state_sets)

    sets, tpt_cg = tpt1.coarse_grain(state_sets)
    for s in sets:
        print(s)

    flux = tpt_cg.flux
    print(flux)

    print(np.shape(flux))

    print(10*flux/(np.max(flux)))

    paths, caps = tpt_cg.pathways(fraction=0.99)

    print(paths)
    print(caps)

    print(np.sum(caps))

    for i in range(len(sets)):
        for j in range(len(sets)):
            if not i == j:
                _, mfpt, _, _ = calc_tpt(msm, sets[i], sets[j])
                print("MFPT from macrostate %d to %d: %f"%(i, j, mfpt))
                print("Inverse MFPT from macrostate %d to %d: %e"%(i, j, 1/mfpt))

