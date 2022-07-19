import glob
import numpy as np
import mdtraj as md

targettopfile="stripped.ShHTL7_DOH-covH.prmtop"
#targettopfile="stripped.AtD14_DOH-covH.prmtop"

#pairs = [[3754,353],[3754,371],[3754,378]] #AtD14
pairs = [[3775,336],[3775,355],[3775,364]] #ShHTL7
labels = ["F28","G29","T30"]

for file in glob.glob('./stripped/*xtc'):
    t = md.load(file, top=targettopfile)
    for j in range(len(labels)):

        d = md.compute_distances(t, [pairs[j]])     	# distance between residue pairs on helix T1 and T2 respectively

        n_frames = t.n_frames

        dis = np.empty([n_frames, 1])

        for i in range(n_frames):
          dis[i,0:1]=d[i][0]

        np.save('./analysis/'+file.split('/')[-1]+'_ftr-covD-%s.npy'%labels[j], dis)
