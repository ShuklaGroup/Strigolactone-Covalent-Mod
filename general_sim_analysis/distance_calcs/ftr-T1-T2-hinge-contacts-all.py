import glob
import numpy as np
import mdtraj as md

#targettopfile="stripped.ShHTL7_DOH-covH.prmtop"
targettopfile="stripped.AtD14_DOH-covH.prmtop"

pairs = [[144,148],[145,149],[146,150]] #AtD14
#pairs = [[146,150],[147,151],[148,152]] #ShHTL7
labels = ["1","2","3"]

for file in glob.glob('./stripped/*round25*xtc'):
    t = md.load(file, top=targettopfile)
    for j in range(len(labels)):

        d = md.compute_contacts(t, [pairs[j]], scheme='ca')[0]     	# distance between residue pairs on helix T1 and T2 respectively
        n_frames = t.n_frames

        dis = np.empty([n_frames, 1])

        for i in range(n_frames):
          dis[i,0:1]=d[i][0]

        np.save('./analysis/'+file.split('/')[-1]+'_ftr-helix-T1-T2-hinge-contact-%s.npy'%labels[j], dis)
