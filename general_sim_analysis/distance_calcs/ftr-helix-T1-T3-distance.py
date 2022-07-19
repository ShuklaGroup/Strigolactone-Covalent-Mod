import glob
import numpy as np
import mdtraj as md

targettopfile="stripped.AtD14_DOH-covH.prmtop"
#targettopfile="stripped.ShHTL7_DOH-covH.prmtop"

for file in glob.glob('./stripped/*round25*xtc'):
    t = md.load(file, top=targettopfile)
    d = md.compute_contacts(t,[[136,166]], scheme='ca')[0]     	# distance between the two residues on helix T3 and T1 respectively (AtD14)
    #d = md.compute_contacts(t,[[138,169]], scheme='ca')[0]     	# distance between the two residues on helix T3 and T1 respectively (ShHTL7)
    n_frames = t.n_frames

    dis = np.empty([n_frames, 1])

    for i in range(n_frames):
      dis[i,0:1]=d[i][0]

    np.save('./analysis/'+file.split('/')[-1]+'_ftr-helix-T1-T3-distance.npy',dis)
