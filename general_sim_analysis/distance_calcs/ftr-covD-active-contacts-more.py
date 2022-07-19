import glob
import numpy as np
import mdtraj as md

targettopfile="stripped.ShHTL7_DOH-covH.prmtop"
#targettopfile="stripped.AtD14_DOH-covH.prmtop"

#pairs = [[3754,2463],[3754,3351],[3754,2678],[3754,2623],[3754,3335],[3754,344],[3754,3769],[3754,1905]] #AtD14
pairs = [[3775,2464],[3775,3363],[3775,2705],[3775,2649],[3775,3350],[3775,327],[3775,3790],[3775,1884]] #ShHTL7
labels = ["A163","S220","T178","F175","V219","G27","L248","F126"]

for file in glob.glob('./stripped/*xtc'):
    t = md.load(file, top=targettopfile)
    for j in range(len(labels)):

        d = md.compute_distances(t, [pairs[j]])     	# distance between residue pairs on helix T1 and T2 respectively

        n_frames = t.n_frames

        dis = np.empty([n_frames, 1])

        for i in range(n_frames):
          dis[i,0:1]=d[i][0]

        np.save('./analysis/'+file.split('/')[-1]+'_ftr-covD-%s.npy'%labels[j], dis)
