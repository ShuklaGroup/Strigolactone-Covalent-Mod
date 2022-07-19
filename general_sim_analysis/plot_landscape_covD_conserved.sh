#!/bin/bash

python -m plot_free_energy --ax_files T1-T3-distance.npy covD-F28 --ax_labels T1-T3\ distance\ \(nm\) D-ring-Y26\ distance\ \(nm\) --ax_lims 0 4 0 3 --savename ShHTL7_T1T3-covD-F28-weighted.png --weights oasis4/clusters_350/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy covD-G29 --ax_labels T1-T3\ distance\ \(nm\) D-ring-G27\ distance\ \(nm\) --ax_lims 0 4 0 3 --savename ShHTL7_T1T3-covD-G29-weighted.png --weights oasis4/clusters_350/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy covD-T30 --ax_labels T1-T3\ distance\ \(nm\) D-ring-T28\ distance\ \(nm\) --ax_lims 0 4 0 3 --savename ShHTL7_T1T3-covD-T30-weighted.png --weights oasis4/clusters_350/msm_weights.pkl

