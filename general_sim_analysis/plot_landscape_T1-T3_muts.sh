#!/bin/bash

python -m plot_free_energy --ax_files T1-T3-distance.npy H247-G121 --ax_labels T1-T3\ distance\ \(nm\) H246-S119\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename ShHTL7_T1T3_dloop-H247-G121-weighted.png --weights oasis4/clusters_350/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy H247-R173 --ax_labels T1-T3\ distance\ \(nm\) H246-Q172\ distance\ \(nm\) --ax_lims 0 4 0 4 --savename ShHTL7_T1T3_dloop-H247-R173-weighted.png --weights oasis4/clusters_350/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy H247-S168 --ax_labels T1-T3\ distance\ \(nm\) H246-S168\ distance\ \(nm\) --ax_lims 0 4 0 4 --savename ShHTL7_T1T3_dloop-H247-S168-weighted.png --weights oasis4/clusters_350/msm_weights.pkl

#python -m plot_free_energy --ax_files T1-T3-distance.npy G121-Q214 --ax_labels T1-T3\ distance\ \(nm\) S119-Q213\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename ShHTL7_T1T3_dloop-G121-Q214-weighted.png --weights oasis_test/msm_weights.pkl

#python -m plot_free_energy --ax_files T1-T3-distance.npy G121-S97 --ax_labels T1-T3\ distance\ \(nm\) S119-S95\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename ShHTL7_T1T3_dloop-G121-S97-weighted.png --weights oasis_test/msm_weights.pkl

#python -m plot_free_energy --ax_files T1-T3-distance.npy R173-D167 --ax_labels T1-T3\ distance\ \(nm\) Q172-D165\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename ShHTL7_T1T3_dloop-R173-D167-weighted.png --weights oasis_test/msm_weights.pkl

#python -m plot_free_energy --ax_files T1-T3-distance.npy R173-P169 --ax_labels T1-T3\ distance\ \(nm\) Q172-E168\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename ShHTL7_T1T3_dloop-R173-P169-weighted.png --weights oasis_test/msm_weights.pkl
