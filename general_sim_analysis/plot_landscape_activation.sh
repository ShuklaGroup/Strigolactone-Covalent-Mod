#!/bin/bash

#python -m plot_free_energy --ax_files T1-T3-distance.npy helical-content-unfolded --ax_labels T1-T3\ distance\ \(nm\) T2\ helical\ content --ax_lims 0 4 0 1 --savename ShHTL7_act_raw_rnd20.png

#python -m plot_free_energy --ax_files T1-T3-distance.npy hinge-contact-1 --ax_labels T1-T3\ distance\ \(nm\) T1-T2\ hinge\ contact --ax_lims 0 4 0 2 --savename ShHTL7_act_raw_hinge_rnd20.png

#python -m plot_free_energy --ax_files T1-T3-distance.npy loop-distance-3 --ax_labels T1-T3\ distance\ \(nm\) D-loop\ distance\ \(nm\) --ax_lims 0 4 0 3 --savename ShHTL7_act_dloop_rnd20.png

python -m plot_free_energy --ax_files T1-T3-distance.npy helical-content-unfolded --ax_labels T1-T3\ distance\ \(nm\) T2\ helical\ content --ax_lims 0 4 0 1 --savename ShHTL7_act_weighted.png --weights oasis4/clusters_350/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy hinge-contact-1 --ax_labels T1-T3\ distance\ \(nm\) T1-T2\ hinge\ contact --ax_lims 0 4 0 2 --savename ShHTL7_act_weighted_hinge.png --weights oasis4/clusters_350/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy loop-distance-3 --ax_labels T1-T3\ distance\ \(nm\) D-loop\ distance\ \(nm\) --ax_lims 0 4 0 3 --savename ShHTL7_act_dloop_weighted.png --weights oasis4/clusters_350/msm_weights.pkl

