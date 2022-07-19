#!/bin/bash

#python -m plot_free_energy --ax_files A216-backbone-distance-1 loop-distance-3 --ax_labels S215-backbone\ distance\ \(nm\) D217-H246\ distance\ \(nm\) --ax_lims 0 2 0 2 --savename ShHTL7_dloop-mut-old1_weighted.png --weights msm_275/msm_weights.pkl
python -m plot_free_energy --ax_files A216-backbone-distance-1 loop-distance-3 --ax_labels S215-backbone\ distance\ \(nm\) D217-H246\ distance\ \(nm\) --ax_lims 0 2 0 2 --savename ShHTL7_dloop-mut-old1_weighted_new.png --weights oasis_test/msm_weights.pkl

#python -m plot_free_energy --ax_files salt*distance-3 loop-distance-3 --ax_labels N216-D165\ distance\ \(nm\) D217-H246\ distance\ \(nm\) --ax_lims 0 3 0 2 --savename ShHTL7_dloop-mut-old2_weighted.png --weights msm_275/msm_weights.pkl
python -m plot_free_energy --ax_files salt*distance-3 loop-distance-3 --ax_labels N216-D165\ distance\ \(nm\) D217-H246\ distance\ \(nm\) --ax_lims 0 3 0 2 --savename ShHTL7_dloop-mut-old2_weighted_new.png --weights oasis_test/msm_weights.pkl

#python -m plot_free_energy --ax_files S220-T215 loop-distance-3 --ax_labels M219-S214\ distance\ \(nm\) D217-H246\ \(nm\) --ax_lims 0 2 0 2 --savename ShHTL7_dloop-mut-old3_weighted.png --weights msm_275/msm_weights.pkl
python -m plot_free_energy --ax_files S220-T215 loop-distance-3 --ax_labels M219-S214\ distance\ \(nm\) D217-H246\ distance\ \(nm\) --ax_lims 0 2 0 2 --savename ShHTL7_dloop-mut-old3_weighted_new.png --weights oasis_test/msm_weights.pkl

