#!/bin/bash


model_dir="model"

trn_ft_file="train.X"
trn_lbl_file="train.y"
prop="inv_prop.txt"

# create the model folder
mkdir -p $model_dir

# training
# Reads training features (in $trn_ft_file), training labels (in $trn_lbl_file), and writes FastXML model to $model_dir
./fastXML_train $trn_ft_file $trn_lbl_file $model_dir $prop -T 9 -s 0 -t 10 -b 1.0 -c 1.0 -m 10 -l 10 -q 1

# testing
# Reads test features (in $tst_ft_file), FastXML model (in $model_dir), and writes test label scores to $score_file
# ./parabel_predict $tst_ft_file $model_dir $score_file -t 3

# performance evaluation 
# matlab -nodesktop -nodisplay -r "cd('$PWD'); addpath(genpath('../Tools')); trn_X_Y = read_text_mat('$trn_lbl_file'); tst_X_Y = read_text_mat('$tst_lbl_file'); wts = inv_propensity(trn_X_Y,0.55,1.5); score_mat = read_text_mat('$score_file'); get_all_metrics(score_mat, tst_X_Y, wts); exit;"

