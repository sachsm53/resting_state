#!/usr/bin/env python

import os,sys
from commando import commando
from commando import writeToLog
import argparse
import re
from datetime import datetime
from subprocess import call
from subprocess import check_output
import csv as csv 
import numpy as np
from scipy import stats

## Extract the Percent Signal Change for the ROI for each subject (1-4-2017)

#logging colors
sectionColor = "\033[94m"
sectionColor2 = "\033[96m"
mainColor = "\033[92m"

#Create ROI in Standard Space 
#1) Derrfuss et al. 
Inferior frontal junction (6/8/44), left: -40,4,32 (Talairach) --> -41,4,33
ACC/pre-SMA (32/6), right: 2,14,42 (Talairach) --> 2,12,46
ACC/SFG (32/9), left: -2, 36, 26 --> -2, 37, 28

#2) Laird et al 2005
anterior cingulate --> (2,16,41) --> 2,14,45
left IFJ --> -44,6,34 --> -45,6,36
rostal cingulate zone posterior --> 2,16,41 --> 2,14,45
rostral cingulate zone anterior (-3,37,25) --> -3,38,27

#1) IFG
fslmaths /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz -mul 0 -add 1 -roi 46 1 64 1 65  1 0 1 /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_point_new -odt float
fslmaths /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_point_new -kernel sphere 8 -fmean /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_sphere_new -odt float
fslmaths /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_sphere_new -bin /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_sphere_new_bin

#2) pre-SMA/SMA
fslmaths /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz -mul 0 -add 1 -roi 46 1 64 1 65  1 0 1 /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_point_new -odt float
fslmaths /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_point_new -kernel sphere 8 -fmean /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_sphere_new -odt float
fslmaths /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_sphere_new -bin /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_sphere_new_bin



#set locations
datafolder = "/Volumes/MusicProject/Longitudinal_study/Functional/"
#roifile = datafolder + "reg_yearsmus/reg_emo_yearsmus_apos_aneg.gfeat/cope1.feat/cluster_mask_zstat2.nii.gz"
#rois = ["cluster5_sphere.nii.gz"]
rois = ['sma_cluster_sphere_new_bin']

copeList = ["cope2.feat","cope3.feat","cope5.feat"]