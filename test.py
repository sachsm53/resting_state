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
#1) IFG
fslmaths /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz -mul 0 -add 1 -roi 46 1 64 1 65  1 0 1 /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_point_new -odt float
fslmaths /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_point_new -kernel sphere 8 -fmean /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_sphere_new -odt float
fslmaths /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_sphere_new -bin /Volumes/MusicProject/Longitudinal_study/Functional/sma_cluster_sphere_new_bin

#2) 


#set locations
datafolder = "/Volumes/MusicProject/Longitudinal_study/Functional/"
#roifile = datafolder + "reg_yearsmus/reg_emo_yearsmus_apos_aneg.gfeat/cope1.feat/cluster_mask_zstat2.nii.gz"
#rois = ["cluster5_sphere.nii.gz"]
rois = ['sma_cluster_sphere_new_bin']

copeList = ["cope2.feat","cope3.feat","cope5.feat"]