#!/usr/bin/env python

#process resting state data 

import os,sys
from commando import commando
from commando import writeToLog
import argparse
import numpy as np
import re
from datetime import datetime
from subprocess import call

#parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--nopre",help="skip all preprocessing steps", action="store_true")
parser.add_argument("--nomc",help="skip motion correction conversion", action="store_true")
parser.add_argument("--notf",help="skip temporal filtering", action="store_true")
parser.add_argument("--noseg",help="skip segmentation", action="store_true")
parser.add_argument("--noreg",help="skip registration copying", action="store_true")
parser.add_argument("--noconfound",help="skip confound timeseries creation", action="store_true")
parser.add_argument("--noseed",help="skip creating seed timeseries", action="store_true")
parser.add_argument("--nosmooth",help="skip spatial smoothing", action="store_true")
parser.add_argument("--noresid",help="skip regressing out confounds", action="store_true")
parser.add_argument("--noglm",help="skip glm seed regression", action="store_true")
args = parser.parse_args()

#set locations
maskfolder = "/Volumes/External/Music_training/restingstate/roi/"
datafolder = "/Volumes/External/Music_training/restingstate/group_A/music/"
datafolder_old = "/Volumes/External/Music_training/restingstate/group_A/"
genericdesign = "/Volumes/External/Music_training/restingstate/resting_seed.fsf"

#set analysis values
numtimepoints = 180
numconfounds = 8
smoothmm = 8	#smoothing sigma fwhm in mm
smoothsigma = smoothmm/2.3548	#convert to sigma
additive = 10000	#value added to distinguish brain from background
brightnessthresh = additive * .75

#logging colors
sectionColor = "\033[94m"
mainColor = "\033[92m"

#subjectList = range(24,26)
#subjectList = range(26,33)
#subjectList = [elem for elem in os.listdir(folder) if "." not in elem] 
subjectList = ["6304AD_SS"]

for subj in subjectList:

	subject = "SF%.2d" % subj

	#command = "cp /Volumes/Thalamus/fMRI/narrative/data/%s/mprage_brain.nii.gz /Volumes/Caudate/fMRI/narrative/data/%s/" % (subject,subject)
	#call(command,shell=True)

	#define files and folders here, in case a particular step is skipped
	
	subjfolder = datafolder + subject + "/"
	subjfolder_old = datafolder_old + subject + "/"
	origdatafile = subjfolder + "RestingState.nii.gz"
	restfolder = subjfolder + "restingstate"
	mcfile = restfolder + "/RestingState_mc.nii.gz"
	filteredfile = restfolder + "/RestingState_filtered.nii.gz"
	t1image = subjfolder + "mprage_brain.nii.gz"
	outef = restfolder + "/example_func_rest"
	outef_brain = restfolder + "/example_func_rest_brain" 
	regdir = restfolder + "/reg_rest"
	designfile = restfolder + "/design.mat"
	logfile = restfolder + "/restanalysis_log.txt"
	reportfile = restfolder + "/report.html"

	print sectionColor + "Working on file: %s%s"  % (origdatafile,mainColor)
	if not os.path.exists(restfolder):
		command = "mkdir %s" % restfolder
		call(command,shell=True)

	#prepare web page report
	timestart= datetime.now()
	timestamp = timestart.strftime('%b %d %G %I:%M%p')
	fsldir = os.environ['FSLDIR']
	writeToLog("<html><head><title>Resting State Analysis Report "+subject+"</title><link REL=stylesheet TYPE=text/css href="+fsldir+"/doc/fsl.css></head><body>",reportfile)
	writeToLog("\n<h1>Resting State Analysis for "+subject+"</h1>Processing started at: "+timestamp+"<br><hr><br>",reportfile)
	call("open " + reportfile,shell=True)

	if not args.nopre:
		
		if not args.nomc:
			#----------------------------------------
			# motion correction
			#----------------------------------------
			print sectionColor + "Motion correction ---------------------------" + mainColor

			command = "mcflirt -in %s -o %s -plots -report -rmsrel -rmsabs" % (origdatafile,mcfile)
			commando(command,logfile)

			command="fsl_tsplot -i %s/RestingState_mc.nii.gz.par -t 'MCFLIRT estimated rotations (radians)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o %s/rot.png"  % (restfolder,restfolder)
			commando(command,logfile)
			command="fsl_tsplot -i %s/RestingState_mc.nii.gz.par -t 'MCFLIRT estimated translations (mm)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o %s/trans.png" % (restfolder,restfolder)
			commando(command,logfile)
			command="fsl_tsplot -i %s/RestingState_mc.nii.gz_abs.rms,%s/RestingState_mc.nii.gz_rel.rms -t 'MCFLIRT estimated mean displacement (mm)' -u 1 -w 640 -h 144 -a absolute,relative -o %s/disp.png" % (restfolder,restfolder,restfolder)
			commando(command,logfile)

			absfile = "%s/RestingState_mc.nii.gz_abs_mean.rms" % restfolder
			relfile = "%s/RestingState_mc.nii.gz_rel_mean.rms" % restfolder
			mean_abs = np.loadtxt(absfile)
			mean_rel = np.loadtxt(relfile)

			writeToLog("<h2>Motion Correction</h2><br>Mean absolute motion: <b>"+str(mean_abs)+"</b><br>",reportfile)
			writeToLog("Mean relative motion: <b>"+str(mean_rel)+"</b><br><br><img src=rot.png><br><br><img src=trans.png><br><br><img src=disp.png><br><br><hr>",reportfile)
			

		if not args.notf:
			#----------------------------------------
			# temporal filtering
			#----------------------------------------
			print sectionColor + "Temporal filtering ---------------------------" + mainColor

			hp_sigma = 50 #10 seconds in TRs, pass signal faster than .01 Hz 
			lp_sigma = 5 #100 seconds in TRs, pass signal slower than .1 Hz
			command = "fslmaths %s -bptf %d %d %s -odt float" % (mcfile,hp_sigma,lp_sigma,filteredfile)
			commando(command,logfile)

			command = "fslmeants -i %s -c 32 32 20 -o %s/origdata_samplets.txt" % (origdatafile,restfolder)
			commando(command,logfile)
			command = "fslmeants -i %s -c 32 32 20 -o %s/filtdata_samplets.txt" % (filteredfile,restfolder)
			commando(command,logfile)

			command="fsl_tsplot -i %s/filtdata_samplets.txt -o %s/filt.png -t 'Bandpass Filtered Data'" % (restfolder,restfolder)
			commando(command,logfile)
			command="fsl_tsplot -i %s/origdata_samplets.txt -o %s/orig.png -t 'Unfiltered Data'" % (restfolder,restfolder)
			commando(command,logfile)

			writeToLog("<h2>Temporal Filtering</h2><br><img src=orig.png><br><br><img src=filt.png><br><br><hr>",reportfile)

		if not args.noseg:
			#----------------------------------------
			# segmentation
			#----------------------------------------
			print sectionColor + "Segmentation ---------------------------" + mainColor

			t1out = subjfolder + "mprage_brain"
			command = "fast -g -o %s %s" % (t1out,t1image)
			commando(command,logfile)

			command="slicer %s %smprage_brain_seg_0 -a %s/CSF.png" % (t1image,subjfolder,restfolder)
			commando(command,logfile)
			command="slicer %s %smprage_brain_seg_1 -a %s/GM.png" % (t1image,subjfolder,restfolder)
			commando(command,logfile)
			command="slicer %s %smprage_brain_seg_2 -a %s/WM.png" % (t1image,subjfolder,restfolder)
			commando(command,logfile)

			writeToLog("<h2>Segmentation</h2><br>CSF:<br><img src=CSF.png><br><br>White matter:<br><img src=WM.png><br><br>Gray matter:<br><img src=GM.png><br><br><hr>",reportfile)

		if not args.noreg:
			#----------------------------------------
			# registration
			#----------------------------------------
			#create reference image 
		
			print sectionColor + "Registration ---------------------------" + mainColor

			command = "fslroi %s %s 90 1" % (origdatafile,outef)
			commando(command,logfile)

			#skull strip it
			command = "bet %s %s" % (outef,outef_brain)
			commando(command,logfile)

			#create a brain mask
			command = "fslmaths %s -bin %s/mask.nii.gz" % (outef_brain,restfolder)
			commando(command,logfile)

			#create a folder for registration files
			if not os.path.exists(regdir):
				command = "mkdir %s" % regdir
				commando(command,logfile)

			#compute the registration to highres image
			outmat = regdir + "/example_func2highres.mat"
			command = "flirt -in %s -ref %s -dof 6 -omat %s -o %s/example_func2highres" % (outef_brain,t1image,outmat,regdir)
			commando(command,logfile)

			#create the other neccessary files, borrowing from reg's already done for other runs
			command = "cp %srun1.feat/reg/highres2standard.mat %s/highres2standard.mat" % (subjfolder_old,regdir)
			commando(command,logfile)
			command = "cp %srun1.feat/reg/standard2highres.mat %s/standard2highres.mat" % (subjfolder_old,regdir)
			commando(command,logfile)
			command = "convert_xfm -omat %s/highres2example_func.mat -inverse %s/example_func2highres.mat" % (regdir,regdir)
			commando(command,logfile)
			command = "convert_xfm -omat %s/example_func2standard.mat -concat %s/highres2standard.mat %s/example_func2highres.mat" % (regdir,regdir,regdir)
			commando(command,logfile)
			command = "convert_xfm -omat %s/standard2example_func.mat -inverse %s/example_func2standard.mat" % (regdir,regdir)
			commando(command,logfile)
			command = "cp %s %s/highres.nii.gz" % (t1image,regdir)
			commando(command,logfile)
			command = "ln -s $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz %s/standard.nii.gz" % (regdir)
			commando(command,logfile)

			command="slicer %s/highres %s/example_func2highres -a %s/reg.png" % (regdir,regdir,restfolder)
			commando(command,logfile)

			writeToLog("<h2>Registration</h2><br><img src=reg.png><br><br><hr>",reportfile)
			
		if not args.noconfound:
			#----------------------------------------
			# create confound timeseries
			#----------------------------------------
			# mprage_brain_seg_0 is CSF
			# mprage_brain_seg_2 is WM

			print sectionColor + "Extracting confound timeseries ---------------------------" + mainColor

			command = "flirt -in %smprage_brain_seg_0 -ref %s/example_func_rest -applyxfm -init %s/highres2example_func.mat -interp nearestneighbour -o %s/rest_CSF.nii.gz" % (subjfolder,restfolder,regdir,restfolder)
			commando(command,logfile)
			command = "flirt -in %smprage_brain_seg_2 -ref %s/example_func_rest -applyxfm -init %s/highres2example_func.mat -interp nearestneighbour -o %s/rest_WM.nii.gz" % (subjfolder,restfolder,regdir,restfolder)
			commando(command,logfile)
			command = "fslmeants -i %s -m %s/rest_CSF.nii.gz -o %s/rest_CSF.txt" % (filteredfile,restfolder,restfolder)
			commando(command,logfile)
			command = "fslmeants -i %s -m %s/rest_WM.nii.gz -o %s/rest_WM.txt" % (filteredfile,restfolder,restfolder)
			commando(command,logfile)


		if not args.noresid:
			#----------------------------------------
			# regress out confounds
			#----------------------------------------
		
			print sectionColor + "Regressing out confounds  ---------------------------" + mainColor

			#first construct design.mat
			mcparamfile = mcfile + ".par"
			WMfile = restfolder + "/rest_WM.txt"
			CSFfile = restfolder + "/rest_CSF.txt"

			mcparams = np.loadtxt(mcparamfile)
			wmparams = np.loadtxt(WMfile,ndmin=2)
			csfparams = np.loadtxt(CSFfile,ndmin=2)
			allconfounds = np.hstack([mcparams,wmparams,csfparams])

			heights = np.max(allconfounds,0)-np.min(allconfounds,0)

			df = open(designfile,'w')
			df.write('/NumWaves\t%d\n' % numconfounds)
			df.write('/NumPoints\t%d\n' % numtimepoints)
			df.write('/PPheights\t')
			for x in range(numconfounds):
				df.write('%f\t' % heights[x])
			df.write('\n\n/Matrix\n')
			for row in range(numtimepoints):
				for column in range(numconfounds):
					df.write('%f\t' % allconfounds[row][column])
				df.write('\n')

			df.close()

			#now use it to regress out confounds
			command = "fsl_glm -i %s --demean -m %s/mask.nii.gz -d %s/design.mat --out_res=%s/RestingState_residuals" % (filteredfile,restfolder,restfolder,restfolder)
			commando(command,logfile)

			#add 1000 to raise brain above background
			command = "fslmaths %s/RestingState_residuals -add %d -mul %s/mask %s/RestingState_residuals -odt float" % (restfolder,additive,restfolder,restfolder)
			commando(command,logfile)

		if not args.nosmooth:
			#----------------------------------------
			# spatial smoothing
			#----------------------------------------	

			print sectionColor + "Smoothing ---------------------------" + mainColor

			command = "fslmaths %s/RestingState_residuals -Tmean %s/mean_func" % (restfolder,restfolder)
			commando(command,logfile)

			command = "susan %s/RestingState_residuals %.3f %.3f 3 1 1 %s/mean_func %.3f %s/RestingState_residuals_smoothed" % (restfolder,brightnessthresh,smoothsigma,restfolder,brightnessthresh,restfolder)
			commando(command,logfile)

			command = "fslmaths %s/RestingState_residuals_smoothed -mas %s/mask %s/RestingState_residuals_smoothed" % (restfolder,restfolder,restfolder)
			commando(command,logfile)

		if not args.noseed:
			#----------------------------------------
			# create seed timeseries
			#----------------------------------------

			print sectionColor + "Extracting seed timeseries ---------------------------" + mainColor

			pmc_mask = maskfolder + "PCC_fox_7.5mm.nii.gz"
			mpfc_mask = maskfolder + "MPFC_fox_7.5mm.nii.gz"
			lp_mask = maskfolder + "LP_fox_7.5mm.nii.gz"
			pmc_mask_out = restfolder + "/pmc_mask.nii.gz"
			mpfc_mask_out = restfolder + "/mpfc_mask.nii.gz"
			lp_mask_out = restfolder + "/lp_mask.nii.gz"

			command = "flirt -in %s -ref %s/example_func_rest -applyxfm -init %s/standard2example_func.mat -o %s" % (pmc_mask,restfolder,regdir,pmc_mask_out)
			commando(command,logfile)
			command = "flirt -in %s -ref %s/example_func_rest -applyxfm -init %s/standard2example_func.mat -o %s" % (mpfc_mask,restfolder,regdir,mpfc_mask_out)
			commando(command,logfile)
			command = "flirt -in %s -ref %s/example_func_rest -applyxfm -init %s/standard2example_func.mat -o %s" % (lp_mask,restfolder,regdir,lp_mask_out)
			commando(command,logfile)

			command = "fslmeants -i %s/RestingState_residuals_smoothed -m %s -o %s/pmc_ts.txt" % (restfolder,pmc_mask_out,restfolder)
			commando(command,logfile)
			command = "fslmeants -i %s/RestingState_residuals_smoothed -m %s -o %s/mpfc_ts.txt" % (restfolder,mpfc_mask_out,restfolder)
			commando(command,logfile)
			command = "fslmeants -i %s/RestingState_residuals_smoothed -m %s -o %s/lp_ts.txt" % (restfolder,lp_mask_out,restfolder)
			commando(command,logfile)

			command = "slicer %s/example_func_rest %s/pmc_mask -x .5 %s/pmc.png" % (restfolder,restfolder,restfolder)
			commando(command,logfile)
			command = "slicer %s/example_func_rest %s/mpfc_mask -x .5 %s/mpfc.png" % (restfolder,restfolder,restfolder)
			commando(command,logfile)
			command = "slicer %s/example_func_rest %s/lp_mask -z .7 %s/lp.png" % (restfolder,restfolder,restfolder)
			commando(command,logfile)

			writeToLog("<h2>Masks</h2><br>PMC:<br><img src=pmc.png><br><br>MPFC:<br><img src=mpfc.png><br><br>LP:<br><img src=lp.png><br><br><hr>",reportfile)


	if not args.noglm:
		#----------------------------------------
		#seed correlation GLM
		#----------------------------------------

		print sectionColor + "GLM analysis ---------------------------" + mainColor

		#seedlist = ["pmc","mpfc","lp"]
		seedlist  = ["pmc"]
		for seed in seedlist:
			inputfile = "%s/RestingState_residuals_smoothed" % restfolder
			outputfile = "%s/%s_seed.fsf" % (restfolder,seed)
			command = "sed -e 's/DEFINESUBJECT/%s/g' -e 's/DEFINEINPUT/%s/g' -e 's/DEFINESEED/%s/g' %s > %s" % (subject,re.escape(inputfile),seed,genericdesign,outputfile)
			commando(command,logfile)

			writeToLog("<h2>Feat analysis</h2><br><a href="+restfolder+"/"+seed+"_seed.feat/report.html>Feat results for "+seed+"</a><br>",reportfile)

			command = "feat %s" % outputfile
			commando(command,logfile)

			command = "cp -r %s %s/%s_seed.feat/reg" % (regdir,restfolder,seed)
			commando(command,logfile)

	#finish report file
	timeend = datetime.now()
	endtimetext = datetime.now().strftime('%b %d %G %I:%M%p')
	elapsed = timeend-timestart
	elapsedminutes = elapsed.total_seconds()/60
	elapsedminutesstr = "%.2f" % elapsedminutes
	writeToLog("<br><hr>Finished at: "+endtimetext+"<br>Elapsed: "+elapsedminutesstr+" minutes<br><br></body></html>",reportfile)
	
	print sectionColor + "Done ---------------------------\n\n" + mainColor 

