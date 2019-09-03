#!/usr/bin/env python

#Pre-process Resting State data and does an ICA

import os,sys
from commando import commando
from commando import writeToLog
import argparse
import numpy as np
import re
from datetime import datetime
from subprocess import call
from subprocess import check_output
import csv
import shutil
import pandas as pd

#parse command line arguments

parser = argparse.ArgumentParser()

parser.add_argument("--nopre",help="skip all preprocessing steps", action="store_true")
parser.add_argument("--nomc",help="skip motion correction conversion", action="store_true")
parser.add_argument("--noscrub",help="skip motion scrubbing", action="store_true")
parser.add_argument("--noprefeat",help="skip FEAT design preprocessing steps", action="store_true")
parser.add_argument("--noaroma",help="skip Melodic and ICA-AROMA", action="store_true")
parser.add_argument("--noseg",help="skip segmentation", action="store_true")
#parser.add_argument("--noreg",help="skip registration copying", action="store_true")
parser.add_argument("--noconfound",help="skip confound timeseries creation", action="store_true")
parser.add_argument("--noresid",help="skip regressing out confounds", action="store_true")
parser.add_argument("--notf",help="skip temporal filtering", action="store_true")
parser.add_argument("--nopost",help="skip normalizing and downsampling for gICA", action="store_true")

args = parser.parse_args()

#set locations
datafolder = '/Volumes/BCI-1/Matt-Emily/ASD_restingstate/'
if not os.path.exists(datafolder):
	datafolder = '/Volumes/BCI-2/Matt-Emily/ASD_restingstate/'
if not os.path.exists(datafolder):
	datafolder = '/Volumes/BCI/Matt-Emily/ASD_restingstate/'

genericdesign = datafolder + "preprocess_design_fieldmap.fsf"
seeddesign = datafolder + "resting_seed_design.fsf"
maskfolder = datafolder + "ROIs/"
stand_image = "/usr/local/fsl/data/standard/MNI152_T1_2mm_brain" 

#set analysis values
smoothmm = 6	#smoothing sigma fwhm in mm
smoothsigma = smoothmm/2.3548	#convert to sigma
additive = 10000	#value added to distinguish brain from background
brightnessthresh = additive * .75

#logging colors
sectionColor = "\033[94m"
sectionColor2 = "\033[96m"
groupColor = "\033[90m"
mainColor = "\033[92m"
pink = '\033[95m'
yellow = '\033[93m'
red = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'

subcount = 0

subjectList = [elem for elem in os.listdir(datafolder) if "LA" in elem]
excludeList = ["LA022_TD","LA009_ASD","LA10000_Template"] #remove LA009 for motion
subjectList = [elem for elem in subjectList if elem not in excludeList]
subjectList.sort()
#subjectList = ["LA106_ASD"]

for subj in subjectList:
	subject = subj
	subcount = subcount + 1
	
	#define files and folders here, in case a particular step is skipped
	subjfolder = datafolder + subject + '/'
	origdatafile = subjfolder + "NIFTI/RestingState.nii.gz"
	if not os.path.exists(origdatafile):
		origdatafile = subjfolder + "NIFTI/RestingState_7min.nii.gz"

	t1image = subjfolder + "NIFTI/T1_brain.nii.gz"
	wholebrain = subjfolder + "NIFTI/T1.nii.gz"
	radfile1 = subjfolder + "TopupCorrection/Fieldmap_1/my_fieldmap_rads.nii.gz"
	magfile1 = subjfolder + "TopupCorrection/Fieldmap_1/my_fieldmap_mag_brain.nii.gz"
	radfile2 = subjfolder + "NIFTI/TopupCorrection/Fieldmap_1/my_fieldmap_rads.nii.gz"
	magfile2 = subjfolder + "NIFTI/TopupCorrection/Fieldmap_1/my_fieldmap_mag_brain.nii.gz"
	radfile3 = subjfolder + "NIFTI/Fieldmap_1/my_fieldmap_rads.nii.gz"
	magfile3 = subjfolder + "NIFTI/Fieldmap_1/my_fieldmap_mag_brain.nii.gz"

	#Make Resting state preprocessing folder
	rsdir = subjfolder + 'preprocess/'
	if not os.path.exists(rsdir):
		print yellow + 'Making %s folder in %s%s' %(rsdir,subject, mainColor)
		os.mkdir(rsdir)

	#Specific all relevant files for specific subject
	designOutput = rsdir + "rest_preprocess_nonlinear.fsf"
	WMfile = rsdir + "rest_WM.txt"
	CSFfile = rsdir + "rest_CSF.txt"
	scrubout = rsdir + "RestingState_MotionOutliers_FD.txt"
	metric_values_text = rsdir + "scrub_metric_values_FD"
	metric_values_plot = rsdir + "scrub_metric_plot_FD"

	logfile = rsdir + "restanalysis_log.txt"
	reportfile = rsdir + "report.html"

	#Updated on 02-07-18
	designfile = rsdir + "design_noscrub.mat"
	residfile = rsdir + "RestingState_residuals_noscrub.nii.gz"
	tffile = rsdir + "RestingState_pregICA.nii.gz"
	addfile = rsdir + "RestingState_pregICA_add.nii.gz"
	nonlinearwarp = rsdir + "RestingState_pregICA_add_standard.nii.gz"
	downsample = rsdir + "RestingState_pregICA_add_stand_iso.nii.gz"


	#Redo files with different high pass filtereing (9-29-17) and on the not scrubbed data (02-07-18)
	tffile_new = rsdir + "RestingState_resid_tf_noscrub.nii.gz"
	addfile_new = rsdir + "RestingState_resid_tf_noscrub_add.nii.gz"
	nonlinearwarp_new = rsdir + "RestingState_resid_tf_noscrub_add_stand.nii.gz"
	downsample_new = rsdir + "RestingState_resid_tf_noscrub_add_stand_iso.nii.gz"

	#SKIP Subjects where this file exists
	group = subj[-2:]
	if group == 'SD':
		groupnum = 1
	elif group == 'TD':
		groupnum = 2
	elif group == 'CD':
		groupnum = 3


	if os.path.exists(downsample_new):
		print subcount,group,groupnum,downsample_new
		continue

	preprocess_featfolder = rsdir + "rest_preprocess_nonlinear_2.feat"
	mc_abs = preprocess_featfolder + "/mc/prefiltered_func_data_mcf_abs.rms"
	mc_rel = preprocess_featfolder + "/mc/prefiltered_func_data_mcf_rel.rms"
	mcparamfile = preprocess_featfolder + "/mc/prefiltered_func_data_mcf.par"
	outef = preprocess_featfolder + "/example_func.nii.gz"
	outef_brain = preprocess_featfolder + "/example_func_brain.nii.gz" 
	regdir = preprocess_featfolder + "/reg"
	preprocess_output = preprocess_featfolder + "/filtered_func_data.nii.gz"

	# if os.path.exists(regdir):
	# 	print red + "Subject %s has an old, bad feat folder. Be careful%s" %(subject,mainColor)
	# 	#shutil.rmtree(preprocess_featfolder)

	icafolder = preprocess_featfolder + "/ICA_AROMA"
	icafiles = icafolder + '/melodic.ica'
	filteredfile = icafolder + "/denoised_func_data_nonaggr.nii.gz"

	#Everything good? Start the work
	print
	print sectionColor + "WORKING ON SUBJECT: %s%s" %(subject, mainColor)

	#Check to make sure subject has all relevant folders
	if os.path.exists(radfile1):
		radfile = radfile1
		magfile = magfile1
		#print yellow + 'Found radfile at %s%s' %(radfile,mainColor)
	elif os.path.exists(radfile2):
		radfile = radfile2
		magfile = magfile2
		#print yellow + 'Found radfile at %s%s' %(radfile, mainColor)
	elif os.path.exists(radfile3):
		radfile = radfile3
		magfile = magfile3
		#print yellow + 'Found radfile at %s%s' %(radfile, mainColor)
	else: 
		print red + '%s is missing fieldmap data: %s%s' %(subject,radfile1,mainColor)
		continue

	if not os.path.exists(magfile):
		print red + '%s is missing fieldmap data: %s%s' %(subject,magfile,mainColor)
		continue


	#Make sure we have the wholebrain, non skullstripped brain as well
	if not os.path.exists(wholebrain):
		print red + '%s T1 wholebrain image is not found or labeled something different%s' %(subject,mainColor)

	if not os.path.exists(t1image):
		print red + '%s T1_brain image is not found or labeled something different%s' %(subject,mainColor)

	if not os.path.exists(origdatafile):
		print red + '%s Restingstate.nii.gz is not found or labeled something different%s' %(subject,mainColor)


	#Check the number of TRs for resting state
	command = 'fslinfo %s' %(origdatafile)
	origresults = check_output(command,shell=True)
	origreport = origresults.split()
	indx = origreport.index('dim4')
	numtimepoints = origreport[indx + 1]
	if numtimepoints != '150':
		if "7min" in origdatafile:
			print yellow + '%s restingstate file has %s timepoints. Will cut at end' %(subject,numtimepoints)
		else: 
			print red + '%s restingstate file has %s timepoints. Please check. Moving on to next participant' %(subject,numtimepoints)
			continue

	
	if not args.nopre:

		# prepare web page report (only if not created yet)
		if not os.path.exists(tffile_new):
			timestart= datetime.now()
			timestamp = timestart.strftime('%b %d %G %I:%M%p')
			fsldir = os.environ['FSLDIR']
			writeToLog("<html><head><title>Resting State Analysis Report "+subject+"</title><link REL=stylesheet TYPE=text/css href="+fsldir+"/doc/fsl.css></head><body>",reportfile)
			writeToLog("\n<h1>Resting State Analysis for "+subject+"</h1>Processing started at: "+timestamp+"<br><hr><br>",reportfile)
			call("open " + reportfile,shell=True)

		#Check to see if the completed file exists, skip anyone who else it
		# if os.path.exists(tffile_new):
		# 	print yellow + "Preprocessed gICA file already completed for %s. Moving on\n%s"  % (subject,mainColor) 
		# 	continue 

	# 	#----------------------------------------

	# 	# Scrubbing with FD - just to get mFD

	# 	#----------------------------------------				
		if not args.noscrub: 
			if not os.path.exists(metric_values_text):
				print sectionColor2 + " Scrubbing for %s to determine mFD\n%s"  % (subject,mainColor)
				command = "fsl_motion_outliers -i %s -o %s --fd --thresh=%s -s %s -p %s -v" % (origdatafile, scrubout, "0.5", metric_values_text, metric_values_plot)
				commando(command, logfile)

			else: 
				print yellow + "FSL Motion Outliers already completed for %s. Moving on\n%s"  % (subject,mainColor)

		#Make sure scrub file exists and if not, then set the number of scrubbing TRs to zero
		if os.path.exists(metric_values_text):
			fds = np.loadtxt(metric_values_text)
			tr = len(fds)
			mfd = np.mean(fds)

		else:
			print red + "No outliers found for %s. Moving on\n%s"  % (subject,mainColor)
			num_cols = 0


		if not args.noprefeat:

		# 	#----------------------------------------

		# 	# Preprocessing steps in FEAT - smoothing (6mm), motion correction, slice timing correction, registration (nonlinear, warp resolution 10mm), BET, NO TEMPORAL FILTERING

		# 	#----------------------------------------	
			
			print sectionColor + "Preprocessing steps in FEAT ---------------------------%s" %(mainColor)

			if not os.path.exists(preprocess_featfolder):
				command = "sed -e 's/DEFINEINPUT/%s/g' -e 's/DEFINEOUTPUT/%s/g' -e 's/DEFINESTRUCT/%s/g' -e 's/DEFINEVOLUME/%s/g' -e 's/DEFINERAD/%s/g' -e 's/DEFINEMAG/%s/g' %s > %s" % (re.escape(origdatafile),re.escape(preprocess_featfolder),re.escape(t1image), numtimepoints, re.escape(radfile), re.escape(magfile),genericdesign,designOutput)
				commando(command, logfile)
				command = "feat %s" % designOutput
				commando(command, logfile)
			else: 
				print yellow + "FEAT Preprocessing already completed for %s. Moving on\n%s"  % (subject,mainColor)

		#--------------------
		# Report on motion data 
		# -------------------
		#read in the data from the motion report file 1
		filename = preprocess_featfolder + '/report_prestats.html'
		textfile = open(filename,'r')
		filetext = textfile.read()
		textfile.close()

		#find absolute motion
		result_ab1 = re.search('absolute=(.*)mm,',filetext)
		motion_ab1 = result_ab1.groups()[0]

		#find relative motion
		result_rel1= re.search('relative=(.*)mm',filetext)
		motion_rel1 = result_rel1.groups()[0]

		##########
		 #Determine if they moved more than 3mm and print out motion
		##########

		counter = 0
		counter2 = 0
		c1 = 0 
		c2 = 0
		with open(mc_abs) as f:
			for row in csv.reader(f):
				counter = counter + 1
				number = float(row[0])
				if number > 3:
					c1 = c1 + 1
					#print red + "%s has absolute motion greater than 3mm. %f at TR = %d%s" %(subject,number,counter, mainColor)


		with open(mc_rel) as f:
			for row in csv.reader(f):
				counter2 = counter2 + 1
				number = float(row[0])
				if number > 3:
					c2 = c2 + 1
					#print red + "%s has relative motion greater than 3mm. %f at TR = %d%s" %(subject,number,counter2, mainColor)

		#print red + "%s\tTRs: %s\tmean FD: %.2f\tAbs Motion: %s\tRel Motion: %s\tTR greater than 3mm mvmt: %s,%s%s" %(subject,numtimepoints,mfd,motion_ab1,motion_rel1,c1,c2,mainColor)
		print sectionColor2 + "Motion Report: Mean FD: %.2f, Absolute: %s, Relative, %s, Spikes1: %s, Spikes2: %s%s" %(mfd,motion_ab1,motion_rel1,c1,c2,mainColor)

		#----------------------------------------

		# MELODIC ICA and ICA-AROMA

		#----------------------------------------

		if not args.noaroma:

			print sectionColor + "ICA-AROMA ---------------------------" + mainColor	

			#Make sure preprocessing feat ran correctly
			if not os.path.exists(preprocess_output):
				print  red + "Preprocess feat not completed correctly. %s does not exist. Moving on to next subject\n%s"  %(preprocess_output,mainColor)
				continue

			if not os.path.exists(icafolder):
				#print red + "ICA-AROMA has not been completed for %s\n%s"  % (subject,mainColor)
				icalocation = "/Volumes/BCI-1/Matt-Emily/ASD_restingstate/ICA-AROMA/ICA_AROMA.py"
				if not os.path.exists(icalocation):
					icalocation = "/Volumes/BCI/Matt-Emily/ASD_restingstate/ICA-AROMA/ICA_AROMA.py"
				if not os.path.exists(icalocation):
					icalocation = "/Volumes/BCI-2/Matt-Emily/ASD_restingstate/ICA-AROMA/ICA_AROMA.py"
				command = "%s -feat %s -out %s" % (icalocation,preprocess_featfolder,icafolder)
				commando(command,logfile)

			else: 
				print yellow + "ICA-AROMA already completed for %s. Moving on\n%s"  % (subject,mainColor) 

		#Check and make sure it completed properly
		if not os.path.exists(icafiles):
			print red + "%s does not have the completed ICA File%s" %(subject,mainColor)
			continue
		# else:
		# 	print yellow + "%s has the completed ICA File%s" %(subject,mainColor)

		if not args.noseg:

		#----------------------------------------

		# segmentation

		#----------------------------------------

			print sectionColor + "Segmentation ---------------------------" + mainColor

			segfile = '%sT1_brain_seg_0.nii.gz' %(rsdir)

			if not os.path.exists(segfile):
				t1out = rsdir + "T1_brain"
				command = "fast -g -o %s %s" % (t1out,t1image)
				commando(command,logfile)
				command="slicer %s %sT1_brain_seg_0 -a %sCSF.png" % (t1image,rsdir,rsdir)
				commando(command,logfile)
				command="slicer %s %sT1_brain_seg_1 -a %sGM.png" % (t1image,rsdir,rsdir)
				commando(command,logfile)
				command="slicer %s %sT1_brain_seg_2 -a %sWM.png" % (t1image,rsdir,rsdir)
				commando(command,logfile)
				writeToLog("<h2>Segmentation</h2><br>CSF:<br><img src=CSF.png><br><br>White matter:<br><img src=WM.png><br><br>Gray matter:<br><img src=GM.png><br><br><hr>",reportfile)
			else: 
				print yellow + "Segmentation already completed for %s. Moving on\n%s"  % (subject,mainColor) 

		if not args.noconfound:

		#----------------------------------------

		# create confound timeseries

		#----------------------------------------

			# CSF is mprage_brain_seg_0 is and WM is mprage_brain_seg_2

			print sectionColor + "Extracting confound timeseries ---------------------------" + mainColor
			
			if not os.path.exists(CSFfile):			
				command = "flirt -in %sT1_brain_seg_0 -ref %s/example_func.nii.gz -applyxfm -init %s/highres2example_func.mat -interp nearestneighbour -o %srest_CSF.nii.gz" % (rsdir,regdir,regdir,rsdir)
				commando(command,logfile)
				command = "flirt -in %sT1_brain_seg_2 -ref %s/example_func.nii.gz -applyxfm -init %s/highres2example_func.mat -interp nearestneighbour -o %srest_WM.nii.gz" % (rsdir,regdir,regdir,rsdir)
				commando(command,logfile)
				command = "fslmeants -i %s -m %srest_CSF.nii.gz -o %srest_CSF.txt" % (filteredfile,rsdir,rsdir)
				commando(command,logfile)
				command = "fslmeants -i %s -m %srest_WM.nii.gz -o %srest_WM.txt" % (filteredfile,rsdir,rsdir)
				commando(command,logfile)
			else: 
				print yellow + "Confound timeseries already created for %s. Moving on\n%s"  % (subject,mainColor)


		if not args.noresid:

		#----------------------------------------

		# Regress out confounds (do not include scrubbing file): 

		#----------------------------------------

			print sectionColor + "Regressing out confounds  ---------------------------" + mainColor

			#first construct design.mat
			if not os.path.exists(designfile):
				mcparams = np.loadtxt(mcparamfile)
				wmparams = np.loadtxt(WMfile,ndmin=2)
				csfparams = np.loadtxt(CSFfile,ndmin=2)
				allconfounds = np.hstack([mcparams,wmparams,csfparams])
				numconfounds = 8
				
				heights = np.max(allconfounds,0)-np.min(allconfounds,0)
				print red + "Number of total regressors: %s%s" %(numconfounds, mainColor)

				df = open(designfile,'w')
				df.write('/NumWaves\t%d\n' % numconfounds)
				df.write('/NumPoints\t%s\n' % numtimepoints)
				df.write('/PPheights\t')

				for x in range(numconfounds):
					df.write('%f\t' % heights[x])
				df.write('\n\n/Matrix\n')
				numtimepoints = int(numtimepoints)

				for row in range(numtimepoints):
					for column in range(numconfounds):
						df.write('%f\t' % allconfounds[row][column])
					df.write('\n')
				df.close()

			#now use it to regress out confounds

			if not os.path.exists(residfile):
				print yellow + "Regressing out confounds for %s\n%s"  % (subject,mainColor) 
				command = "fsl_glm -i %s --demean -m %s/mask.nii.gz -d %s --out_res=%s" % (filteredfile,preprocess_featfolder,designfile,residfile)
				commando(command,logfile)
			else: 
				print yellow + "Confounds already regressed out alraady for %s. Moving on\n%s"  % (subject,mainColor) 

		if not args.notf:

		#----------------------------------------

		# temporal filtering - Redo everything after this without scrubbing

		#----------------------------------------

		#Check TR here and change if not 2 
			if os.path.exists(residfile):
				command = 'fslinfo %s' %(residfile)
				results = check_output(command,shell=True)
				report = results.split()
				trindx = report.index('pixdim4')
				tr = float(report[trindx+1])

				trorigindx = origreport.index('pixdim4')
				tr_orig = float(origreport[trorigindx+1])
				if tr_orig != tr:
					print red + "%s exists, but TR has been changed from %d to %d%s" %(residfile,tr,tr_orig,mainColor)
					command = "fslcpgeom %s %s" %(origdatafile,residfile)
					print yellow + 'Copying parameters to residfile'
					commando(command,logfile)
			else:
				print red + "Residual file not found. Go back to preprocessing\n%s" %(mainColor)
				continue


			print sectionColor + "Temporal filtering ---------------------------" + mainColor

			if not os.path.exists(tffile_new):
				#hp_sigma = 50 #10 seconds in TRs, pass signal faster than .01 Hz 
				hp_sigma = 25 #10 seconds in TRs, pass signal faster than .01 Hz (CHANGED ON 9/29/17)
				#lp_sigma = 5 #100 seconds in TRs, pass signal slower than .1 Hz (TR = 2seconds * 5 TR = 10 seconds (0.1HZ))
				lp_sigma = -1
				command = "fslmaths %s -bptf %d %d %s -odt float" %(residfile,hp_sigma,lp_sigma,tffile_new)
				commando(command,logfile)
				#print command

				command = "fslmeants -i %s -c 32 32 20 -o %sorigdata_samplets.txt" % (origdatafile,rsdir)
				commando(command,logfile)
				command = "fslmeants -i %s -c 32 32 20 -o %sfiltdata_samplets.txt" % (tffile_new,rsdir)
				commando(command,logfile)

				command="fsl_tsplot -i %sfiltdata_samplets.txt -o %sfilt.png -t 'Bandpass Filtered Data'" % (rsdir,rsdir)
				commando(command,logfile)
				command="fsl_tsplot -i %sorigdata_samplets.txt -o %sorig.png -t 'Unfiltered Data'" % (rsdir,rsdir)
				commando(command,logfile)
				writeToLog("<h2>Temporal Filtering</h2><br><img src=orig.png><br><br><img src=filt.png><br><br><hr>",reportfile)

			else: 
				print yellow + "Temporal filtering already completed for %s. Moving on\n%s"  % (subject,mainColor)


		# #finish report file
		# timeend = datetime.now()
		# endtimetext = datetime.now().strftime('%b %d %G %I:%M%p')
		# elapsed = timeend-timestart
		# elapsedminutes = elapsed.total_seconds()/60
		# elapsedminutesstr = "%.2f" % elapsedminutes
		# writeToLog("<br><hr>Finished at: "+endtimetext+"<br>Elapsed: "+elapsedminutesstr+" minutes<br><br></body></html>",reportfile)

		# print sectionColor + "Done all preprocessing ---------------------------\n\n" + mainColor 

	if not args.nopost:

		#----------------------------------------

		# Prepare data for gICA by normalizing and downsampling

		#----------------------------------------

		print sectionColor + "Postprocessing ---------------------------" + mainColor

		#First check Temporal Filtered file to make sure TR is correct 
		if os.path.exists(tffile_new):
			command = 'fslinfo %s' %(tffile_new)
			results = check_output(command,shell=True)
			report = results.split()
			trindx = report.index('pixdim4')
			tr = float(report[trindx+1])

			trorigindx = origreport.index('pixdim4')
			tr_orig = float(origreport[trorigindx+1])
			if tr_orig != tr:
				print red + "%s exists, but TR has been changed from %d to %d%s" %(tffile_new,tr,tr_orig,mainColor)
		else:
			print red + "Temporal filtered preprocessed file not found. Go back to preprocessing\n%s" %(mainColor)
			continue

		#Cut to time off to be the righ length
		if numtimepoints != '150':
			tffile_7min = rsdir + "RestingState_resid_tf_noscrub_7min.nii.gz"
			os.rename(tffile_new,tffile_7min)
			print yellow + '%s restingstate file was 7minutes. Cutting from %s to 150' %(numtimepoints,subject)
			command = "fslroi %s %s 0 150" % (tffile_7min,tffile_new)
			commando(command,logfile)


		#Add 1000 to raise brain above background
		# os.remove(addfile_new)
		# print red + "Deleting old %s\n%s"  % (addfile_new,mainColor)
		if not os.path.exists(addfile_new): 
			print sectionColor + "Adding background to image ---------------------------" + mainColor
			command = "fslmaths %s -add %d -mul %s/mask %s -odt float" % (tffile_new,additive,preprocess_featfolder,addfile_new)
			commando(command,logfile)
		else: 
			print yellow + "Additing background to preprocessed file already completed. Checking TR...%s"  % (mainColor)

		#Check TR
		command = 'fslinfo %s' %(addfile_new)
		results = check_output(command,shell=True)
		report = results.split()
		trindx = report.index('pixdim4')
		tr_add = float(report[trindx+1])
		if tr_add != 2:
			print red + '%s TR is not 2, its %d!%s' %(addfile_new,tr_add,mainColor)

		# Convert preprocessed file to standard space using nonlinear warp (Applying matrix to warp to standard space)
		# os.remove(nonlinearwarp_new)
		# print red + "Deleting old %s\n%s"  % (nonlinearwarp_new,mainColor)
		if not os.path.exists(nonlinearwarp_new):
			print sectionColor + "Warping to standard ---------------------------" + mainColor
			command = "flirt -in %s -ref %s -applyxfm -init %s/example_func2standard.mat -interp nearestneighbour -o %s" % (addfile_new,stand_image,regdir,nonlinearwarp_new)
			commando(command,logfile)
		else: 
			print yellow + "Standard nonlinear warp already applied to preprocessed file already completed. Checking TR...%s"  % (mainColor)

		#Check TR
		command = 'fslinfo %s' %(nonlinearwarp_new)
		results = check_output(command,shell=True)
		report = results.split()
		trindx = report.index('pixdim4')
		tr_stand = float(report[trindx+1])
		if tr_stand != 2:
			print red + '%s TR is not 2, its %d!%s' %(nonlinearwarp_new,tr_stand,mainColor)

		# Downsample to 4x4x4
		# os.remove(downsample_new)
		# print red + "Deleting old %s\n%s"  % (downsample_new,mainColor)
		if not os.path.exists(downsample_new):
			print sectionColor + "Downsampling to 4x4x4mm ---------------------------" + mainColor
			command = 'flirt -in %s -ref %s -applyisoxfm 4 -interp nearestneighbour -out %s' %(nonlinearwarp_new, stand_image,downsample_new)
			commando(command,logfile)
		else: 
			print yellow + "Downsampling of preprocessed file already completed. Checking TR...%s"  % (mainColor)

		#Check TR
		command = 'fslinfo %s' %(downsample_new)
		results = check_output(command,shell=True)
		report = results.split()
		trindx = report.index('pixdim4')
		tr_ds = float(report[trindx+1])
		if tr_ds != 2:
			print red + '%s TR is not 2, its %d!%s' %(downsample_new,tr_ds,mainColor)
		print 

	# if os.path.exists(downsample_new):
	# 	print '%s\t%d\t%s' %(subject,subcount,downsample_new)




