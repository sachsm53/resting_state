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

from subprocess import check_output



#parse command line arguments

parser = argparse.ArgumentParser()

parser.add_argument("--nopre",help="skip all preprocessing steps", action="store_true")

parser.add_argument("--notf",help="skip temporal filtering", action="store_true")

parser.add_argument("--noseg",help="skip segmentation", action="store_true")

parser.add_argument("--noreg",help="skip registration copying", action="store_true")

parser.add_argument("--noconfound",help="skip confound timeseries creation", action="store_true")

parser.add_argument("--noseed",help="skip creating seed timeseries", action="store_true")

parser.add_argument("--noresid",help="skip regressing out confounds", action="store_true")

parser.add_argument("--noglm",help="skip glm seed regression", action="store_true")

args = parser.parse_args()



#set locations

maskfolder = "/Volumes/External/Music_training/restingstate/roi/"

datafolder = "/Volumes/External/Music_training/restingstate/groupA_Y3/Music/"

genericdesign = "/Volumes/External/Music_training/restingstate/resting_seed_music.fsf"



#set analysis values

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

#subjectList = [elem for elem in os.listdir(datafolder) if "." not in elem] 

#excludeList = ["c6747AD-HS", "c6618AD_JG", "c6673AD-RG", "c6720ad-bp", "c6757AD-RL", "c6765ad-jc", "c6771AD-AV", "m6373AD-RM", "m6529AD-HBJ", "s6659ad-am", "s6814AD-GH", "c6783AD-RL", "m6304AD_SS", "m6325AD-AS", "m6348AD-OM", "c6892AD-IB", "c6901ad-sh"]

#subjectList = [elem for elem in subjectList if elem not in excludeList]

#subjectList = ["02_7924AD_HS","03_7926AD_RG","05_7980AD_BP","07_7993AD_JS"] #Control subjects, clipped and ICA2 removed
#subjectList = ["03_7926AD_RG","05_7980AD_BP","07_7993AD_JS"] 
#subjectList = ['09_7879AD_AS','10_7867AD_MV'] #Music subjects, clipped and ICA2 removed
#subjectList = ['11_7964AD_RL'] # control subject, clipped, no ICA2 removed
#subjectList = ['15_7919AD_NR'] # music subject, clipped, no ICA2 removed
#subjectList = ['08_7949AD_JC','19_8042AD_LE'] # control subject, not clipped, ICA1 removed
subjectList = ['01_7878AD_AP','04_7861AD_MA','11_7864AD_MG']


print subjectList

for subj in subjectList:



	subject = subj



	#command = "cp /Volumes/Thalamus/fMRI/narrative/data/%s/mprage_brain.nii.gz /Volumes/Caudate/fMRI/narrative/data/%s/" % (subject,subject)

	#call(command,shell=True)



	#define files and folders here, in case a particular step is skipped


	subjfolder = datafolder + subject + "/"

	#origdatafile = subjfolder + "RestingState_clip_denoised.nii.gz"

	origdatafile = subjfolder + "RestingState_denoised.nii.gz"

	restfolder = subjfolder + "restingstate_clip_ica"

	filteredfile = restfolder + "/RestingState_tf.nii.gz"

	t1image = subjfolder + "mprage_brain.nii.gz"

	icafolder = subjfolder + "ICA_raw.ica"

	outef = icafolder + "/example_func"

	outef_brain = icafolder + "/example_func_brain" 

	regdir = icafolder + "/reg"

	designfile = restfolder + "/design.mat"

	logfile = restfolder + "/restanalysis_log.txt"

	reportfile = restfolder + "/report.html"

	command = 'fslinfo %s' %(origdatafile)

	results = check_output(command,shell=True)

	numtimepoints = results.split()[9]

	print numtimepoints



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


		if not args.notf:

			#----------------------------------------

			# temporal filtering

			#----------------------------------------

			print sectionColor + "Temporal filtering ---------------------------" + mainColor



			hp_sigma = 50 #10 seconds in TRs, pass signal faster than .01 Hz 

			lp_sigma = 5 #100 seconds in TRs, pass signal slower than .1 Hz

			command = "fslmaths %s -bptf %d %d %s -odt float" % (origdatafile,hp_sigma,lp_sigma,filteredfile)

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

			command = 'fslinfo %s/RestingState_tf.nii.gz' %(restfolder)
			results = check_output(command,shell=True)
			numtimepoints = results.split()[9]
			print numtimepoints

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

		

			#skull strip example_func

			command = "bet %s %s" % (outef,outef_brain)

			commando(command,logfile)


			#create a brain mask

			command = "fslmaths %s -bin %s/mask.nii.gz" % (outef_brain,restfolder)

			commando(command,logfile)


			print sectionColor + "Inverse Warping ---------------------------" + mainColor


			#create a folder for registration files

			if not os.path.exists(regdir):

				command = "mkdir %s" % regdir

				commando(command,logfile)



			#create the other neccessary files, borrowing from reg's already done for other runs

			#compute the registration to highres image

			outmat = regdir + "/example_func2highres.mat"

			command = "flirt -in %s -ref %s -dof 6 -omat %s -o %s/example_func2highres" % (outef_brain,t1image,outmat,regdir)

			commando(command,logfile)


			#Do highres to standard registration
			outmat = regdir + "/highres2standard.mat"

			command = "flirt -in %s -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain -dof 12 -omat %s -o %s/highres2standard" % (t1image,outmat,regdir)

			commando(command,logfile)
		

			command = "convert_xfm -omat %s/standard2highres.mat -inverse %s/highres2standard.mat" % (regdir,regdir)
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



			command = "flirt -in %smprage_brain_seg_0 -ref %s/example_func -applyxfm -init %s/highres2example_func.mat -interp nearestneighbour -o %s/rest_CSF.nii.gz" % (subjfolder,icafolder,regdir,restfolder)

			commando(command,logfile)

			command = "flirt -in %smprage_brain_seg_2 -ref %s/example_func -applyxfm -init %s/highres2example_func.mat -interp nearestneighbour -o %s/rest_WM.nii.gz" % (subjfolder,icafolder,regdir,restfolder)

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

			mcparamfile = os.path.join(subjfolder, "ICA_raw.ica/mc/","prefiltered_func_data_mcf.par")
			#mcparamfile = os.path.join(subjfolder, "RestingState_mc.par")

			WMfile = restfolder + "/rest_WM.txt"

			CSFfile = restfolder + "/rest_CSF.txt"



			mcparams = np.loadtxt(mcparamfile)

			wmparams = np.loadtxt(WMfile,ndmin=2)

			csfparams = np.loadtxt(CSFfile,ndmin=2)

			allconfounds = np.hstack([mcparams,wmparams,csfparams])



			heights = np.max(allconfounds,0)-np.min(allconfounds,0)



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

			command = "fsl_glm -i %s --demean -m %s/mask.nii.gz -d %s/design.mat --out_res=%s/RestingState_residuals.nii.gz" % (filteredfile,restfolder,restfolder,restfolder)

			commando(command,logfile)



			#add 1000 to raise brain above background

			command = "fslmaths %s/RestingState_residuals.nii.gz -add %d -mul %s/mask %s/RestingState_residuals.nii.gz -odt float" % (restfolder,additive,restfolder,restfolder)

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



			command = "flirt -in %s -ref %s/example_func -applyxfm -init %s/standard2example_func.mat -o %s" % (pmc_mask,icafolder,regdir,pmc_mask_out)

			commando(command,logfile)

			command = "flirt -in %s -ref %s/example_func -applyxfm -init %s/standard2example_func.mat -o %s" % (mpfc_mask,icafolder,regdir,mpfc_mask_out)

			commando(command,logfile)

			command = "flirt -in %s -ref %s/example_func -applyxfm -init %s/standard2example_func.mat -o %s" % (lp_mask,icafolder,regdir,lp_mask_out)

			commando(command,logfile)



			command = "fslmeants -i %s/RestingState_residuals.nii.gz -m %s -o %s/pmc_ts.txt" % (restfolder,pmc_mask_out,restfolder)

			commando(command,logfile)

			command = "fslmeants -i %s/RestingState_residuals.nii.gz -m %s -o %s/mpfc_ts.txt" % (restfolder,mpfc_mask_out,restfolder)

			commando(command,logfile)

			command = "fslmeants -i %s/RestingState_residuals.nii.gz -m %s -o %s/lp_ts.txt" % (restfolder,lp_mask_out,restfolder)

			commando(command,logfile)



			command = "slicer %s/example_func %s/pmc_mask -x .5 %s/pmc.png" % (icafolder,restfolder,restfolder)

			commando(command,logfile)

			command = "slicer %s/example_func %s/mpfc_mask -x .5 %s/mpfc.png" % (icafolder,restfolder,restfolder)

			commando(command,logfile)

			command = "slicer %s/example_func %s/lp_mask -z .7 %s/lp.png" % (icafolder,restfolder,restfolder)

			commando(command,logfile)



			writeToLog("<h2>Masks</h2><br>PMC:<br><img src=pmc.png><br><br>MPFC:<br><img src=mpfc.png><br><br>LP:<br><img src=lp.png><br><br><hr>",reportfile)





	if not args.noglm:

		#----------------------------------------

		#seed correlation GLM

		#----------------------------------------



		print sectionColor + "GLM analysis ---------------------------" + mainColor



		seedlist = ["pmc","mpfc","lp"]

		#seedlist  = ["pmc"]

		for seed in seedlist:

			inputfile = "%s/RestingState_residuals.nii.gz" % restfolder

			outputfile = "%s/%s_seed.fsf" % (restfolder,seed)

			command = "sed -e 's/DEFINESUBJECT/%s/g' -e 's/DEFINEINPUT/%s/g' -e 's/DEFINEVOLUME/%s/g' -e 's/DEFINESEED/%s/g' %s > %s" % (subject,re.escape(inputfile),numtimepoints,seed,genericdesign,outputfile)

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



