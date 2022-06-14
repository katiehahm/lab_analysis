Gait analysis folders (matlab/python): start after 1st com mtg [8/19/21]
------------------------------------------------
For experiments after amendment for testing different gait analysis. Assumes the sensors are synchronized
by the motion capture system. Raw data is stored in Lab/Analysis/Experiment2 folder. Processed data (.mat workspaces
and excel files) are stored in Lab/Analysis/ProcessedExperiment2

Experiment Process
===============================================
Rename files
convert_to_db
extract_straight_paths
hist_steptime
python scripts

Matlab
================================================
convert_to_db.m (8/19/21)
	takes in raw data and converts the information to matlab workspace and an excel spreadsheet for python
	These are the variables:
		coordinates X,Y
		right foot (heel strike time, heel strike magnitude, toe off time)
		left foot (heel strike time, heel strike magnitude, toe off time)
		arrival time index
		peak index, peak magnitude
		peak-to-arrival time
		previous coordinates, if any, if not, (NaN,NaN)
	To use:
	run a section at a time. check that the time lengths are similar
	check that Lheel and Rheel peaks are properly placed
	if not, change findpeaks parameters
	then save file 

find_lpf_thresh.m (8/19/21)
	Maps out the frequency response of a footfall and just noise clip to see where the lpf threshold should be
		
clip_fsr_fromMocap.m (8/30/21)
	clips the fsr data where the 9th channel peaks and dips.
	Start of mocap means value dips < -1
	End of mocap means value rises > -1

clip_pcb_fromMocap.m (8/30/21)
	clips the pcb data where the 6th channel peaks and dips.
	Start of mocap means value rises > 1
	End of mocap means value dips < 1

convertMocap.m (8/30/21)
	converts table form of mocap data to array form

extract_straight_paths.m (9/9/21), heavily rewritten 11/16/21
	Uses the changes in sign in x-directions of mocap to extract out turning points
	Results in more consistent data since only straight path walking is used
	Saves all original variables with turns extracted
	Additional variable walk_episodes saved that indicates start (-1) and stop (1)
	Saved as mat workspace 'filename_extract_straight_paths'

distance_btw_coordinates.m (9/14/21)
	Used by extract_straight_paths.m to calculate distance traveled by foot between two heel strikes
	To extract turning points

findimpacts_fsr_compact.m (9/15/21)
	Same as findimpacts_fsr but only returns one variable that contains all the information
	Impacts:
		1. heel start idx
		2. heel pk index
		3. heel peak mag
		4. toe end idx
		5. toe pk idx
		6. toe pk mag
		7. 1 = right, 0 = left
	adjustable distance between peaks with findpeaks method to find all peaks
	delete peaks that don't have data go to "zero" between
	find actual start of heel strike with aic_pick
	delete 1st peak of both bc usually not clean
	find toe times
	sort, plot

hammering_analysis.m (10/14/21)
	To analyze the hammering data to see if time of arrival is consistent with different 
	types of impacts

hist_steptime.m (10/19/21)
	takes processed data and outputs the histogram of its distribution,
	mean, std, bimodality
	used to analyze step time variability

real_steptime.m (10/21/21)
	to get ground truth of actual step times between left vs right foot
	to compare against the prediction from hist_steptime

steptime_analysis.m (10/24/21)
	steptime analysis other than the histogram plotting
	section 1 looks at direct plotting to observe clustering
	Deletes differences datapoint if the value is > 2x mean step time bc this is outlier

manual_fix_fsr/pcb/mocap.m (11/5/21)
	takes in wrong and right arrays and corrects the features
	manually change values based on the plots

fixing_subj3_fsr_insole.m (11/5/21)
	use this file when fsr conversion failed and need to convert csv file to mat file

trying_snob_steptime.m (11/9/21)
	tried different distributions of mixture models

groundrxnforce_analysis.m (11/10/21)
	(probably) ignore first section "accumulating all forces"
	energy extraction:
		defined a hard noise_thresh cutoff for each sensor 
		window is current arrival idx to next, or to the end of the dataset
		envelope of abs of this window while doing movmean. Find first instance where this dips below thresh
		if none dips below, keep increasing thresh until it does
		sum the window until this dip
		save to variable energy (#impacts x 4)
		this energy var is appended to all processed data files

analysis_11_18_21.m (11/19/21)
	code to analyze data in folder 11/18/21

convert_to_db_wAccel.m (11/22/21)
	similar to convert to db but for data containing accelerometer data
	Used to process data collected on 11/21/21
	assumes fsr = [fsr, accx, accy, accz].

fixing_subj3_fsr_insole.m (11/22/21)
	this file is used when fsr conversion fails
	so there is both .mat and .xlsx file at the same time

find_accel_impacts.m (11/22/21)
	Used to find the peak acceleration for each impact and its indeces

experiment3_convert2db.m (12/13/21)
	First file to run for experiments 12/13/21-12/16/21 (right before winter break)
	Modified from convert to db with accel file

experiment3_analysis.m (12/13/21)
	Part 2 of experiment 3 run
	Real step time, step time GMM, localization, GRF, 

experiment3_makefigures.m (12/21/21)
	scripts to make figures in paper

armax_TA_estimation.m (1/12/22)
	Uses armax coeffs and gets energy, etc. from armax results to estimate TA

experiment4_multitesting.m (2/8/22)
	2 people walking data processing

experiment4_analysis.m (3/1/22)
	to analyze, tbd. contains spectrogram, fft, energy, etc. methods
	most of these don't work. Move onto _footfalldetection.m, these work

crc_analysis.m (3/10/22)
	analyze data collected at crc to see relation between TA and GRF

experiment4_footfalldetection.m (3/24/22)
	used to detect when footfalls happen
	separate doc from exp4_analysis because it was getting too long
	refer to blue notebook to see pseudocode

recursive_stepID.m (3/30/22)
	recursively assigns person ID to list of step times for a single walking segment
	1st round of step times inputs come from experiment4_footfalldetection, where 
	cwt, etc is used to detect step times

experiment4_localizationgrf.m (4/5/22)
	to use after experiment4_footfalldetection.m
	use the step times from previous file to localize and get grf/ta values

delsys_csv2mat.m (4/6/22)
	converts delsys data csv file to mat file
	to use from files after delsys software update

experiment4_allprocessingcompiled.m (4/18/22)
	put the entire processing into sections into this file
	start with raw .mat workspaces and do step time, gmm, etc.

dfs_IDsequence.m (4/20/22)
	uses DFS/binary tree approach to search through all combinations of step times
	outputs a small list of possible step time combinations
	hopefully faster than doing uniqueperms

notrecursive_stepID_limp_whole.m (5/19/22)
	cuts est step time segments so that they are max 22 in length
	processes it using uniqueperms
	parent code (experiment4_allprocessingcompiled.m) uses while loop
	to replace the recursive nature. This runs faster.

https://www.mathworks.com/matlabcentral/fileexchange/24462-wiener-filter-for-noise-reduction-and-speech-enhancement

experiment4_processing1.m (6/10/22)
	initial filtering to sync, extract impacts, ground truth data, eliminate walk edges, wiener filter
experiment4_processing2.m (6/10/22)
	uses cwt to get training data using ground truth
experiment4_processing3.m (6/10/22)
	uses cwt to get testing data, use decision tree to classify impacts, manually clean edges
experiment4_processing4.m (6/10/22)
	store success rate etc, perform GMM, create localization csv
experiment4_processing5.m (6/10/22)
	makes csv for python localization, TA estimation, kmeans

Data
===============================================
Hammering: 9/30
	two takes of hammering on the same spot of the floor (used a side of wrench)
	hammering hard vs soft. Try to see if the aic picker is consistent
	regardless of impact hardness. First take had the impacts too close together, so second take.
	Hammering point was close to sensor3, marked by pink tape on floor.

Experiment2:
	weight is suitcase carry of weight
	box is carrying cardboard box

ProcessedData:
	data as a result of convert_to_db
	two versions - matlab and excel

11_18_21:
	Did testing to figure out a more time limp intervention and also test hammer hits at different magnitudes
	Last dataset in this folder contains accelerometer data from fsr
	Details are at the back of green notebook

12-2-21:
	Testing different knee braces

Python
================================================
GMM_derivation.py (2/28/22)



magnitude_predict.py (8/19/21)
	Predicts the magnitude of heel strike based on data

to run: anaconda prompt and type python main.py

11-30-21_localization_dataset.csv 
	xcoord, ycoord, arrival idx, peak mag, energy, prev x, prev y
	for 11-21-21 dataset


12-1-21_localization_dataset.csv 
	xcoord, curr mag/prev mag, curr energy/prev energy, prev x
	for 11-21-21 dataset

12-1-21_localization_dataset2.csv 
	xcoord, arrival index, curr mag, curr energy, curr mag/prev mag, curr energy/prev energy, prev x
	for 11-21-21 dataset

12-6-21_localization_dataset3.csv
	first column indicates 0 for train, 1 for train
	to use for recursive localization

12-6-21_trackinglocalization_1/2/3/4quartile.csv
	Splits train/test to 0 or 1 value in 1st column
	Uses 25% of data for test, 4 different files
	Used to make predicted location dataset for GRF analysis

12-7-21_grf_dataset_predictedloc.csv
	Contains grf features based on predicted locations
	Doesn't quite work, the rmse is 2.8 rather than 0.6
	But if this works decently, then this is the worse case scenario