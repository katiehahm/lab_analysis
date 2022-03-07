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