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

extract_straight_paths.m (9/9/21)
	Uses the left and right directions of mocap to extract out turning points
	Results in more consistent data since only straight path walking is used
	Save resulting data matrix into both excel and matlab files. 
	Data matrix:
		1-4. arrival index
		5-8. peak index
		9-12. peak mag
		13. coordinate x
		14. coordinate y
		15. heel peak magnitude
		16. toe end time - heel start time
		17. previous x coordinate
		18. previous y coordinate
		19. edge (-1,0,....0,1) to indicate edges of walking segment paths

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

Python
================================================
magnitude_predict.py (8/19/21)
	Predicts the magnitude of heel strike based on data