Gait analysis folders (matlab/python): start after 1st com mtg [8/19/21]
------------------------------------------------
For experiments after amendment for testing different gait analysis. Assumes the sensors are synchronized
by the motion capture system. Raw data is stored in Lab/Analysis/Experiment2 folder. Processed data (.mat workspaces
and excel files) are stored in Lab/Analysis/ProcessedExperiment2

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

distance_btw_coordinates.m (9/14/21)
	Used by extract_straight_paths.m to calculate distance traveled by foot between two heel strikes
	To extract turning points


Python
================================================
magnitude_predict.py (8/19/21)
	Predicts the magnitude of heel strike based on data