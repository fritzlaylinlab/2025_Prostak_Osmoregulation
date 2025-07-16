## script to process cell wall and contractile vacuole count data 
## script will process data and return results in blinded and unblinded format
## any file names cannot have any spaces


## This script will take a blinded .csv of 
## linescan, total cell intensity, and pump counts per cell per time point and
## will normalized the line scan by % along the line, and background subtract, divide the intensity values
## will do the same intensity normalization with the total cell intensity 
## will also consolidate all the pump data per cell
## will determine the time it takes for a cell to stop pumping after the cell wall reaches a given % above background


import pandas as pd
from datetime import datetime
import re

##########################################################################################
######################### BEGIN USER DEFINED INPUT #######################################

## No change needed: sets current date of when the script is run
current_date = datetime.now().strftime('%Y-%m-%d')


## the file containing cell wall and pump data for processing
## input data file MUST be a .csv with the following columns:
	## file	exp_date	cell	Total_pumps	IGNORE	area	cfw_mean	cfw_min	cfw_max	distance_um	cfw_intensity
data_inFile = "/Users/laboratory/Documents/Fritz_Laylin_Lab/Bd_osmoreg/cv-turgor_transition/analyzed_data/blinded_spreadsheets/20250311_all_updated_wall-pumps.csv"

## the file containing the average background intensity per file
## input background intensity file MUST be a .csv with the following columns:
	## blinded_file	IGNORE	area	cfw_mean	cfw_min	cfw_max
background_file = "/Users/laboratory/Documents/Fritz_Laylin_Lab/Bd_osmoreg/cv-turgor_transition/analyzed_data/blinded_spreadsheets/Updated-Cell_wall_and_pumps_backgrounds.csv"

## the file containing the blinding key
## input blinding key MUST be a .csv with the following columns:
	##OG_file	blinded
##OG_file MUST be in the following format:
	##[YYYMMDD]_[string]_##min[string]_###_[string]
key_file = "/Users/laboratory/Documents/Fritz_Laylin_Lab/Bd_osmoreg/cv-turgor_transition/analyzed_data/blinding_keys/20250311_all_blindingKeys.csv"

## directory where the processed blinded data will go
## DO NOT include final /
blinded_output_directory = "/Users/laboratory/Documents/Fritz_Laylin_Lab/Bd_osmoreg/cv-turgor_transition/analyzed_data/processed/blinded_still"

## directory where the processed and unblinded data will go
## DO NOT include final /
unblinded_output_directory = "/Users/laboratory/Documents/Fritz_Laylin_Lab/Bd_osmoreg/cv-turgor_transition/analyzed_data/processed"

## the names of the output files
## full path will be "{output_directory}_{current_date}_{outfile}.csv"
## blinded files will have _B between {outfile} and .csv
lineScan_outfile = "norm_lineScan"
pump_int_outfile = "norm_pump_int"
avg_outfile = "rep_avgerages"
pumpStopTiming_outfile = "pumpStopTiming"
counts_outfile = "timing_counts"

## define the cutoff for determining at which point the cell is considered to have a cell wall
## this value will be used to find the first time point in which total cell intensity
## is equal to or brighter than this percentage value of the background
threshold_pct = 20

## define the cutoff for pumps/min 
## this value will be used to find the first time point in which pumps/min is equal to or fewer than this in a cell
## the time point determined here will be used to calculate the time from 
## when a cell is considered to have a cell wall until the cell has equal to or fewer this pump/min than this value
theshold_pumps = 0

timing_column_name = f'''time_until_below{theshold_pumps}'''

## the highest time point in the data, used for when a cell does not reach below the defined value pumps/min in the time frame
max_time = 40


######################### END USER DEFINED INPUT #########################################
##########################################################################################

## initialize list to store values for making the mapping dataframe
mapping_data = []

## process the blinding key into a dictionary for mapping timepoint and experiment date info 
with open(key_file, "r") as f:
	for line in f:
		line = line.rstrip()
		## store the line in a list using the "," separator
		line_list = line.split(",")

		# Skip header line
		if "OG_file" in line:
			continue

		# Extract relevant values from the line
		og_name = line_list[0]
		blinded_name = line_list[1]

		# Extract time and date using regex
		time_match = re.search(r'_(\d+)min', og_name)
		date_match = re.search(r'(\d+)', og_name)

		# Ensure matches were found before accessing groups
		time = time_match.group(1) if time_match else None
		date = date_match.group(1) if date_match else None

		# Append extracted data to the mapping list
		mapping_data.append({
			'file': blinded_name,
			'time': time,
			'date': date
		})


## Convert list of dictionaries to DataFrame
mapping_df = pd.DataFrame(mapping_data)
print(f'''
Mapping key populated
''')


## Display the resulting DataFrame
# print(mapping_df)


## import the data file as a pandas dataframe
data_df = pd.read_csv(data_inFile)
# print(data_df)

## import the background data file as a pandas dataframe
bckgd_df =  pd.read_csv(background_file)
# print(bckgd_df)

## Merge data_df with bckgd_df based on file/blinded_file columns
data_df = data_df.merge(bckgd_df[['blinded_file', 'cfw_mean']], 
				left_on='file', 
				right_on='blinded_file', 
				how='left')

## Rename the new cfw_mean column to 'background'
data_df.rename(columns={'cfw_mean_y': 'background'}, inplace=True)

# Drop the redundant 'blinded_file' column
data_df.drop(columns=['blinded_file'], inplace=True)

## Rename existing cfw_mean column to avoid confusion
data_df.rename(columns={'cfw_mean_x': 'cfw_mean'}, inplace=True)

# print(data_df)
print(f'''
Background intensities added
''')


## group the data by file then cell
## each group represents one cell in a given file
indv_cells = data_df.groupby(["file", "cell"])


## create a new dataframe to have only normalized line scan data
norm_line_df = pd.DataFrame()

## create a new tuple to add pump and intensity data, once for each cell per image
int_pump_data = []


## normalize data and pull values for each cell
for name, group in indv_cells:

	file = name[0]
	cell = name[1]
	
	## get the background intensity value for each image
	bckgd = group["background"].iloc[0]
	
	
	## get the line scan line length
	max_dist = group["distance_um"].max()
# 	print(max_dist)
	
	## find the max intensity value along the line
	max_int = group["cfw_intensity"].max()
# 	print(max_int)
	
	
	## normalize the intensities along the line as a percentage of the max
	## add normalized intensity values to a new column
	norm_line_pct = (group["cfw_intensity"]/max_int)*100
	group["norm_cfw_pct"] = norm_line_pct
	
	## normalize the intensities along the line by subtracting or dividing the background
	## add normalized intensity values to table
	norm_line_sub = group["cfw_intensity"] - bckgd
	norm_line_div = group["cfw_intensity"]/bckgd
	group["norm_cfw_sub"] = norm_line_sub
	group["norm_cfw_div"] = norm_line_div
	
	## normalize the distance along the line as a percentage of the length
	## add normalized distance to a new column
	norm_dist = (group["distance_um"]/max_dist)*100
	group["norm_dist"] = norm_dist
	
	## take only the given columns for line scan data and store them in a new dataframe
	line_info = group[["file", "cell", "norm_dist", "cfw_intensity", "norm_cfw_sub", "norm_cfw_div", "norm_cfw_pct"]]
	
	## add the normalized line scan info for each cell per file to the normalized line scan dataframe
	norm_line_df = pd.concat([norm_line_df, line_info], ignore_index = True)
	
	## get the pump count for each cell per image
	pumps = group["Total_pumps"].iloc[0]
	
	## get the average intensity of each cell per image
	mean_int = group["cfw_mean"].iloc[0]
	

	## subtract or divide the background from the mean intensity
	## add normalized intensity to new column
	norm_int_sub = mean_int - bckgd
	norm_int_div = mean_int/bckgd
	group["norm_int_sub"] = norm_int_sub
	group["norm_int_div"] = norm_int_div 
	
	
	## calculate what percentage brighter the average cell intensity is than background
	pct_brighter = (norm_int_sub/ bckgd)*100
	group["pct_brighter"] = pct_brighter
	
	## take only the given columns for average cell intensity and pump counts
	## store them in a tuple
	int_pump_data.append((file, cell, pumps, mean_int, norm_int_sub, norm_int_div, pct_brighter))
	


# Convert the collected pump and intensity data into a DataFrame
int_pump_df = pd.DataFrame(int_pump_data, columns=['file', 'cell', 'total_pumps', 'mean_int', 'norm_int_sub', 'norm_int_div', 'pct_brighter'])


# print(norm_line_df)
# print(int_pump_df)
print(f'''
Values normalized
''')
	
## group the pump and intensity dataframe by file
indv_files = int_pump_df.groupby("file")

rep_data = []

## calculate replicate averages 
for name, group in indv_files:
	
	## calculate the average number of pumps per minute per file
	avg_pumps = group["total_pumps"].mean()
	
	## calculate the average intensity of the cells per file, both raw and norm
	avg_intR = group["mean_int"].mean()
	avg_intNS = group["norm_int_sub"].mean()
	avg_intND = group["norm_int_div"].mean()
	avg_pctBright = group["pct_brighter"].mean()
	
	
	rep_data.append((name, avg_pumps, avg_intR, avg_intNS, avg_intND, avg_pctBright )) 

## convert to a dataframe
rep_avgs = pd.DataFrame(rep_data, columns=['file', 'avg_pumps', 'avg_intR', 'avg_intNS', 'avg_intND', 'avg_pctBright'])


# print(rep_avgs)

print(f'''
Replicate averages calculated
''')
	
## save the blinded dataframes as .csvs to the desired directory, dated for the current date
norm_line_df.to_csv(f'''{blinded_output_directory}/{current_date}_{lineScan_outfile}_B.csv''', index=False)
int_pump_df.to_csv(f'''{blinded_output_directory}/{current_date}_{pump_int_outfile}_B.csv''', index=False)
rep_avgs.to_csv(f'''{blinded_output_directory}/{current_date}_{avg_outfile}_B.csv''', index=False)


## Print out an aesthetically pleasing completion message for blinded files
## Calculate the length of the separator based on the output_directory length
separator1 = "~" * (len(blinded_output_directory)+ 20)  # Add padding for aesthetics

print(f'''
{separator1}
Blinded files written to:
{blinded_output_directory}
{separator1}
''')


## assign each file a timepoint and experiment date by merging with the blinding key 
norm_line_df = norm_line_df.merge(mapping_df, on='file', how='left')
int_pump_df = int_pump_df.merge(mapping_df, on='file', how='left')
rep_avgs = rep_avgs.merge(mapping_df, on='file', how='left')

# print(norm_line_df)
# print(int_pump_df)
# print(rep_avgs)


## group by experiment date then by cell
cell_df = int_pump_df.groupby(["date", "cell"])

## create a new tuple to add pump and intensity data, once for each cell per image
timing_dataT = []

## for each cell per replicate, 
## find the first time the cell is above the desired intensity threshold
## and find the first time the cell is below the desired pumping threshold
for name, group in cell_df:
	
	date = name[0]
	cell = name[1]
	
	## for each cell, keep only the rows where the cell is brighter than the threshold
	filtered = group[group['pct_brighter'] > threshold_pct]
# 	print(filtered)

	## determine the first time point where the cell is brighter than the threshold
	first_time_brighter = filtered.iloc[0]['time'] if not filtered.empty else None
	
	## for each cell, keep only the rows where the total pumps are equal to or below the threshold
	pumps_below_threshold = group[group['total_pumps'] <= theshold_pumps]
# 	print(f'''\n{no_pumps}\n\n''')

	## determine the first time point where the cell in equal to or below the pump threshold
	first_time_below = pumps_below_threshold.iloc[0]['time'] if not pumps_below_threshold.empty else None
	
	## if the cell has time point where it is below the pump threshold,
	## calculate the difference between this time and 
	## the first time the cell was above the intensity threshold
	if first_time_below:
		time_until_below = int(first_time_below) - int(first_time_brighter)
	
	## if a cell does not have a time point where it is below the pump threshold,
	## indicate the time it takes to get below the pump threshold
	## from the first time the cell was above the intensity threshold
	## is above the max time imaging took place
	else:
		time_until_below = f'''>{max_time}'''
	

	## take the variables and store them in a tuple
	timing_dataT.append((date, cell, first_time_brighter, first_time_below, time_until_below))


	
## Convert the collected pump and intensity data into a DataFrame
timing_data = pd.DataFrame(timing_dataT, columns=['date', 'cell', 'time_brighter', 'first_time_below', timing_column_name])

# print(timing_data)


# Count occurrences of each unique difference value
counts_df = timing_data[timing_column_name].value_counts().reset_index()

# Rename columns
counts_df.columns = [timing_column_name, 'count']

# Sort values if needed (optional)
counts_df = counts_df.sort_values(by=timing_column_name, key=lambda x: pd.to_numeric(x, errors='coerce').fillna(float('inf')))

# Convert the column to numeric, replacing ">40" with a reasonable number (e.g., 41)
# timing_data[timing_column_name] = timing_data[timing_column_name].replace({f'''>{max_time}''': max_time+2}).astype(float)  # Use time+1 for binning

# print(timing_data)



## save the unblinded dataframes as .csvs to the desired directory, dated for the current date
norm_line_df.to_csv(f'''{unblinded_output_directory}/{current_date}_{lineScan_outfile}.csv''', index=False)
int_pump_df.to_csv(f'''{unblinded_output_directory}/{current_date}_{pump_int_outfile}.csv''', index=False)
rep_avgs.to_csv(f'''{unblinded_output_directory}/{current_date}_{avg_outfile}.csv''', index=False)
# timing_data.to_csv(f'''{unblinded_output_directory}/{current_date}_{pumpStopTiming_outfile}.csv''', index=False)
counts_df.to_csv(f'''{unblinded_output_directory}/{current_date}_{counts_outfile}.csv''', index=False)


## Print out an aesthetically pleasing completion message for unblinded
## Calculate the length of the separator based on the output_directory length
separator2 = "~" * (len(unblinded_output_directory)+ 20)  # Add padding for aesthetics

print(f'''
{separator2}
Unblinded files written to:
{unblinded_output_directory}
{separator2}
''')



