#For the analysis of the colocalisation of mRNAs with cytoplasmic Processing
#bodies (P-bodies). Data for this script should contain mRNAs tagged with probes
#of one colour, and P-bodies marked by a different fluorescent marker.

#mRNAs are assigned to their nearest P-body, filtered based on the radius of a P-body
#and then the P-body colocalisation index is calculated (see paper methods).

import pandas as pd
import numpy as np
import scipy
from scipy.spatial import distance_matrix
from scipy.optimize import linear_sum_assignment

#REQUIRED USER INPUTS
#voxels of specific image (float)
xy_voxels =
z_voxels =

#input data file names and locations (.loc3 files from Airlocalize) (strings)
singles_data = ''
TSs_data = ''
pbodies_data = ''

#location to output data (string)
output_dir = ''

#read in .loc3 files for singles, TSs and pbodies
singles = pd.read_csv(singles_data, sep="  ", header=None, engine='python')
TSs = pd.read_csv(TSs_data, sep="  ", header=None, engine='python')
pbodies = pd.read_csv(pbodies_data, sep="  ", header=None, engine='python')

#DEFINE FUNCTIONS TO BE USED IN THE SCRIPT

#function to remove TSs from singles data
def remove_TSs(singles, transcription_sites):
    just_singles = singles[~singles.isin(transcription_sites)].dropna()
    return(just_singles)

#function to drop all non relevant columns
def get_pos(singles):
    singles_pos = singles.drop(columns=['intensity', 'distance', 'index', 'object_ID', 'level_0'])
    return(singles_pos)

#function to get distances between ends and P-bodies and assign end to nearest P-body
def get_distances_pbodies(singles, pbodies, singles_loc, pbodies_loc):
    distances = distance_matrix(pbodies.values, singles.values)
    min_distance = np.amin(distances, axis=0)
    min_distance = pd.Series(min_distance)
    min_index = np.argmin(distances, axis=0)
    min_index = pd.Series(min_index)
    pbodies_ids = pbodies_loc.iloc[min_index]
    pbodies_ids = pbodies_ids.drop(columns=['level_0', 'x', 'y', 'z', 'intensity', 'distance', 'index'])
    pbodies_ids.columns = ['pbody_ID']
    singles_results = pd.DataFrame(np.vstack((singles_loc['object_ID'],pbodies_ids['pbody_ID'], min_distance)).T, columns=['object_ID', 'pbody_object', 'distance_to_pbody'])
    return(singles_results)

#function to filter the end-P body pairs for ends that are within the radius of a P-body
def get_coloc_pbodies(singles):
    coloc_singles = singles[singles.distance_to_pbody <= 0.2]
    non_coloc_singles = singles[singles.distance_to_pbody > 0.2]
    return(coloc_singles, non_coloc_singles)

#function to output data as summary statistics
def get_stats(coloc_singles, non_coloc_singles, pbodies, results):
    total_singles = len(coloc_singles.index) + len(non_coloc_singles.index)
    total_pbodies = len(pbodies)
    coloc_singles_num = len(coloc_singles)
    CI_200 = (coloc_singles_num/(total_singles*total_pbodies))
    results = results.append({'Total_singles':total_singles, 'Total_pbodies':total_pbodies, 'Coloc_singles':coloc_singles_num, 'CI':CI_200}, ignore_index=True)
    return(results)

#MAIN BODY OF CODE
data = pd.DataFrame([], columns=['Total_singles', 'Total_pbodies', 'Coloc_singles', 'CI'], dtype=object)

singles.columns = ['x', 'y', 'z', 'intensity', 'distance']
singles = singles[singles.intensity != -8]
singles['x'] = singles['x']*xy_voxels
singles['y'] = singles['y']*xy_voxels
singles['z'] = singles['z']*z_voxels
singles['index'] = singles.index

TSs.columns = ['x', 'y', 'z', 'intensity', 'distance']
TSs = TSs[TSs.intensity != -8]
TSs['x'] = TSs['x']*xy_voxels
TSs['y'] = TSs['y']*xy_voxels
TSs['z'] = TSs['z']*z_voxels

pbodies.columns = ['x', 'y', 'z', 'intensity', 'distance']
pbodies = pbodies[pbodies.intensity != -8]
pbodies['x'] = pbodies['x']*xy_voxels
pbodies['y'] = pbodies['y']*xy_voxels
pbodies['z'] = pbodies['z']*z_voxels
pbodies['index'] = pbodies.index

#give each object (here P bodies therefore p1 p2 etc.) a unique identifier to make the results more readable
pbodies['object_ID'] = ['p' + str(i) for i in range(0, pbodies.shape[0])]
pbodies.reset_index(inplace=True)

singles = remove_TSs(singles, TSs)
singles.reset_index(inplace=True)
singles['object_ID'] = ['s' + str(i) for i in range(0, singles.shape[0])]
print(singles)
print(pbodies)

singles_pos = get_pos(singles)
pbodies_pos = get_pos(pbodies)

# get distances and assign mRNA to closest P-body
singles = get_distances_pbodies(singles_pos, pbodies_pos, singles, pbodies)

# filter for those within the radius of their nearest P-body (then colocalised)
coloc_singles_pbodies, non_coloc_singles_pbodies = get_coloc_pbodies(singles)

# create summary df of numbers
summary_df = get_stats(coloc_singles_pbodies, non_coloc_singles_pbodies, pbodies, data)
print(summary_df)

#OUTPUT RESULTS
summary_df.to_csv(output_dir + '/pbody_coloc_summary_data.csv')
