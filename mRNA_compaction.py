#For the analysis of colocalisation of mRNA ends and the distance between the
#ends of an mRNA. Experimental data for this script should use two sets of probes,
#one targeted to each end of the mRNA, tagged with different fluorophores.

#The mRNA ends are paired, filtered using a distance threshold and then the
#distances between the ends of the paired mRNAs are determined.

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
start_singles_data = ''
start_TSs_data = ''
end_singles_data = ''
end_TSs_data = ''

#location to output data (string)
output_dir = ''

#read in .loc3 files from Airlocalize for single mRNAs and TSs
start_singles = pd.read_csv(start_singles_data, sep="  ", header=None, engine='python')
start_TSs = pd.read_csv(start_TSs_data, sep="  ", header=None, engine='python')
end_singles = pd.read_csv(end_singles_data, sep="  ", header=None, engine='python')
end_TSs = pd.read_csv(end_TSs_data, sep="  ", header=None, engine='python')

#DEFINE FUNCTIONS TO BE USED IN THE SCRIPT

#function to remove TSs from single mRNAs data
def remove_TSs(singles, transcription_sites):
    just_singles = singles[~singles.isin(transcription_sites)].dropna()
    return(just_singles)

#function to drop all non relevant columns
def get_pos(singles):
    singles_pos = singles.drop(columns=['intensity', 'distance', 'index', 'object_ID', 'level_0'])
    return(singles_pos)

#function to find the distances between pairs of start and end signals
def get_distances_rnas(start_singles_pos, end_singles_pos, start_singles, end_singles):
    distances = distance_matrix(start_singles_pos.values, end_singles_pos.values)

    #find pairs of colocalised spots minimising distance between all spots
    row_ind, col_ind = scipy.optimize.linear_sum_assignment(distances)

    #make a dataframe of the distances between pairs
    pair_distances = distances[row_ind, col_ind]
    end_singles_ids = end_singles.iloc[col_ind]
    end_singles_ids = end_singles_ids.drop(columns=['level_0', 'x', 'y', 'z', 'intensity', 'distance', 'index'])
    end_singles_ids.columns = ['end_object_ID']
    start_singles_ids = start_singles.iloc[row_ind]
    start_singles_ids = start_singles_ids.drop(columns=['level_0', 'x', 'y', 'z', 'intensity', 'distance', 'index'])
    start_singles_ids.columns = ['start_object_ID']
    distances_df = pd.DataFrame(data=pair_distances, columns=['distance'])
    final_df = pd.DataFrame(np.vstack((start_singles_ids['start_object_ID'],end_singles_ids['end_object_ID'], pair_distances)).T, columns=['start_object', 'end_object', 'rna_distance'])
    final_df['rna_ID'] = ['r' + str(i) for i in range(0, final_df.shape[0])]
    return(final_df)

#function to filter the paired ends for those that are close enough to be considered ends of the same RNA
def get_colocs_rnas(rna_distances):
    coloc_df = rna_distances[rna_distances.rna_distance <= 0.3]
    non_coloc_df = rna_distances[rna_distances.rna_distance > 0.3]
    return(coloc_df, non_coloc_df)

#MAIN BODY OF CODE

end_singles.columns = ['x', 'y', 'z', 'intensity', 'distance']
end_singles = end_singles[end_singles.intensity != -8]
end_singles['x'] = end_singles['x']*xy_voxels
end_singles['y'] = end_singles['y']*xy_voxels
end_singles['z'] = end_singles['z']*z_voxels
end_singles['index'] = end_singles.index

end_TSs.columns = ['x', 'y', 'z', 'intensity', 'distance']
end_TSs = end_TSs[end_TSs.intensity != -8]
end_TSs['x'] = end_TSs['x']*xy_voxels
end_TSs['y'] = end_TSs['y']*xy_voxels
end_TSs['z'] = end_TSs['z']*z_voxels

start_singles.columns = ['x', 'y', 'z', 'intensity', 'distance']
start_singles = start_singles[start_singles.intensity != -8]
start_singles['x'] = start_singles['x']*xy_voxels
start_singles['y'] = start_singles['y']*xy_voxels
start_singles['z'] = start_singles['z']*z_voxels
start_singles['index'] = start_singles.index

start_TSs.columns = ['x', 'y', 'z', 'intensity', 'distance']
start_TSs = start_TSs[start_TSs.intensity != -8]
start_TSs['x'] = start_TSs['x']*xy_voxels
start_TSs['y'] = start_TSs['y']*xy_voxels
start_TSs['z'] = start_TSs['z']*z_voxels

start_singles = remove_TSs(start_singles, start_TSs)
start_singles.reset_index(inplace=True)
#give each object a unique identifier
start_singles['object_ID'] = ['s' + str(i) for i in range(0, start_singles.shape[0])]

end_singles = remove_TSs(end_singles, end_TSs)
end_singles.reset_index(inplace=True)
end_singles['object_ID'] = ['e' + str(i) for i in range(0, end_singles.shape[0])]

start_singles_pos = get_pos(start_singles)
end_singles_pos = get_pos(end_singles)

#get distances between mRNA ends
rna_distances = get_distances_rnas(start_singles_pos, end_singles_pos, start_singles, end_singles)

#find colocalised mRNA ends
coloc_singles, non_coloc_singles = get_colocs_rnas(rna_distances)
print(coloc_singles)

#create a dataframe of the summary statistics (mean and median) for compaction of mRNAs inside and outside P bodies
summary_stats = {'Mean_compaction': [coloc_singles.rna_distance.mean()],
'Median_compaction': [coloc_singles.rna_distance.median()]}
summary_stats_df = pd.DataFrame(data=summary_stats)
print(summary_stats_df)

#OUTPUT RESULTS
coloc_singles.to_csv(output_dir + '/coloc_singles.csv')
summary_stats_df.to_csv(output_dir + '/compaction_summary_stats.csv')
