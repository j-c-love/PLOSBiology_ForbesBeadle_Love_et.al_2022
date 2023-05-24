#For the analysis of the colocalisation of whole mRNAs, and lone mRNA ends,
#with P-bodies. This script requires staining of mRNAs with two sets of probes,
#one targeted to each end of the mRNA, tagged with different fluorophores, as
#well as P-bodies marked by a different fluorescent marker.

#mRNAs are paired to identify whole mRNAs and lone ends, and the
#colocalisation of each of these species with P-bodies is assessed.

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
pbodies_data = ''

#location to output data (string)
output_dir = ''

#read in .loc3 files from Airlocalize for single mRNAs, TSs and pbodies
start_singles = pd.read_csv(start_singles_data, sep="  ", header=None, engine='python')
start_TSs = pd.read_csv(start_TSs_data, sep="  ", header=None, engine='python')
end_singles = pd.read_csv(end_singles_data, sep="  ", header=None, engine='python')
end_TSs = pd.read_csv(end_TSs_data, sep="  ", header=None, engine='python')
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

#function to filter the end-P body pairs for ends that are within the radius of a P body (more stringent filter as one end)
def get_coloc_pbodies(singles):
    coloc_singles = singles[singles.distance_to_pbody <= 0.15]
    non_coloc_singles = singles[singles.distance_to_pbody > 0.15]
    return(coloc_singles, non_coloc_singles)

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

pbodies.columns = ['x', 'y', 'z', 'intensity', 'distance']
pbodies = pbodies[pbodies.intensity != -8]
pbodies['x'] = pbodies['x']*xy_voxels
pbodies['y'] = pbodies['y']*xy_voxels
pbodies['z'] = pbodies['z']*z_voxels
pbodies['index'] = pbodies.index

#give each object (here P-bodies therefore p1 p2 etc.) a unique identifier
pbodies['object_ID'] = ['p' + str(i) for i in range(0, pbodies.shape[0])]
pbodies.reset_index(inplace=True)

start_singles = remove_TSs(start_singles, start_TSs)
start_singles.reset_index(inplace=True)
start_singles['object_ID'] = ['s' + str(i) for i in range(0, start_singles.shape[0])]

end_singles = remove_TSs(end_singles, end_TSs)
end_singles.reset_index(inplace=True)
end_singles['object_ID'] = ['e' + str(i) for i in range(0, end_singles.shape[0])]

start_singles_pos = get_pos(start_singles)
end_singles_pos = get_pos(end_singles)
pbodies_pos = get_pos(pbodies)

# GET COLOCALISED MRNAS
rna_distances = get_distances_rnas(start_singles_pos, end_singles_pos, start_singles, end_singles)
coloc_singles, non_coloc_singles = get_colocs_rnas(rna_distances)
start_coloc = coloc_singles.start_object.tolist()
end_coloc = coloc_singles.end_object.tolist()

# GET LONE STARTS AND ENDS
lone_start_singles = start_singles[~start_singles['object_ID'].isin(start_coloc)]
lone_end_singles = end_singles[~end_singles['object_ID'].isin(end_coloc)]

# GET LOCALISATION WITH P BODIES
start_singles = get_distances_pbodies(start_singles_pos, pbodies_pos, start_singles, pbodies)
end_singles = get_distances_pbodies(end_singles_pos, pbodies_pos, end_singles, pbodies)
coloc_start_singles_pbodies, non_coloc_start_singles_pbodies= get_coloc_pbodies(start_singles)
coloc_end_singles_pbodies, non_coloc_end_singles_pbodies= get_coloc_pbodies(end_singles)

#FILTER LONE STARTS AND ENDS FOR THOSE IN P BODIES
start_in_pbodies_list = coloc_start_singles_pbodies.object_ID.tolist()
end_in_pbodies_list = coloc_end_singles_pbodies.object_ID.tolist()

lone_start_in_pbody = lone_start_singles[lone_start_singles['object_ID'].isin(start_in_pbodies_list)]
lone_end_in_pbody = lone_end_singles[lone_end_singles['object_ID'].isin(end_in_pbodies_list)]

#GET WHOLE RNAs WITH EITHER END IN P-BODY
start_in_pbodies_list = coloc_start_singles_pbodies.object_ID.tolist()
end_in_pbodies_list = coloc_end_singles_pbodies.object_ID.tolist()
coloc_singles_start_in_pbody = coloc_singles[coloc_singles['start_object'].isin(start_in_pbodies_list)]
coloc_singles_end_in_pbody = coloc_singles[coloc_singles['end_object'].isin(end_in_pbodies_list)]
coloc_pbody_either_end = coloc_singles_start_in_pbody.append(coloc_singles_end_in_pbody)
coloc_pbody_either_end.drop_duplicates(inplace=True)

lone_ends_summary_stats = {'probe_set': ['start', 'end'], 'total_singles' :[len(start_singles), len(end_singles)],'num_lone_singles': [len(lone_start_singles), len(lone_end_singles)],'num_lone_singles_in_pbodies': [len(lone_start_in_pbody), len(lone_end_in_pbody)],'proportion_in_pbodies': [(len(lone_start_in_pbody)/len(lone_start_singles)), (len(lone_end_in_pbody)/len(lone_end_singles))]}
lone_ends_summary_stats_df = pd.DataFrame(data=lone_ends_summary_stats)
whole_summary_stats = {'total_whole' :[len(coloc_singles)],'whole_in_Pbodies': [len(coloc_pbody_either_end)],'proportion_in_pbodies': [(len(coloc_pbody_either_end)/len(coloc_singles))]}
whole_summary_stats_df = pd.DataFrame(data=whole_summary_stats)

print(lone_ends_summary_stats_df)
print(whole_summary_stats_df)

#OUTPUT RESULTS
lone_ends_summary_stats_df.to_csv(output_dir + '/lone_ends_in_Pbodies_summary_data.csv')
whole_summary_stats_df.to_csv(output_dir + '/whole_mRNAs_in_Pbodies_summary_data.csv')
