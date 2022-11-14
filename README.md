# PLOSBiology_ForbesBeadle_Love_2022
Code to accompany "Single molecule imaging and modelling of mRNA decay dynamics in the Drosophila embryo", PLOS Biology, 2022

These scripts are for the analysis of Single molecule fluorescence in situ hybridisation (smFISH) imaging data.

Scripts included are as follows:
- mRNA_compaction.py : For the analysis of colocalisation of mRNA ends and the distance between the ends of an mRNA. Experimental data for this script should use two sets of probes, one targeted to each end of the mRNA, tagged with different fluorophores. The mRNA ends are paired, filtered using a distance threshold and then the distances between the ends of the paired mRNAs are determined.

- Pbody_colocalisation.py : For the analysis of the colocalisation of mRNAs with cytoplasmic Processing bodies (P-bodies). Data for this script should contain mRNAs tagged with probes of one colour, and P-bodies marked by a different fluorescent marker. mRNAs are assigned to their nearest P-body, filtered based on the radius of a P-body and then the P-body colocalisation index is calculated (see paper methods).

- lone_ends_whole_mRNAs_in_Pbodies.py : For the analysis of the colocalisation of whole mRNAs, and lone mRNA ends, with P-bodies. This script requires staining of mRNAs with two sets of probes, one targeted to each end of the mRNA, tagged with different fluorophores, as well as P-bodies marked by a different fluorescent marker. mRNAs are paired to identify whole mRNAs and lone ends, and the colocalisation of each of these species with P-bodies is assessed as above.

All scripts require user input of filepaths to input data, voxel dimensions of the specific images analysed and filepath for the output summary data.

For further details of methods, please see the above paper.
