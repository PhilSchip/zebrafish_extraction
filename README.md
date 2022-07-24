# Supplementary Material

The code as well as the raw data provided is part of the supplementary material of the paper "Evaluation of Extraction Methods for Untargeted Metabolomic Studies for Future Applications in a Zebrafish Tuberculosis Model" by Philip Schippers et al.

# Centwaveopt_modified

Contains a slightly modified centWaveOpt algorithm from the paper "Automated Optimization of XCMS Parameters for Improved Peak Picking of LC-MS Data using the Coefficient of Variation and Parameter Sweeping for Untargeted Metabolomic", Drug Testing and Analysis, 2018, by Sascha K. Manier et al.
The original version can be found here: https://github.com/saskema/centWaveOpt

# metaboLib_modified

Contains a slightly modified version of the metaboLib repository by Sascha K. Manier.
The original version can be found here: https://github.com/saskema/metaboLib

# Raw_Data

Contains the raw data in mzXML format, sorted by used column and polarity and further by sample groups.

# Inclusion_Lists

Contains the inclusion lists used for identifying molecules of interest (MoInt) or internal standards (IS).
For MoInt the inclusion lists are sorted by column and polarity, while for IS only a single identical inclusion list was used.

# Evaluation_Script

Contains the used R.Script for the evaluation of the raw data.

## How to Use

For each column/polarity combination:

Copy the Centwaveopt script and the raw data of the QC samples in the same folder and run the script for optimized parameters.

Copy the optimized parameters, the raw data of all samples, the metaboLib repository, the respective inclusion lists and the evaluation script into the same folder, change the parameters inside the evaluation script accordingly and run it.

The following datasets are generated:
  
### In the subfolder "Statistics_Results":

"Evaluation":             a PDF containing the results of the statistical evaluation (ANOVA, Heatmap, t-SNE, PCA and PC-DFA)

"Significant_features":   a CSV document, containing the data of each significant feature  
                          (name, pvalue, column, polarity, annotation of isototpes/adducts, pcgroup and absolute loading in the PC-DFA)
  
### In the subfolder "Target_Evaluation":
Extracted Ion chromatograms of each detected MoInt or IS, depending on the evaluation
  
### In the original folder:
"peaklist_CAMERA": a CSV document, containing the data of each feature  
(m/z and retention time values, peak count, peak areas and annotation of isotopes/adducts)
  
  
