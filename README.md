# Ecoacoustic_indicators

Bold scripts can be run using the Rproj workspace and this repository, while italicized scripts can be run using data found in the associated Zenodo data repository with small, local pathway setups. Scripts in non-stylized text require adaptation to an HPC or local environment.
0_acoustic_index_calculations
- 0_acoustic_index_calculations.R & 0_acoustic_index_calculations.sh: calculates the 15 acoustic indices here using an HPC environment
- 1_completion_check.R: gather full path names needed to run acoustic index calculations 

1_data_processing
- *0_cnn_by_site_predictions.R*: gathers by-min ABGQI predictions to site csvs
- *1a_abgqi_aggregation.R*: averages site-level ABGQI statistics
- *1b_acoustic_index_aggregation.R*: averages site-level acoustic index statistics

2_soundscape_modeling
- **0_by_site_GAM.R**: uses by-site acoustic index csvs and ABGQI csv for Eq.1 generalized additive model (GAM) selction
- **1_GAM_result_analysis.R**: uses GAM models and slope summaries to analyze model results. Generates GAM (Eq.1) model summary txt files, Figure 2, and Table 3 content.
- **2_ACI_interference_correction.R**: uses GAM models and slope summaries, Zenodo by-site ABGQI and acoustic indices to generate non-Interference summary to generate Figures 3 and 4. Can be run using only the GitHub repo with the precomputed no-Interference csvs provided.

3_richness_modeling
- **0_bird_species_richness.R**: combines acoustic index, ABGQI, bird species richness to analyze Eqs 2,3,4,5
- **1_bird_species_richness_morning.R**: combines acoustic index, ABGQI, bird species richness to analyze Eqs 2,3,4,5 for targeted 4 a.m. to 12 p.m.

Figures and tables in the primary manuscript generated from scripts in this repository:
Figures:
- Fig 2: 1_GAM_result_analysis.R
- Fig 3: 2_ACI_interference_correction.R
- Fig 4: 2_ACI_interference_correction.R
- Fig 5: 1_bird_species_richness_morning.R

Tables:
- Table 1: 1_GAM_result_analysis.R
- Table 2: 1_GAM_result_analysis.R
- Table 3: 1_GAM_result_analysis.R
- Table 4: 0_bird_species_richness.R and 1_bird_species_richness_morning.R