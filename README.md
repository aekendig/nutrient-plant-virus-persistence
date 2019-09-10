# nutrients-plant-viruses

### Code and data associated with the manuscript "Nutrients mediate plant virus interactions through within-host density and transmission"

### Contents
- code: descriptions below
- output: outputs of R scripts that are used as inputs to other R scripts (also in code directory) or as figures in manuscript
- data: descriptions below (will be made publicly available on https://portal.edirepository.org)
- metadata: EML-formatted metadata for each data file, see methods.txt for details on data collection

|code                                     |desription                                                                                        |
|:----------------------------------------|:-------------------------------------------------------------------------------------------------|
|coinfection_correlation_figure.R         |code to create figure of correlation between PAV and RPV in coinfection                           |
|concentration_analysis.R                 |code to analyze density of viruses within source plants                                           |
|concentration_figure_results.R           |code to create figure of treatment effects on virus density                                       |
|concentration_transmission_figure.R      |code to to create figure of the relationship between virus density and transmission               |
|eml_metadata.R                           |code to create metadata                                                                           |
|infection_analysis.R                     |code to analyze infection status of source plants                                                 |
|lacroix_concentration_priors.R           |code to create figure of treatment effects on infection status                                    |
|lacroix_transmission_priors.R            |code to derive priors for virus density analysis                                                  |
|qPCR_raw_data_processing.R               |code to derive priors for transmission analysis                                                   |
|supplementary_infection_figure_results.R |code to process raw qPCR data files                                                               |
|supplementary_prior_comparison_figure.R  |code to create supplementary figure of comparison between models with and without priors          |
|supplementary_truncated_dataset_figure.R |code to create supplementary figure of comparison between models with full and truncated datasets |
|transmission_analysis.R                  |code to analyze transmission                                                                      |
|transmission_figure_results.R            |code to create figure of treatment effects on transmission                                        |

|data                         |description                                                                                                                                                     |
|:----------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------|
|group_01_qPCR_data.csv       |group 01 qPCR data for source plants                                                                                                                            |
|group_02_qPCR_data.csv       |group 02 qPCR data for source plants                                                                                                                            |
|group_03_qPCR_data.csv       |group 03 qPCR data for source plants                                                                                                                            |
|group_04_qPCR_data.csv       |group 04 qPCR data for source plants                                                                                                                            |
|group_05_qPCR_data.csv       |group 05 qPCR data for source plants                                                                                                                            |
|group_06_qPCR_data.csv       |group 06 qPCR data for source plants                                                                                                                            |
|group_07_qPCR_data.csv       |group 07 qPCR data for source plants                                                                                                                            |
|group_08_qPCR_data.csv       |group 08 qPCR data for source plants                                                                                                                            |
|group_09_qPCR_data.csv       |group 09 qPCR data for source plants                                                                                                                            |
|group_10_qPCR_data.csv       |group 10 qPCR data for source plants                                                                                                                            |
|group_11_qPCR_data.csv       |group 11 qPCR data for source plants                                                                                                                            |
|group_12_qPCR_data.csv       |group 12 qPCR data for source plants                                                                                                                            |
|group_13_qPCR_data.csv       |group 13 qPCR data for source plants                                                                                                                            |
|group_14_qPCR_data.csv       |group 14 qPCR data for source plants                                                                                                                            |
|group_15_qPCR_data.csv       |group 15 qPCR data for source plants                                                                                                                            |
|group_16_qPCR_data.csv       |group 16 qPCR data for source plants                                                                                                                            |
|group_17_qPCR_data.csv       |group 17 qPCR data for source plants                                                                                                                            |
|group_18_qPCR_data.csv       |group 18 qPCR data for source plants                                                                                                                            |
|group_19_qPCR_data.csv       |group 19 qPCR data for source plants                                                                                                                            |
|group_20_qPCR_data.csv       |group 20 qPCR data for source plants                                                                                                                            |
|group_21_qPCR_data.csv       |group 21 qPCR data for source plants                                                                                                                            |
|group_22_qPCR_data.csv       |group 22 qPCR data for source plants                                                                                                                            |
|group_23_qPCR_data.csv       |group 23 qPCR data for source plants                                                                                                                            |
|group_24_qPCR_data.csv       |group 24 qPCR data for source plants                                                                                                                            |
|group_25_qPCR_data.csv       |group 25 qPCR data for source plants                                                                                                                            |
|group_26_qPCR_data.csv       |group 26 qPCR data for source plants                                                                                                                            |
|group_27_qPCR_data.csv       |group 27 qPCR data for source plants                                                                                                                            |
|group_28_qPCR_data.csv       |group 28 qPCR data for source plants                                                                                                                            |
|group_29_qPCR_data.csv       |group 29 qPCR data for source plants                                                                                                                            |
|group_30_qPCR_data.csv       |group 30 qPCR data for source plants                                                                                                                            |
|group_31_qPCR_data.csv       |group 31 qPCR data for source plants                                                                                                                            |
|group_32_qPCR_data.csv       |group 32 qPCR data for source plants                                                                                                                            |
|group_33_qPCR_data.csv       |group 33 qPCR data for source plants                                                                                                                            |
|group_34_qPCR_data.csv       |group 34 qPCR data for source plants                                                                                                                            |
|group_35_qPCR_data.csv       |group 35 qPCR data for source plants                                                                                                                            |
|group_36_qPCR_data.csv       |group 36 qPCR data for source plants                                                                                                                            |
|group_38_qPCR_data.csv       |group 38 qPCR data for source plants                                                                                                                            |
|group_39_qPCR_data.csv       |group 39 qPCR data for source plants                                                                                                                            |
|group_40_qPCR_data.csv       |group 40 qPCR data for source plants                                                                                                                            |
|group_41_qPCR_data.csv       |group 41 qPCR data for source plants                                                                                                                            |
|group_42_qPCR_data.csv       |group 42 qPCR data for source plants                                                                                                                            |
|group_43_qPCR_data.csv       |group 43 qPCR data for source plants                                                                                                                            |
|group_44_qPCR_data.csv       |group 44 qPCR data for source plants                                                                                                                            |
|group_45_qPCR_data.csv       |group 45 qPCR data for source plants                                                                                                                            |
|qPCR_group_comments.csv      |comments about controls and standards for each qPCR group                                                                                                       |
|sample_exp_molc_data.csv     |sample experiment molecular data - includes experimental treatments, measurements taken during harvesting, and molecular analysis information for source plants |
|source_plant_chlorophyll.csv |chlorophyll measurements for source plants                                                                                                                      |
|transmission_data.csv        |experimental treatments and molecular analysis information for receiving plants in transmission trials                                                          |

