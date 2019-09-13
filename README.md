# nutrients-plant-viruses

### Code and data associated with the manuscript "Nutrients mediate plant virus interactions through within-host density and transmission"

### Contents
- code: descriptions below
- output: outputs of R scripts that are used as inputs to other R scripts (also in code directory) or as figures in manuscript
- data: descriptions below (publicly available: https://doi.org/10.6073/pasta/01e7bf593676a942f262623710acba13)
- metadata: EML-formatted metadata for each data file, see methods.txt for details on data collection

|code                                     |desription                                                                                        |
|:----------------------------------------|:-------------------------------------------------------------------------------------------------|
|coinfection_correlation_figure.R         |code to create figure of correlation between PAV and RPV in coinfection                           |
|concentration_analysis.R                 |code to analyze density of viruses within source plants                                           |
|concentration_figure_results.R           |code to create figure of treatment effects on virus density                                       |
|concentration_transmission_figure.R      |code to to create figure of the relationship between virus density and transmission               |
|infection_analysis.R                     |code to analyze infection status of source plants                                                 |
|lacroix_concentration_priors.R           |code to derive priors for virus density analysis                                                  |
|lacroix_transmission_priors.R            |code to derive priors for transmission analysis                                                   |
|qPCR_raw_data_processing.R               |code to process raw qPCR data files                                                               |
|supplementary_infection_figure_results.R |code to create supplementary figure of treatment effects on infection status                      |
|supplementary_prior_comparison_figure.R  |code to create supplementary figure of comparison between models with and without priors          |
|supplementary_truncated_dataset_figure.R |code to create supplementary figure of comparison between models with full and truncated datasets |
|transmission_analysis.R                  |code to analyze transmission                                                                      |
|transmission_figure_results.R            |code to create figure of treatment effects on transmission                                        |

|data                              |description                                                                                                                                                     |
|:---------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------|
|qPCR_group_comments.csv           |comments about controls and standards for each qPCR group                                                                                                       |
|sample_exp_molc_data.csv          |sample experiment molecular data - includes experimental treatments, measurements taken during harvesting, and molecular analysis information for source plants |
|source_plant_chlorophyll_data.csv |chlorophyll measurements for source plants                                                                                                                      |
|source_plant_qPCR_data.csv        |qPCR data for source plants                                                                                                                                     |
|transmission_data.csv             |experimental treatments and molecular analysis information for receiving plants in transmission trials                                                          |

## Citation
Kendig A. E., E. T. Borer, E. N. Boak, T. C. Picard, E. W. Seabloom. 2019. Soil nitrogen and phosphorus effects on plant virus density, transmission, and species interactions. Environmental Data Initiative. https://doi.org/10.6073/pasta/01e7bf593676a942f262623710acba13. Dataset accessed 9/13/2019.

## License
see https://doi.org/10.6073/pasta/01e7bf593676a942f262623710acba13 