### Negative Binomial Mixture Model for Identifica-tion of Noise in Antibody-Antigen Specificity Predictions from Single-Cell Data
Wasdin et al., Bioinformatics Advances, 2024.

-------------- 

Code for manuscript demonstrating use of model fitting and predictions on LIBRA-seq dataset with SARS-CoV-2 binding B cells. To fit to new data, adapt code in fig5_model_fitting.ipynb.

Notebooks included:

* fig1_dist_plots.ipynb: distributions of raw data shown in Figure 2
* preprocessing_pipeline.ipynb: basic filtering of raw data and separation of VRC01 cells
* fig4_model_comparison.ipynb: fitting of different types of mixture models for Figure 4
* fig5_model_fitting.ipynb: fitting model to Donor 2 UMI counts, with plots for Figure 5

These notebooks use data from /data directory and import functions from /MM_functions

Python libraries needed include:
* Numpy, Pandas, Seaborn, Matplotlib
* Scipy
* Levenshtein
* Biopython
* Statsmodels
