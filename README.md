# Dementia_progression

This repository contains all the code used to replicate, review, and extend the work done by Williams et al. (https://doi.org/10.1080/01621459.2019.1594831) to study the effect of ageing on dementia. It being part of my work to credit the Research Comprehensive Exam of the PhD program from the Department of Statistical Sciences at the University of Toronto.
The code was used to conduct simulation studies (see the report above) regarding the continuous-time hidden Markov model (HMM) developed by Williams et al..

Te repository is organized through three folders:

- "New_work" contains all of our original work made for this project. Some of the scripts contain code adapted from the work by Williams et al.. When that is the case, the part of the script that was adapted is indicated. 

- "Supplemental" contains results like those shown in Appendix 2 of the report for every parameter in the dementia progression model. 

- "TimeCapsule_Code" contains all the work done by Williams et al., as made publicly available in https://jonathanpw.github.io/research.html . All the material contained there is either unchanged, or minimally adapted to be used for this study. 

The scripts were used for the study as follows: 

- "Simulate_cav.r" (New_work) was used to generate all synthetic CAV data-sets.

- "fit_cav_models.r (New_work) was used to fit all CAV models whose results are shown in Figures 3 and 4 from the report.

- "delayed_enrollment_plots.r (New_work) was used to produce plots as the ones shown in Figures 3 and 4 from the report.

- "Simulate.r" (TimeCapsule_Code) was used to create the MCSA synthetic data-set used for Section 4.4 in the report.

- "RunFile.r" (imeCapsule_Code) was used to implement the dementia progression model.

- "OutFile_figures.r (TimeCapsule_Code) was used to produce the Figures in Section 4.4, and the material contained in "Supplemental".

- "cav.stan" and "mcsa.stan" (New_work) contain the CAV and MCSA models respectively for implementation using Stan.

- "fit_ext_models.r (New_work) was used to produce the results from Section 5.2 in the report.
