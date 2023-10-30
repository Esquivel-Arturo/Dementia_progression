This directory contains the code to reproduce the results in the paper, A Bayesian 
Approach to Multi-State Hidden Markov Models: Application to Dementia Progression.  
Please contact Jonathan Williams at williams.jonathan1@mayo.edu for any help or 
questions.

See workflow.sh for line-by-line Unix command line code for reproducing the results.

synthetic_MCSA_data.rda has been simulated from a modified version of the script 
file Simulate.r to closely resemble the real MCSA data set.  The real MCSA data set 
is not provided due to HIPAA restrictions.

Transitions_heat_movie.mp4 is a movie to visualize the evolution of the transition 
rates in the state space over all integer ages from 50-90.





install.packages("~/u/esqfuente/TimeCapsule_Code/hmm_1.1-6.tar.gz", repos = NULL, type = "source")

/TimeCapsule_Code/hmm_1.1-6.tar.gz

MoonriseFoiledCrucialAlarm

Probar distintas priors para un conjunto de datos, ver si lleva a los outliers observados 
Tratar de implementar b-splines, ya sea para el efecto de los covariates en las transiciones o en las emisiones.

-Checar resultados de corrida actual, checar coverage
-try a different seed for dataset 1
-maybe storing the models instead of just the results?
- maybe contact them to ask about mile bias 


for seed in {1..5}
do
Rscript Simulate.r $seed Simulation_output/
Rscript RunFile.r $seed Simulation Bayesian no_scale 3 Simulation_output/
done




Compressed sensing 