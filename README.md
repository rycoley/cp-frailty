cp-frailty
==========

compound poisson frailty model

Simulation and application code for Coley and Brown (2015) “Estimating
effectiveness in HIV prevention trials with a hierarchical compound
Poisson frailty model”

R code for application to HPTN 035 data contained in "hptn035-strat.R" and "hptn035-strat-source.R". Unfortunately, the data for this application is not publicly available. The application code here is an extension of the simulation code in that it also accommodates (1) interval-censored event times and (2) a more complicated baseline hazard specification and prior. Details for the baseline hazard and prior are given in Coley and Brown (2015). Code generally follows parameter and covariates notation used in paper. 



Simulation code can be found in the simulations folder. Here is a general explanation of what is contained in this folder.

cpe-cpt.R, cpe-to-run.R, cpe-source.R, cpe-cpt.sh
These scripts reproduce the simulation with results shown in Figure 1 and Table 1. As detailed in Coley and Brown (2015), this simulation varies the proportion at risk while maintaining the same values for other model parameters. “Correct, broadly informative” priors were used, that is the expected value of the prior equaled the true data-generating value for each parameter. While informative, priors were sufficiently board to allow for a wide range of values in posterior samples.

cpe-cpt1.R, cpe-to-run.R, cpe-source.R, cpe-cpt1.sh
These scripts reproduce the prior mispecifcation simulation described in (I) of Section 3.1. An incorrect prior is placed on P(Risk). Results are given in Table 2.


cpe-cpt2.R, cpe-to-run.R, cpe-source.R, cpe-cpt2.sh
These scripts reproduce the prior mispecifcation simulation described in (II) of Section 3.1. An incorrect prior is placed on the baseline hazard. Results are given in Table 2.

cpe-cpt3.R, cpe-to-run3.R, cpe-source3.R, cpe-cpt3.sh
These scripts reproduce the prior mispecifcation simulation described in (III) of Section 3.1. Non-informative priors are placed on all model parameters. Results are given in Table 2.


Batch scripts (.sh) were used to run simulations on a Sun Grid Engine (SGE) cluster. Event rate (ER) and proportion not at risk (PNR) were exported to all jobs in the qsub command using -v, for e.g. qsub -v ER=0.06 PNR=0.75 cpe-cpt.sh. Starting seeds for each simulation were sent to jobs via task id, see “-t 1-250 -V” in batch scripts.

Batch scripts will need to be edited by the user to contain local working directories.