#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/dissert/cp-sims

module load R/3.1

#$ -cwd

#set environment variables

#output files
#$ -o cpe-cpt.out
#$ -e cpe-cpt.err

#request memory
#$ -l h_vmem=10G

#send seeds
#$ -t 1-250
#$ -V

# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH /home/bst/student/rcoley/dissert/cp-sims/cpe-cpt.R cpe-cpt.Rout

#clean up temp files?

