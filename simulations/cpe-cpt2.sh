#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/dissert/cp-sims

module load R/3.1

#$ -cwd

#set environment variables

#output files
#$ -o cpe-cpt2.out
#$ -e cpe-cpt2.err

#request memory

#send seeds
#$ -t 1-250
#$ -V

# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH --vanilla /home/bst/student/rcoley/dissert/cp-sims/cpe-cpt2.R cpe-cpt2.Rout

#clean up temp files?

