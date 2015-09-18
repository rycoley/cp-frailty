#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/dissert/cp-sims

module load R/3.1

#$ -cwd

#set environment variables

#output files
#$ -o cpe-cpt1.out
#$ -e cpe-cpt1.err

#request memory

#send seeds
#$ -t 1-250
#$ -V

# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH --vanilla /home/bst/student/rcoley/dissert/cp-sims/cpe-cpt1.R cpe-cpt1.Rout

#clean up temp files?

