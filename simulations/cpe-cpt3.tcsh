#!/bin/tcsh
#

#set working directory.
cd /home/bst/student/rcoley/dissert

#$ -cwd

#export environment variables to job in qsub function using -v command 

#output and error files #name these in qsub function to reflect environment variables, e.g. PNR and ER


#request nodes
#$ -pe local 20


# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH --vanilla /home/bst/student/rcoley/dissert/cpe-cpt3.R cpe-cpt3-75-06.Rout


