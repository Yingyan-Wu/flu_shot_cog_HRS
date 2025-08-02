#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o /u/home/y/yingyanw/flu_cog_HRS/joblogs/joblog.$JOB_ID.$TASK_ID #creates a file called joblog.jobidnumber to write to. 
#$ -j y 
#$ -l h_rt=24:00:00,h_data=12G #requests 24 hours, 12GB of data (per core)
#$ -pe shared 5 #requests 5 cores
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system 
# Notify when
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
# submit array job:
# run an array of n jobs; this corresponds to $SGE_TASK_ID, nrow of scenario_tib
#$ -t 1-10:1
## 

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.2.2 #loads R/4.2.2 for use 
export OMP_NUM_THREADS=5 #uses max 1 threads (needs to match -pe shared)


echo "======"
echo SGE_TASK_ID=$SGE_TASK_ID      
R CMD BATCH --no-save --no-restore "--args scenario_num=$SGE_TASK_ID"  /u/home/y/yingyanw/flu_cog_HRS/scripts/5.1_sens_positive_control_cluster_TMLE_superlearner.R /u/home/y/yingyanw/flu_cog_HRS/console_output/output.$JOB_ID.$SGE_TASK_ID
echo R CMD BATCH --no-save --no-restore "--args scenario_num=$SGE_TASK_ID"  /u/home/y/yingyanw/flu_cog_HRS/scripts/5.1_sens_positive_control_cluster_TMLE_superlearner.R /u/home/y/yingyanw/flu_cog_HRS/console_output/output.$JOB_ID.$SGE_TASK_ID