#!/bin/bash

#$ -S /bin/bash #defines bash as the shell for execution
#$ -N admix #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -o ~/log #standard output file or directory (joined with error because of -j y)
#$ -q archgen.q #queue
#$ -pe smp 4 #needs X CPU cores
#$ -l h_vmem=10G #request XGb of memory
#$ -V # load personal profile
#$ -t 1-5 # array job length (here: amount of Ks * amount of reps)
#$ -tc 20 # number of concurrently running tasks in array

## Set ADMIXTURE run parameters. Output folder, input .bed file, minimum and maximum K values, and Number of CPUs per admixture run.
fn0="/mnt/archgen/users/schmid/paagen/playground/admixpops_test_data/admixture_test" ## The intended result directory. It will be created, with one subfolder per K (which itself contains one folder per replicate and one 'Logs' folder with the ADMIXTURE logfiles).
bedFile="/mnt/archgen/users/schmid/paagen/playground/admixpops_test_data/mbutihanfrench_merged/mbutihanfrench_merged.bed" ## The input .bed file you ran ADMIXTURE on.
Kmin='3' ## The minimum number of Ks you want to run
Kmax='3' ## The maximum number of Ks you want to run
Reps='5' ## The number of replicates to run for each K value. We normally use 5 replicates per K.
NumCPUs='4' ## The number of CPUs you wish to use. Make sure this number matches the '-c' option for SBATCH above.

date

## Create the bash arrays which we will iterate through with each array job.
AllKs=($(seq ${Kmin} ${Kmax}))
AllReps=($(seq 1 ${Reps}))
 
## Use the SGE_TASK_ID to iterate over the Ks and Reps in the correct manner (only iterate over AllKs once for every full iteration over AllReps)
i=$((SGE_TASK_ID - 1))
CurrentK=${AllKs[`expr ${i} / ${#AllReps[@]}`]}
CurrentRep=${AllReps[`expr ${i} % ${#AllReps[@]}`]}
 
## Make necessary output directories if needed.
mkdir -p ${fn0}/
cd ${fn0}
mkdir -p ${CurrentK}/Logs
mkdir -p ${CurrentK}/${CurrentRep}

## Finally, run ADMIXTURE.
cd ${CurrentK}/${CurrentRep}
admixture --supervised -j${NumCPUs} -s ${RANDOM} --cv ${bedFile} ${CurrentK} 1>${fn0}/${CurrentK}/Logs/K${CurrentK}_${CurrentRep}.log

date

exit 0
