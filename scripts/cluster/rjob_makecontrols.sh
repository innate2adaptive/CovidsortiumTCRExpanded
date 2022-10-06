#!/bin/bash -l

# Example jobscript to run a single core R job

# Request ten minutes of wallclock time (format hours:minutes:seconds).
# Change this to suit your requirements.
#$ -l h_rt=6:0:0
# Request 1 gigabyte of RAM. Change this to suit your requirements.
#$ -l mem=32G

# Set the name of the job. You can change this if you wish.
#$ -N makeControls

# Set the working directory to somewhere in your scratch space.  This is
# necessary because the compute nodes cannot write to your $HOME
# NOTE: this directory must exist.
# Replace "<your_UCL_id>" with your UCL user ID
#$ -wd /home/regmili/Scratch/CovidsortiumTCRExpanded/

# Load the R module and run your R program
module -f unload compilers mpi gcc-libs
module load r/recommended

cd /home/regmili/Scratch/CovidsortiumTCRExpanded/
R --no-save < /home/regmili/Scratch/CovidsortiumTCRExpanded/scripts/controls_long.R > /home/regmili/myR_job.out

# Preferably, tar-up (archive) all output files to transfer them back
# to your space. This will include the R_output file above.
# tar zcvf $HOME/Scratch/R_output/files_from_job_$JOB_ID.tgz $TMPDIR

# Make sure you have given enough time for the copy to complete!
#$ -m beas #send emails at beginning, end, suspension or abortion of job