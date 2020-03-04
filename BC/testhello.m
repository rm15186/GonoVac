%% SUBMISSION SCRIPT NOT MATLAB
%script that tells blue crystal what to do, calls my run avg, outputs some stuff? has all the commands?

#!/bin/bash
#! Sample PBS submission script for a MATLAB job
#! Requesting resource (processors and wall-clock time): #! nodes=1:ppn=1 indicates a single processor.
#! nodes=1:ppn=16 would request a whole node for BCp3. #! 02:30:00 indicates 02 hours and 30 minutes
#PBS -l nodes=1:ppn=1,walltime = 00:30:00
#! change the working directory (default is home directory) 
cd $PBS_O_WORKDIR
#! Record some useful job details in the output file echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This jobs runs on the following nodes: 
echo `cat $PBS_NODEFILE | uniq`
#! add the MATLAB module (as per BCp3) 
module add apps/matlab-r2019a 
#! NB have limited MATLAB to a single thread 
options="-nodesktop -noFigureWindows -singleCompThread"
#! Run MATLAB in batch mode 
matlab $options -r hello3 
#! time ./my-serial-program
#! it says hello2 is undefined