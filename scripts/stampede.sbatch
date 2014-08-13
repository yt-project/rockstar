#!/bin/bash
#SBATCH -J your-job-name
#SBATCH -o your-log-file
#For a 128 (8 node) run:
#One rockstar process will be started on each node; each Rockstar process will
#start 16 subprocesses
#SBATCH -N 8 -n 8
#SBATCH -p normal
#SBATCH -t your-time-limit
#SBATCH --mail-user=your.address@email.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH -A your-account
set -x

workdir=/path/to/your/working/directory
exe=/path/to/rockstar

cd $workdir

echo Entering $(pwd)

rm auto-rockstar.cfg
$exe -c stampede.cfg >& server.dat &
perl -e 'sleep 1 while (!(-e "auto-rockstar.cfg"))'


ibrun $exe -c auto-rockstar.cfg

