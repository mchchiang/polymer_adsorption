#!/bin/bash

e_wall_start=$1
e_wall_end=$2
e_wall_inc=$3
e_bead_start=$4
e_bead_end=$5
e_bead_inc=$6
run_start=$7
run_end=$8
run_inc=${9}
dir=${10}

e_wall=$(python -c "print '%.2f' % ($e_wall_start)")
e_bead=$(python -c "print '%.2f' % ($e_bead_start)")
run=$run_start

max_jobs=10
cmd=()
jobid=0

while (( $(bc <<< "$e_wall<=$e_wall_end") ))
do
    e_bead=$(python -c "print '%.2f' % ($e_bead_start)")
    while (( $(bc <<< "$e_bead<=$e_bead_end") ))
    do
	for nbeads in 50 100 200 300 400 500 1000
	do
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		echo "Creating files for N = $nbeads EW = $e_wall EB = $e_bead run = $run"
		cmd[$jobid]="bash homopolymer_adsorption_grafted.sh $nbeads $e_wall $e_bead $run $dir"
		jobid=$(bc <<< "$jobid + 1")
	    done
	done
	e_bead=$(python -c "print '%.2f' % ($e_bead + $e_bead_inc)")
    done
    e_wall=$(python -c "print '%.2f' % ($e_wall + $e_wall_inc)")
done

# Parallel runs

total_jobs=$jobid
jobid=0

while (( $(bc <<< "$jobid < $total_jobs") ))
do
    for (( i=0; i<$max_jobs && $jobid < $total_jobs; i++))
    do
	echo "${cmd[jobid]} &"
	${cmd[jobid]} &
	jobid=$(bc <<< "$jobid + 1")
    done
    wait
done
