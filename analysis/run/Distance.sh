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
in_dir=${10}
out_dir=${11}

if [ ! -d $out_dir ]; then
mkdir -p $out_dir
fi

e_wall=$(python -c "print '%.2f' % ($e_wall_start)")
e_bead=$(python -c "print '%.2f' % ($e_bead_start)")
run=$run_start

distance_exe="../bin/exe/Distance"

max_jobs=8
cmd=()
jobid=0

t_start=100000
t_end=200000
t_inc=1000


while (( $(bc <<< "$e_wall<=$e_wall_end") ))
do
    e_bead=$(python -c "print '%.2f' % ($e_bead_start)")
    while (( $(bc <<< "$e_bead<=$e_bead_end") ))
    do
	for nbeads in 50 100 200 300 400 500 1000 2000
	do
	    nbeadsp1=$(bc <<< "$nbeads+1") # To include the anchoring bead
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		name="hom-ads-gra_N_${nbeads}_EW_${e_wall}_EB_${e_bead}_run_${run}"
		pos_file="${in_dir}/position/pos_${name}.dat"
		out_file="${in_dir}/siminfo/run_${name}.out"
		dist_file="${out_dir}/dist_${name}.dat"
		
		# Determine the box size
		lx=$(awk '{if(NR==10){printf("%.13f",$2-$1)}}' $out_file)
		ly=$(awk '{if(NR==11){printf("%.13f",$2-$1)}}' $out_file)
		lz=$(awk '{if(NR==12){printf("%.13f",$2-$1)}}' $out_file)
		echo $lx $ly $lz
		
		cmd[$jobid]="$distance_exe $nbeadsp1 $lx $ly $lz $t_start $t_end $t_inc $pos_file $dist_file"
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
