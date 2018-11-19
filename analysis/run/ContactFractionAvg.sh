#!/bin/bash

#nbeads=$1
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

# Average scripts
avg_py="../src/TimeAverage.py"
multi_avg_py="../src/AverageMultiFiles.py"

t_start=100000
t_end=200000
t_inc=1000


while (( $(bc <<< "$e_wall<=$e_wall_end") ))
do
    e_bead=$(python -c "print '%.2f' % ($e_bead_start)")
    while (( $(bc <<< "$e_bead<=$e_bead_end") ))
    do
	echo "Doing EW = $e_wall EB = $e_bead"
	combined_avg_file="${out_dir}/wall-frac_hom-ads-gra_EW_${e_wall}_EB_${e_bead}_avg.dat"
	> $combined_avg_file
	for nbeads in 50 100 200 300 400 500 1000
	do
	name="hom-ads-gra_N_${nbeads}_EW_${e_wall}_EB_${e_bead}"
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    frac_file="${in_dir}/wall-frac_${name}_run_${run}.dat"
	    avg_file="${out_dir}/wall-frac_${name}_run_${run}_avg.dat"
	    python $avg_py 0 1 $t_start $t_end $t_inc $frac_file $avg_file
	    data=$(cat $avg_file)
	done
	multi_avg_file="${out_dir}/wall-frac_${name}_avg.dat"
	python $multi_avg_py -1 0 -1 -1 $multi_avg_file "${out_dir}/wall-frac_${name}_run"*_avg.dat
	echo $nbeads $(cat $multi_avg_file) >> $combined_avg_file
	rm "${out_dir}/wall-frac_${name}_run"*_avg.dat
	rm $multi_avg_file
	done
	e_bead=$(python -c "print '%.2f' % ($e_bead + $e_bead_inc)")
    done
    e_wall=$(python -c "print '%.2f' % ($e_wall + $e_wall_inc)")
done
