#!/bin/bash

e_wall_start=$1
e_wall_end=$2
e_wall_inc=$3
e_bead_start=$4
e_bead_end=$5
e_bead_inc=$6
in_dir=$7
out_dir=$8

if [ ! -d $out_dir ]; then
mkdir -p $out_dir
fi

e_wall=$(python -c "print '%.2f' % ($e_wall_start)")
e_bead=$(python -c "print '%.2f' % ($e_bead_start)")
run=$run_start

binder_py="../src/BinderCumulant.py"

t_start=0
t_end=200000
t_inc=1000

bead_types=2
wall_dist=2.0

for nbeads in 50 100 200 300 400 500 1000 2000
do
    combine_file="${out_dir}/binder_hom-ads-gra_N_${nbeads}.dat"
    e_wall=$(python -c "print '%.2f' % ($e_wall_start)")
    while (( $(bc <<< "$e_wall<=$e_wall_end") ))
    do
	e_bead=$(python -c "print '%.2f' % ($e_bead_start)")
	while (( $(bc <<< "$e_bead<=$e_bead_end") ))
	do
	    echo "Doing N = $nbeads EW = $e_wall EB = $e_bead"
	    name="hom-ads-gra_N_${nbeads}_EW_${e_wall}_EB_${e_bead}"
	    bin_file="${out_dir}/bin_${name}.dat"
	    python $binder_py 0 1 $t_start $t_end $t_inc $bin_file "${in_dir}/wall-frac_${name}_run"*.dat
	    echo $e_wall $e_bead $(cat $bin_file) >> $combine_file
	    rm $bin_file
	    e_bead=$(python -c "print '%.2f' % ($e_bead + $e_bead_inc)")
	done
	e_wall=$(python -c "print '%.2f' % ($e_wall + $e_wall_inc)")
    done
done
