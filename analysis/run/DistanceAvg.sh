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
multi_avg_py="../src/AverageMultiFiles.py"

while (( $(bc <<< "$e_wall<=$e_wall_end") ))
do
    e_bead=$(python -c "print '%.2f' % ($e_bead_start)")
    while (( $(bc <<< "$e_bead<=$e_bead_end") ))
    do
	echo "Doing EW = $e_wall EB = $e_bead"
	for nbeads in 2000
	do
	    name="hom-ads-gra_N_${nbeads}_EW_${e_wall}_EB_${e_bead}"
	    r2_avg_file="${out_dir}/r2_${name}_avg.dat"
	    rg_avg_file="${out_dir}/rg_${name}_avg.dat"
	    dist_avg_file="${out_dir}/dist_${name}_avg.dat"
	    python $multi_avg_py 0 1 -1 -1 $r2_avg_file "${out_dir}/dist_${name}_run"*.dat
	    python $multi_avg_py 0 2 -1 -1 $rg_avg_file "${out_dir}/dist_${name}_run"*.dat
	    paste -d" " $r2_avg_file $rg_avg_file | awk '{print $1,$2,$3,$4,$6,$7,$8}' > $dist_avg_file
	    rm $r2_avg_file $rg_avg_file
	done
	e_bead=$(python -c "print '%.2f' % ($e_bead + $e_bead_inc)")
    done
    e_wall=$(python -c "print '%.2f' % ($e_wall + $e_wall_inc)")
done
