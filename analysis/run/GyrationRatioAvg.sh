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

# Average scripts
avg_py="../src/TimeAverage.py"
multi_avg_py="../src/AverageMultiFiles.py"

t_start=100000
t_end=200000
t_inc=1000

for nbeads in 50 100 200 300 400 500 1000
do
    e_wall=$(python -c "print '%.2f' % ($e_wall_start)")
    combined_avg_file="${out_dir}/gyr_hom-ads-gra_N_${nbeads}_avg.dat"
    > $combined_avg_file
    while (( $(bc <<< "$e_wall<=$e_wall_end") ))
    do
	e_bead=$(python -c "print '%.2f' % ($e_bead_start)")
	while (( $(bc <<< "$e_bead<=$e_bead_end") ))
	do
	    echo "Doing N = $nbeads EW = $e_wall EB = $e_bead"
	    name="hom-ads-gra_N_${nbeads}_EW_${e_wall}_EB_${e_bead}"
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		# Compute the average of the gyration tensor components
		thermo_file="${in_dir}/thermo/thermo_${name}_run_${run}.dat"
		rg0_file="${out_dir}/rg0_${name}_run_${run}.dat"
		rg_file="${out_dir}/rg_${name}_run_${run}_avg.dat"
		rgh_file="${out_dir}/rgh_${name}_run_${run}_avg.dat"
		rgv_file="${out_dir}/rgv_${name}_run_${run}_avg.dat"
		awk '{print $1,$5,0.5*($6+$7),$8}' $thermo_file > $rg0_file
		python $avg_py 0 1 $t_start $t_end $t_inc $rg0_file $rg_file
		python $avg_py 0 2 $t_start $t_end $t_inc $rg0_file $rgh_file
		python $avg_py 0 3 $t_start $t_end $t_inc $rg0_file $rgv_file
	    done
	    rg_avg_file="${out_dir}/rg_${name}_avg.dat"
	    rgh_avg_file="${out_dir}/rgh_${name}_avg.dat"
	    rgv_avg_file="${out_dir}/rgv_${name}_avg.dat"
	    gyr_avg_file="${out_dir}/gyr_${name}_avg.dat"
	    python $multi_avg_py -1 0 2 -1 $rg_avg_file "${out_dir}/rg_${name}_run"*_avg.dat
	    python $multi_avg_py -1 0 2 -1 $rgh_avg_file "${out_dir}/rgh_${name}_run"*_avg.dat
	    python $multi_avg_py -1 0 2 -1 $rgv_avg_file "${out_dir}/rgv_${name}_run"*_avg.dat
	    paste -d" " $rg_avg_file $rgh_avg_file $rgv_avg_file > $gyr_avg_file
	    rm "${out_dir}/rg0_${name}_run"*.dat
	    rm "${out_dir}/rg_${name}_run"*_avg.dat
	    rm "${out_dir}/rgh_${name}_run"*_avg.dat
	    rm "${out_dir}/rgv_${name}_run"*_avg.dat
	    rm $rg_avg_file $rgh_avg_file $rgv_avg_file
	    
	    # Compute the ratio and its error
	    ratio=$(awk '{
rgv=$7**0.5; 
rgh=$4**0.5;
rgverr=$9;
rgherr=$6;
ratio=rgv/rgh;
ratioerr=((rgverr**2/rgv**2+rgherr**2/rgh**2)**0.5)*ratio;
printf("%.5f %.5f",ratio,ratioerr)
}' $gyr_avg_file)
	    
	    echo $e_wall $e_bead $ratio >> $combined_avg_file
	    rm $gyr_avg_file
	    e_bead=$(python -c "print '%.2f' % ($e_bead + $e_bead_inc)")
	done
	e_wall=$(python -c "print '%.2f' % ($e_wall + $e_wall_inc)")
    done
done
