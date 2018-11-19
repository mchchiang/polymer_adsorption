#!/bin/bash

# An input script for setting parameters in the LAMMPS config file

# Read in and set parameters
nbeads=$1         # number of beads
e_wall=$2         # wall energy
e_bead=$3         # bead energy
run=$4            # trial number
run_dir=$5        # run directory

# Formatting input arguments
e_wall=$(python -c "print '%.2f' % ($e_wall)")
e_bead=$(python -c "print '%.2f' % ($e_bead)")

# Determine the box size from volume fraction
vol_frac=0.001
L=$(python -c "import math; print '%.13f' % (math.pi*$nbeads/6./$vol_frac)**(1./3.)")
lx=$L
ly=$L
lz=$L
buffer=2 # Buffer from the wall for the initial config

xlo=$(python -c "print -$lx/2.0")
xhi=$(python -c "print $lx/2.0")
ylo=$(python -c "print -$ly/2.0")
yhi=$(python -c "print $ly/2.0")
zlo=$(python -c "print -$lz/2.0")
zhi=$(python -c "print $lz/2.0")

# Other energy parameters
e_harm=1000
e_soft=100

# Init config (default is random walk)
gen_polymer_exe="../bin/exe/GenGraftedRandomWalk"

# Make execution directory
sim_name="hom-ads-gra_N_${nbeads}_EW_${e_wall}_EB_${e_bead}_run_${run}"
run_dir="${run_dir}/EW_${e_wall}_EB_${e_bead}/N_${nbeads}/${sim_name}"

if [ ! -d $run_dir ]; then
mkdir -p $run_dir
fi

# Set output file names
init_file="init_${sim_name}.in"
restart_file="restart_${sim_name}"
prep_outfile="prep_${sim_name}.lammpstrj"
prep_endfile="prep_${sim_name}.out"
run_outfile="run_${sim_name}.lammpstrj"
run_endfile="run_${sim_name}.out"
map_file="${sim_name}.lammpsmap"

# Generate polymer
${gen_polymer_exe} $nbeads $lx $ly $lz $buffer "${run_dir}/${init_file}" "${run_dir}/${map_file}"

max_seed=100000

restart_freq=1000

# A function for generating random numbers
function get_rand(){
    rand=$(python -c "import random, sys; print random.randint(0,$max_seed)")
    echo $rand
}

# Prep/equilibration times
thermo_printfreq=1000
prep_printfreq=1000
run_printfreq=1000
prep_seed=$(get_rand)
run_seed=$(get_rand)
prep_time_1=2000 # 2000 run using soft + harmonic
prep_time_2=4000 # 4000 run using soft + fene
prep_time_3=4000 # 4000 run using lj/cut + fene
run_time=200000  # simulation run time  

delta_t=0.01       # time step size in Brownian time units

# Interaction energy parameters
# LJ potentials
sigma=1.0
cutoff=$(python -c "print '%.13f' % (1.8*$sigma)")

# Normalisation (ensure minimum of potential is actually epsilon)
norm=$(python -c "print '%.13f' % (1.0 + 4.0*(($sigma/$cutoff)**12-($sigma/$cutoff)**6))")

e_wall_norm=$(python -c "print '%.13f' % ($e_wall/$norm)")
e_bead_norm=$(python -c "print '%.13f' % ($e_bead/$norm)")

# Convert all time values to simulation time (i.e. rescale by delta t)
restart_freq=$(bc <<< "$restart_freq/$delta_t")
prep_time_1=$(bc <<< "$prep_time_1/$delta_t")
prep_time_2=$(bc <<< "$prep_time_2/$delta_t")
prep_time_3=$(bc <<< "$prep_time_3/$delta_t")
run_time=$(bc <<< "$run_time/$delta_t")
thermo_printfreq=$(bc <<< "$thermo_printfreq/$delta_t")
prep_printfreq=$(bc <<< "$prep_printfreq/$delta_t")
run_printfreq=$(bc <<< "$run_printfreq/$delta_t")

# Create the lammps command file based on template
lammps_file="${sim_name}.lam"
qsub_file="qsub_${sim_name}.sh"
file="${run_dir}/${lammps_file}"
#qsub_file="${run_dir}/${qsub_file}"
log_file="${sim_name}.log"

# Choose template depending on the type of wall used
cp homopolymer_adsorption_grafted.lam $file

# Replace macros in template with input values
sed -i -- "s/INIT_FILE/${init_file}/g" $file
sed -i -- "s/RESTART_FILE/${restart_file}/g" $file

sed -i -- "s/XLO/${xlo}/g" $file
sed -i -- "s/XHI/${xhi}/g" $file
sed -i -- "s/YLO/${ylo}/g" $file
sed -i -- "s/YHI/${yhi}/g" $file
sed -i -- "s/ZLO/${zlo}/g" $file
sed -i -- "s/ZHI/${zhi}/g" $file

sed -i -- "s/RESTART_FREQ/${restart_freq}/g" $file

sed -i -- "s/THERMO_PRINTFREQ/${thermo_printfreq}/g" $file
sed -i -- "s/PREP_PRINTFREQ/${prep_printfreq}/g" $file
sed -i -- "s/RUN_PRINTFREQ/${run_printfreq}/g" $file

sed -i -- "s/PREP_SEED/${prep_seed}/g" $file
sed -i -- "s/PREP_TIME_1/${prep_time_1}/g" $file
sed -i -- "s/PREP_TIME_2/${prep_time_2}/g" $file
sed -i -- "s/PREP_TIME_3/${prep_time_3}/g" $file
sed -i -- "s/RUN_SEED/${run_seed}/g" $file
sed -i -- "s/RUN_TIME/${run_time}/g" $file

sed -i -- "s/PREP_OUTFILE/${prep_outfile}/g" $file
sed -i -- "s/PREP_ENDFILE/${prep_endfile}/g" $file

sed -i -- "s/RUN_OUTFILE/${run_outfile}/g" $file
sed -i -- "s/RUN_ENDFILE/${run_endfile}/g" $file

sed -i -- "s/DELTA_T/${delta_t}/g" $file

sed -i -- "s/EHARM/${e_harm}/g" $file
sed -i -- "s/ESOFT/${e_soft}/g" $file

sed -i -- "s/EWALL/${e_wall_norm}/g" $file
sed -i -- "s/EBEAD/${e_bead_norm}/g" $file

sed -i -- "s/SIGMA/${sigma}/g" $file
sed -i -- "s/CUTOFF/${cutoff}/g" $file
