# Loading the molecule
mol new prot.pdb

# seperating all the protein chains.

set protein [atomselect top protein]
set chains [lsort -unique [$protein get pfrag]]
foreach chain $chains {
set sel [atomselect top "pfrag $chain"]
$sel writepdb prot_frag${chain}.pdb
}

# Making PSF from PDB using psfgen
package require psfgen 
topology ../toppar/top_all36_carb.rtf    
topology ../toppar/top_all36_lipid.rtf   
topology ../toppar/top_all36_prot.rtf    
topology ../toppar/top_all36_cgenff.rtf
topology ../toppar/top_all36_na.rtf      
topology ../toppar/top_interface.rtf 
topology ../toppar/toppar_water_ions_namd.str

foreach chain $chains {
segment $chain {pdb prot_frag${chain}.pdb}  
pdbalias residue HIS HSE
pdbalias residue HID HSD 
pdbalias residue HIP HSP
pdbalias atom ILE CD1 CD
coordpdb prot_frag${chain}.pdb $chain      
guesscoord  
}
writepdb prot_new.pdb
writepsf prot_new.psf

mol delete all          ;# clearing the queue.

# Centering the protein.

mol new prot_new.psf
mol addfile prot_new.pdb

set allatoms [atomselect top all]
$allatoms moveby [vecinvert [measure center $allatoms]]
puts [measure center $allatoms]
$allatoms writepdb centered.pdb       ;# writing coordinates of a centered protein.
$allatoms delete

mol delete all

# Solvating the protein

#mol new prot_new.psf
#mol addfile centered.pdb

package require solvate
solvate prot_new.psf centered.pdb -t 30 -o solvated

#mol delete all

# Adding ions 

#mol new solvated.psf
#mol addfile solvated.pdb

package require autoionize
autoionize -psf solvated.psf -pdb solvated.pdb -neutralize -o ionized   ;# This neutralises the system with Na and Cl 

#mol delete all

# Calculation of the cellBravisvectors
mol new ionized.psf
mol addfile ionized.pdb

set allatoms [atomselect top all]
set min_max [measure minmax $allatoms]
set vec [vecsub [lindex $min_max 1] [lindex $min_max 0]]
set origin [measure center $allatoms]
puts $vec
puts $origin
$allatoms delete

# Writing the configuration file.

puts -nonewline "For running md Enter the temperature : "
flush stdout
set temperature [gets stdin]

puts -nonewline "Enter outputName : "
flush stdout
set outputName [gets stdin]

set fo [open "config.conf" a+]
puts $fo "############################################################
ADJUSTABLE PARAMETERS                         
#############################################################

structure          ionized.psf
coordinates        ionized.pdb
set restart        0
firsttimestep      0

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm	    on
parameters          ../toppar/par_all36m_prot.prm 
parameters          ../toppar/par_all36_carb.prm 
parameters          ../toppar/par_all36_lipid.prm 
parameters          ../toppar/par_all36_cgenff.prm 
parameters          ../toppar/par_all36_na.prm 
parameters          ../toppar/par_interface.prm 
parameters          ../toppar/toppar_water_ions_namd.str 
temperature         $temperature


# Force-Field Parameters
exclude             scaled1-4
1-4scaling          1.0
cutoff              12.0
switching           on
switchdist          10.0
pairlistdist        14.0


# Integrator Parameters
timestep            2.0  ;# 2fs/step
rigidBonds          all  ;# needed for 2fs steps
nonbondedFreq       1
fullElectFrequency  2  
stepspercycle       10


# Constant Temperature Control
langevin            on    ;# do langevin dynamics
langevinDamping     1     ;# damping coefficient (gamma) of 1/ps
langevinTemp        $temperature
langevinHydrogen    off    ;# don't couple langevin bath to hydrogens


# Periodic Boundary Conditions
cellBasisVector1     [lindex $vec 0]   0.0     0.0
cellBasisVector2      0.0   [lindex $vec 1]    0.0
cellBasisVector3      0.0    0.0    [lindex $vec 2]
cellOrigin           $origin

wrapAll             on


# PME (for full-system periodic electrostatics)
PME                 yes
PMEGridSpacing      1.0

#manual grid definition
#PMEGridSizeX        45
#PMEGridSizeY        45
#PMEGridSizeZ        48


# Constant Pressure Control (variable volume)
useGroupPressure      yes ;# needed for rigidBonds
useFlexibleCell       no
useConstantArea       no

langevinPiston        on
langevinPistonTarget  1.01325 ;#  in bar -> 1 atm
langevinPistonPeriod  100.0
langevinPistonDecay   50.0
langevinPistonTemp    $temperature

# Output
outputName          $outputName

restartfreq         10000     ;# 500steps = every 1ps
dcdfreq             5000
xstFreq             5000
outputEnergies      1000
outputPressure      1000

# Fixing atom constraint ( set pdb beta column to 1 )                           
#fixedAtoms         on                                                           
#fixedAtomsFile     modified.pdb                                                 
#fixedAtomsCol      B        

# Colvars
#colvars            on
#cv configfile      colconf.in
#############################################################
## EXTRA PARAMETERS                                        ##
#############################################################


#############################################################
## EXECUTION SCRIPT                                        ##
#############################################################

# Minimization
minimize            10000
reinitvels          $temperature

run 15000000 ; #100ns 
"
close $fo

# writing SBATCH script

set f1 [open "md.sh" a+]
puts $f1 "#!/bin/bash

#SBATCH --job-name=$outputName            # Job name 
#SBATCH -N 1                    # number of nodes
#SBATCH -n 1                    # number of tasks (default: allocates 1 core per task)
#SBATCH -t 7-00:00:00           # time in d-hh:mm:s
#SBATCH -p publicgpu            # partition 
#SBATCH -q wildfire             # QOS
#SBATCH --gres=gpu:1
#SBATCH -o slurm.%j.out         # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err         # file to save job's STDERR (%j = JobId)
#SBATCH --export=NONE           # Purge the job-submitting shell environment

# Always purge modules to ensure consistent environments
module purge    
# Load required modules for job's environment
module load namd/2.13b1-cuda 

#command

namd2 config.conf > $outputName.log
"
close $f1

puts " There you have it, go run md..............."

