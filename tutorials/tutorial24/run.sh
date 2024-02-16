#!/bin/bash


NPROC=12
mpirun="mpirun -np $NPROC"
mpirun=""

PWpath=""
WANpath="../../"


pwcmd="$mpirun  ${PWpath}pw.x "
pw2wancmd="$mpirun  ${PWpath}pw2wannier90.x"
wan90cmd="${WANpath}wannier90.x"
postw90cmd="$mpirun  ${WANpath}postw90.x"


cp input/* .

cp Te.win-0 Te.win

$pwcmd < Te.scf  > Te.out.scf

rm -r Te.save/wfc*.dat

$pwcmd  < Te.nscf  > Te.out.nscf

$wan90cmd -pp Te

$pw2wancmd < Te.pw2wan > pw2wan.out

rm Te.wfc*
rm Te.save/wfc*.dat

$wan90cmd Te

cat >Te.win <<EOF
gyrotropic=true
gyrotropic_task=-C-dos-D0-Dw-K
fermi_energy_step=0.0025
fermi_energy_min=5.8
fermi_energy_max=6.2
gyrotropic_freq_step=0.0025
gyrotropic_freq_min=0.0
gyrotropic_freq_max=0.1
gyrotropic_smr_fixed_en_width=0.01
gyrotropic_smr_max_arg=5
gyrotropic_degen_thresh=0.001
gyrotropic_box_b1=0.2 0.0 0.0
gyrotropic_box_b2=0.0 0.2 0.0
gyrotropic_box_b3=0.0 0.0 0.2
gyrotropic_box_center=0.33333 0.33333 0.5
gyrotropic_kmesh=50 50 50
EOF

cat input/Te.win >> Te.win
$postw90cmd Te

cp Te.win Te.win-gyrotropic
mv Te.wpout Te.wpout-gyrotropic

cat >Te.win <<EOF
gyrotropic=true
gyrotropic_task=-NOA
fermi_energy=5.95
gyrotropic_freq_step=0.0025
gyrotropic_freq_min=0.0
gyrotropic_freq_max=0.3
gyrotropic_smr_fixed_en_width=0.01
gyrotropic_smr_max_arg=5
gyrotropic_band_list=4-9
gyrotropic_kmesh=50 50 50
EOF

cat input/Te.win >> Te.win
$postw90cmd Te
cp Te.win Te.win-NOA
mv Te.wpout Te.wpout-NOA


