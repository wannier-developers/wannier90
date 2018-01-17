#!/bin/bash


calcQE=0



PWpath=""
WANpath="../../../"


wan90cmd="${WANpath}wannier90.x"
postw90cmd="${WANpath}postw90.x"
pwcmd="${PWpath}pw.x"
pw2wancmd="${PWpath}pw2wannier90.x"
if ((calcQE==1))
then
    echo "calcualting QE"
    cp input/* .
    $pwcmd < Te.scf  > Te.out.scf
    rm -r Te.save/wfc*.dat
    $pwcmd  < Te.nscf  > Te.out.nscf
    $wan90cmd -pp Te
    $pw2wancmd < Te.pw2wan > pw2wan.out
    rm Te.wfc*
    rm Te.save/wfc*.dat


else
    echo "using pre-calculated wannier files"
    ln -s wannier_files/* .
    cp input/Te.win .
fi


$wan90cmd Te

$postw90cmd Te

exit

cat >Te.win <<EOF
gyrotropic=true
gyrotropic_task=-C-dos-D0-Dw-K-NOA
fermi_energy_step=2
fermi_energy_min=2
fermi_energy_max=10
gyrotropic_freq_step=0.05
gyrotropic_freq_min=0.0
gyrotropic_freq_max=0.1
gyrotropic_smr_fixed_en_width=0.1
gyrotropic_smr_max_arg=5
gyrotropic_degen_thresh=0.001
gyrotropic_box_b1=0.2 0.0 0.0
gyrotropic_box_b2=0.0 0.2 0.0
gyrotropic_box_b3=0.0 0.0 0.2
gyrotropic_box_center=0.33333 0.33333 0.5
gyrotropic_kmesh=5 5 5
EOF

cat input/Te.win >> Te.win
