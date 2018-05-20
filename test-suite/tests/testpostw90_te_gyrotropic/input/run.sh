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

