#!/bin/bash

AMethod=FEALL

#/home/tquah/FreeEnergyDomainSize/Prototype/FE_Fluct.py -d "C_*/CL/chi*" -f "model1_operators.dat" -kl 0 2 -kt C chi -n 100 -ed "/home/tquah/DataStore/2Period/" -e CL_DISPhase.dat -T --OP -ep -eft --OPFile "model1_OrientationPersistenceOP.dat" -at ${AMethod} -w 1000
for f in *; do
    echo $f
    if [ -d "$f" ]; then
       
       ~/FreeEnergyDomainSize/Prototype/FixCorruptFiles.py -d "${f}/C_*/CHI_*/A_0.0/ZETA_0.0/F*/DISPHASE/" -i operators.dat -o operators_clean.dat
       ~/FreeEnergyDomainSize/Prototype/FE_Fluct.py -d "${f}/C_*/CHI_*/A_*/ZETA_*/F*/DISPHASE/" -f operators_clean.dat -kl 1 2 3 4 5 -kt C CHI A ZETA F -ed . -e ${f}_FE.dat -ep -eft -at ${AMethod} -n 100

    fi
done


