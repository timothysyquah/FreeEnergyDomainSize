#!/bin/bash

/home/tquah/FreeEnergyDomainSize/Prototype/FE_Fluct.py -d "C*/CL/chi*" -f "model1_operators.dat" -kl 0 2 -kt C chi -n 100 -ed "/home/tquah/DataStore/" -e CL_DISPhase.dat
/home/tquah/FreeEnergyDomainSize/Prototype/FE_Fluct.py -d "C*/CL/chi*" -f "model2_operators.dat" -kl 0 2 -kt C chi -n 100 -ed "/home/tquah/DataStore/" -e CL_LAMPhase.dat
/home/tquah/FreeEnergyDomainSize/Prototype/FE_Fluct.py -d "C*/PartialFTS/chi*" -f "model1_operators.dat" -kl 0 2 -kt C chi -n 100 -ed "/home/tquah/DataStore/" -e PFTS_DISPhase.dat
/home/tquah/FreeEnergyDomainSize/Prototype/FE_Fluct.py -d "C*/PartialFTS/chi*" -f "model2_operators.dat" -kl 0 2 -kt C chi -n 100 -ed "/home/tquah/DataStore/" -e PFTS_LAMPhase.dat
