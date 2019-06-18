gfortran -O3  \
   src/Modules.f90 \
   src/ReadPSF.f90 \
   src/ReadPDB.f90 \
   src/ReadParm.f90 \
   src/GenBondList.f90 \
   src/GenAngleList.f90 \
   src/EnergyBonds.f90 \
   src/EnergyAngles.f90 \
   src/ForceBonds.f90 \
   src/ForceAngles.f90 \
   src/Minimize.f90 \
   src/WritePDB.f90 \
   src/WriteEnergy.f90 \
   src/MolEnergy.f90 \
   -o MolEnergy
