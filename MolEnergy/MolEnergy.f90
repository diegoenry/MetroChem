program molenergy
use mol

!verbose=.true.
verbose=.false.

call ReadPSF
call ReadPDB
call ReadParm
call GenBondList
call GenAngleList
call EnergyAngles

stop

!call EnergyBonds

open(2,file='out.pdb')
call WritePDB(1)

do step=2,50
call ForceBonds
call WriteEnergy
call WritePDB
call Steepest

enddo


close(2)

end program molenergy

