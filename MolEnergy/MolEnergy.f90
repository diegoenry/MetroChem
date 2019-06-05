program molenergy
use mol

call ReadPSF
call ReadPDB
call ReadParm
call GenBondList

!call EnergyBonds

open(2,file='out.pdb')
call WritePDB(1)

do l=2,50
call ForceBonds
call Steepest
call WritePDB(l)
enddo


close(2)

end program molenergy
