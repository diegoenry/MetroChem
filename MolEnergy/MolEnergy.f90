program molenergy
use mol


call ReadPSF
call ReadPDB
call ReadParm
call GenBondList

call EnergyBonds

end program molenergy
