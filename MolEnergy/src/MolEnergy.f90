program molenergy
use mol

verbose=.true.
!verbose=.false.

call ReadPSF
call ReadPDB
call ReadParm
call GenBondList
call GenAngleList

open(2,file='out.pdb')
call WritePDB

do step=1,200

! Zera as forcas 
fx = 0.0 ! Global, mover para fora dessa funcao
fy = 0.0 ! Global, mover para fora dessa funcao
fz = 0.0 ! Global, mover para fora dessa funcao

! Zera os potenciais
EPOT   = 0.0d0  ! Global, deve vir fora dessa funcao
EBOND  = 0.0d0
EANGLE = 0.0d0

call ForceBonds

call ForceAngles

call Minimize

!call WriteEnergy

call WritePDB
enddo


close(2)

end program molenergy

