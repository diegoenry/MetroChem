program WaterOptimize
    character(len=30) 	:: line   		!< String temporario.
    character(len=30) 	:: pdb(3) 		!< Primeiros 30 caracteres de cada linha do arquivo PDB.
	real 				:: coord(3,3)   !< Coordenadas do sistema     
	real 				:: force(3,3)   !< Forcas do sistema     
	integer 			:: step         !< Passo da minimizacao

! Read input file
call ReadPDB(coord,pdb)

! Open output file
open(2,file='out.pdb')

! BIG STEP LOOP
DO istep=1,100

force=0.0d0

call BondForce(coord,force)

call AngleForce(coord,force)

call Minimize(coord,force)

call WritePDB(coord,pdb,istep)

! ---------------------------------
! Write out result
!write(*,'(3f15.3)') angle_deg,energy,angle_force

ENDDO

close(2)

end program WaterOptimize



