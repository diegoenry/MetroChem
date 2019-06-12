!------------------------------------------------------------------------------
!        MetroConf - Water Optimizer
!------------------------------------------------------------------------------
! TITLE         : Water Optimizer
! PROJECT       : MetroConf
! MODULE        : WritePDB
! AFFILIATION   : Universidade Federal de Juiz de Fora
! DATE          : qui jun  6 17:32:46 -03 2019
! REVISION      : V 1.0
!> @author
!> Diego Enry Barreto Gomes
!
!> @brief
!> Escritor do arquivo PDB
!
!> @param[in] coord
!> @param[in] pdb
!---------------------------------------------------------------------------

subroutine WritePDB(coord,pdb,step)

    character(len=30),intent(in):: pdb(3) 		!< Primeiros 30 caracteres de cada linha do arquivo PDB.
	real,intent(in)  			:: coord(3,3)   !< Coordenadas do sistema       
	integer,intent(in) 			:: step         !< Passo da minimizacao
	
write(2,'("MODEL",i5)') step
do i=1,3
    write(2,'(a30,3f8.3)') pdb(i),coord(i,1),coord(i,2),coord(i,3)
enddo
write(2,'("ENDMDL")')

end subroutine writepdb
