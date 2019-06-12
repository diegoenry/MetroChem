!------------------------------------------------------------------------------
!        MetroConf - Water Optimizer
!------------------------------------------------------------------------------
! TITLE         : Water Optimizer
! PROJECT       : MetroConf
! MODULE        : ReadPDB
! AFFILIATION   : Universidade Federal de Juiz de Fora
! DATE          : qui jun  6 17:32:46 -03 2019
! REVISION      : V 1.0
!> @author
!> Diego Enry Barreto Gomes
!
!> @brief
!> Leitor do arquivo PDB
!
!> @param[out] coord
!> @param[out] pdb
!---------------------------------------------------------------------------

subroutine ReadPDB(coord,pdb)
    character(len=30) 				:: line   		!< String temporario.
    character(len=30),intent(out) 	:: pdb(3) 		!< Primeiros 30 caracteres de cada linha do arquivo PDB.
	real,intent(out) 				:: coord(3,3)   !< Coordenadas do sistema       
		
    ! Abre o arquivo de entrada.
    open(1,file='water.pdb')
    
    ! Go to ATOM record
    do while (index(line,'HETATM')==0)
        read(1,*) line
    end do
    backspace(1)

    ! Read Water
    do i=1,3
        read(1,'(a30,3f8.3)')  pdb(i),coord(i,1),coord(i,2),coord(i,3)
    enddo
    
    close(1)
end subroutine readpdb
