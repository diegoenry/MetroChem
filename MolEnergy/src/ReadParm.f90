subroutine ReadParm
use mol

character(len=32) :: line
integer	:: trash

open(1,file='butane.frcmod')


! Based on ReadPSF we need a routine to set:
! the NUMBER of different atomtypes
! the NUMBER of different bonds
! the NUMBER of different angles
! the NUMBER of different dihedrals

forcefield%num_atom_types=2
forcefield%num_bond_types=2
forcefield%num_angle_types=3
forcefield%num_dihedral_types=5


allocate(forcefield%atom_types    (forcefield%num_atom_types))
allocate(forcefield%atom_epsilon  (forcefield%num_atom_types))
allocate(forcefield%atom_sigma    (forcefield%num_atom_types))

allocate(forcefield%bond_types    (forcefield%num_bond_types))
allocate(forcefield%bond_k        (forcefield%num_bond_types))
allocate(forcefield%bond_0        (forcefield%num_bond_types))

allocate(forcefield%angle_types   (forcefield%num_angle_types))
allocate(forcefield%angle_k       (forcefield%num_angle_types))
allocate(forcefield%angle_0       (forcefield%num_angle_types))

allocate(forcefield%dihedral_types(forcefield%num_dihedral_types))
allocate(forcefield%dihedral_k    (forcefield%num_dihedral_types))
allocate(forcefield%dihedral_0    (forcefield%num_dihedral_types))
allocate(forcefield%dihedral_y    (forcefield%num_dihedral_types))

! Read BOND parameters -------------------------------------
do while (index(line,'BOND')==0)
    read(1,'(A32)',end=100) line
end do

do i=1,forcefield%num_bond_types
read(1,*)   forcefield%bond_types(i), &
            forcefield%bond_k(i),     &
            forcefield%bond_0(i) 
enddo


! Read ANGLE parameters -------------------------------------
do while (index(line,'ANGLE')==0)
   read(1,'(A32)',end=100) line
end do

do i=1,forcefield%num_angle_types
read(1,*)   forcefield%angle_types(i), &
            forcefield%angle_k(i),     &
            forcefield%angle_0(i) 
enddo


! water
close(1)

! Read Dihedral ANGLE parameters -------------------------------------
do while (index(line,'DIHE')==0)
   read(1,'(A32)',end=100) line
end do

do i=1,forcefield%num_dihedral_types
read(1,*)  forcefield%dihedral_types(i), &
            forcefield%dihedral_k(i),     &
            forcefield%dihedral_0(i),     &
            forcefield%dihedral_y(i)
enddo

100 continue
close(1)
end subroutine ReadParm
