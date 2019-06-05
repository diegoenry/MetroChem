subroutine read_mol2()
use mol

character(len=72) :: line
character(len=8)  :: trash

open(1,file='butane.mol2')

! Go to molecule
do while (index(line,'@<TRIPOS>MOLECULE')==0)
    read(1,*) line
end do

! Read molecule name
read(1,*) mol_info%mol_name

! Read num_atoms and nbonds
read(1,*) mol_info%num_atoms, mol_info%num_bonds

! Allocate atom array to store coordinates, and charges
allocate(coordinates%atom_id   (mol_info%num_atoms ))
allocate(coordinates%atom_name (mol_info%num_atoms ))
allocate(coordinates%atom_type (mol_info%num_atoms ))
allocate(coordinates%x         (mol_info%num_atoms ))
allocate(coordinates%y         (mol_info%num_atoms ))
allocate(coordinates%z         (mol_info%num_atoms ))
allocate(coordinates%subst_id  (mol_info%num_atoms ))
allocate(coordinates%subst_name(mol_info%num_atoms ))
allocate(coordinates%charge    (mol_info%num_atoms ))


! Go to ATOM record
do while (index(line,'@<TRIPOS>ATOM')==0)
    read(1,*) line
end do

do i=1,mol_info%num_atoms
read(1,*) coordinates%atom_id(i),   &
          coordinates%atom_name(i), &
          coordinates%x(i),         &
          coordinates%y(i),         &
          coordinates%z(i),         &
          coordinates%atom_type(i), &
          coordinates%subst_id(i),  &
          coordinates%subst_name(i),&
          coordinates%charge(i)  
enddo 
end subroutine read_mol2


subroutine write_mol2()
use mol

write(*,"(A)") "@<TRIPOS>MOLECULE"

! Read molecule name
write(*,*) trim(mol_info%mol_name)

! Read num_atoms and nbonds
write(*,'(5i5)') mol_info%num_atoms, mol_info%num_bonds,1,0,0

write(*,"(A)") "SMALL"
write(*,"(A)") "NO_CHARGES"


100 format(i5,3x,a5,f8.4,3x,f8.4,3x,f8.4,5x,a5,i4,7x,a6,1x,f8.4)
do i=1,mol_info%num_atoms

write(*,100) coordinates%atom_id(i),   &
             coordinates%atom_name(i), &
             coordinates%x(i),         &
             coordinates%y(i),         &
             coordinates%z(i),         &
             coordinates%atom_type(i), &
             coordinates%subst_id(i),  &
             coordinates%subst_name(i),&
             coordinates%charge(i)
enddo

end subroutine write_mol2

