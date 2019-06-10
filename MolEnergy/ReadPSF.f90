subroutine ReadPSF
use mol

character(len=32) :: line
open(1,file="ethane.psf")

! Read Atoms -----------------------------------------------
do while (index(line,'NATOM')==0)
    read(1,'(A32)') line
end do
backspace(1)


read(1,*) molecule%num_atoms

allocate(molecule%id     (molecule%num_atoms))
allocate(molecule%name   (molecule%num_atoms))
allocate(molecule%type   (molecule%num_atoms))
allocate(molecule%mass   (molecule%num_atoms))
allocate(molecule%charge (molecule%num_atoms))
allocate(molecule%segid  (molecule%num_atoms))
allocate(molecule%resid  (molecule%num_atoms))
allocate(molecule%resname(molecule%num_atoms))



! Forces
allocate(fx(molecule%num_atoms))
allocate(fy(molecule%num_atoms))
allocate(fz(molecule%num_atoms))



do i=1,molecule%num_atoms
read(1,*)   molecule%id(i),      &
            molecule%segid(i),   &
            molecule%resid(i),   &
            molecule%resname(i), &
            molecule%name(i),    &
            molecule%type(i),    &
            molecule%charge(i),  &
            molecule%mass(i)
enddo



! Read Bonds -----------------------------------------------
do while (index(line,'NBOND')==0)
    read(1,'(A32)') line
end do
backspace(1)

read(1,*) molecule%num_bonds
allocate( molecule%bonds (molecule%num_bonds*2))
allocate( molecule%bond_k(molecule%num_bonds))
allocate( molecule%bond_0(molecule%num_bonds))

read(1,*) molecule%bonds


! Read Angles -----------------------------------------------
do while (index(line,'NTHETA')==0)
    read(1,'(A32)',end=100) line
end do
backspace(1)

read(1,*) molecule%num_angles

allocate( molecule%angles(molecule%num_angles*3))
read(1,*) molecule%angles



! Read Dihedrals  -----------------------------------------------
do while (index(line,'NPHI')==0)
   read(1,'(A32)',end=100) line
end do
backspace(1)

read(1,*) molecule%num_dihedrals
allocate (molecule%dihedrals(molecule%num_dihedrals*4))
read(1,*) molecule%dihedrals



100 continue
close(1)





end subroutine ReadPSF
