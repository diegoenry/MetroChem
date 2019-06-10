subroutine ReadPDB
use mol
character(len=32) :: line


!open(1,file='test.pdb')
open(1,file='ethane.pdb')

! Go to ATOM record
do while (index(line,'ATOM')==0)
    read(1,*) line
end do
backspace(1)


allocate(molecule%pdb (molecule%num_atoms))
allocate(molecule%x   (molecule%num_atoms))
allocate(molecule%y   (molecule%num_atoms))
allocate(molecule%z   (molecule%num_atoms))


do i=1,molecule%num_atoms
    read(1,'(a30,3f8.3)') molecule%pdb(i),molecule%x(i),molecule%y(i),molecule%z(i)

enddo

close(1)

end subroutine ReadPDB

