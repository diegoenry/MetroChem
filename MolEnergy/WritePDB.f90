subroutine WritePDB
use mol

write(2,"(A6,i5)") "MODEL ",step
do i=1,molecule%num_atoms
    write(2,'(a30,3f8.3)') molecule%pdb(i),molecule%x(i),molecule%y(i),molecule%z(i)
enddo
write(2,"(A6)") "ENDMDL"    

end subroutine WritePDB

