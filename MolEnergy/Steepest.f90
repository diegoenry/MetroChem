subroutine Steepest
use mol
    
do i=1,molecule%num_atoms
    molecule%x(i)=molecule%x(i)+(fx(i)*0.001)
    molecule%y(i)=molecule%y(i)+(fy(i)*0.001)
    molecule%z(i)=molecule%z(i)+(fz(i)*0.001)
enddo

end subroutine Steepest
