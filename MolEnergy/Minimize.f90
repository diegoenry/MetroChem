subroutine Minimize
use mol

real*8 :: n,fa,fb,fc
    
do i=1,molecule%num_atoms
    n = sqrt(fx(i)*fx(i) + fy(i)*fy(i) + fz(i)*fz(i))
    molecule%x(i)=molecule%x(i)+((fx(i)/n)*0.01)
    molecule%y(i)=molecule%y(i)+((fy(i)/n)*0.01)
    molecule%z(i)=molecule%z(i)+((fz(i)/n)*0.01)
enddo

end subroutine Minimize

