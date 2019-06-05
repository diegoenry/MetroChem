subroutine ForceBonds
use mol

integer           :: bond_i,bond_j  ! Atom index i, j
real :: bond_k
real :: bond_0
!tmp coordinates
real :: xi
real :: yi
real :: zi
real :: xj
real :: yj
real :: zj
!distances
real :: dx
real :: dy
real :: dz
! Instantaneo
real :: fxi
real :: fyi
real :: fzi

real :: fxj
real :: fyj
real :: fzj


fx=0.0d0
fy=0.0d0
fz=0.0d0

EBOND=0.0d0

do i=1,molecule%num_bonds
    
    ! Save to local variables
    bond_i=bond_list%bond_i(i)
    bond_j=bond_list%bond_j(i)
    bond_k=bond_list%bond_k(i)
    bond_0=bond_list%bond_0(i)

    xi=molecule%x(bond_i)    
    yi=molecule%y(bond_i)
    zi=molecule%z(bond_i)
    
    xj=molecule%x(bond_j)
    yj=molecule%y(bond_j)
    zj=molecule%z(bond_j)
    
    
    ! Calcula as distancias
    dx=xi-xj
    dy=yi-yj
    dz=zi-zj
    rij= sqrt((dx)**2 + (dy)**2 + (dz)**2 )
    
    ! Computa o potential
    bond_potential=-bond_k * (rij - bond_0)**2
    
    ! Acumulador do potential
    EBOND=EBOND+bond_potential

    fxi= - bond_k * (dx) / 2
    fyi= - bond_k * (dy-bond_0) / 2
    fzi= - bond_k * (dz-bond_0) / 2

    
    fxj=-fxi
    fyj=-fyi
    fzj=-fzi

    ! Acumulador da força
    fx(bond_i)=fx(bond_i) + fxi
    fy(bond_i)=fy(bond_i) + fyi
    fz(bond_i)=fz(bond_i) + fzi

    fx(bond_j)=fx(bond_j) + fxj
    fy(bond_j)=fy(bond_j) + fyj
    fz(bond_j)=fz(bond_j) + fzj
    
enddo
write(*,*) 'EBOND = ',EBOND,'R = ',rij,'fxi = ',fxi,'fxj = ',fxj,bond_list%bond_k(1)



end subroutine ForceBonds
