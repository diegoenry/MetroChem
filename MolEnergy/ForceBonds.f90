subroutine ForceBonds
use mol

!, only : bond_list%bond_i,&
!               bond_list%bond_j,&
!               bond_list%bond_k,&
!               bond_list%bond_0,&
!               molecule%x,&
!               molecule%y,&
!               molecule%z,&
!               EBOND

! Indice dos atomos e parametros de ligacoes
integer :: bond_i   ! Atom index i
integer :: bond_j   ! Atom index j
real    :: bond_k   ! Constante de mola (Kb)
real    :: bond_0   ! Distancia de equilibrio (b0)

! Coordenadas temporarias
real    :: xi
real    :: yi
real    :: zi
real    :: xj
real    :: yj
real    :: zj

! Distancias temporarias
real    :: dx
real    :: dy
real    :: dz
real    :: rij
real    :: dij

! Instantaneo
real    :: fxi
real    :: fyi
real    :: fzi

real    :: fxj
real    :: fyj
real    :: fzj


! init
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
    
    dij = rij - bond_0
    
    ! Computa o potential
    bond_potential=-bond_k * (dij)**2
    
    ! Acumulador do potential
    EBOND=EBOND+bond_potential

    ! Computa a forca
    force = - bond_k * (dij ) / 2
    
    ! Decompoe para os eixos
    fxi= force * dx
    fyi= force * dy
    fzi= force * dz

    ! Transfere para o atomo "j", com o sinal inverso.
    fxj=-fxi
    fyj=-fyi
    fzj=-fzi

    ! Acumula a força para os atomos "i" e "j"
    fx(bond_i)=fx(bond_i) + fxi
    fy(bond_i)=fy(bond_i) + fyi
    fz(bond_i)=fz(bond_i) + fzi

    fx(bond_j)=fx(bond_j) + fxj
    fy(bond_j)=fy(bond_j) + fyj
    fz(bond_j)=fz(bond_j) + fzj
    
enddo

end subroutine ForceBonds
