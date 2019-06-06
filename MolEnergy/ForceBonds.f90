subroutine ForceBonds
use mol

! Indice dos atomos e parametros de ligacoes
integer :: b1   ! Atom index i
integer :: b2   ! Atom index j
real    :: bk   ! Constante de mola (Kb)
real    :: b0   ! Distancia de equilibrio (b0)

! Coordenadas temporarias
real    :: x1
real    :: y1
real    :: z1

real    :: x2
real    :: y2
real    :: z2

! Distancias temporarias
real    :: dx
real    :: dy
real    :: dz
real    :: r12
real    :: d12

! Forcas temporarias
real 	:: force
real    :: fx1
real    :: fy1
real    :: fz1
real    :: fx2
real    :: fy2
real    :: fz2


do i=1,molecule%num_bonds
    
    ! Copia indices (b1 e b2) e coordenadas (x1) dos objetos para variaveis locais
    b1=bond_list%bond_i(i)
    b2=bond_list%bond_j(i)
    bk=bond_list%bond_k(i)
    b0=bond_list%bond_0(i)

    x1=molecule%x( b1 )
    y1=molecule%y( b1 )
    z1=molecule%z( b1 )
    
    x2=molecule%x( b2 )
    y2=molecule%y( b2 )
    z2=molecule%z( b2 )
    
! Calcula as distancias --------------------------------------------
    ! O fortran eh esperto e subtrai o vetor completo
    !    ao invez de ter que fazer para cada valor, como abaixo.
    !
    ! dx(1) = xi(1) - xj(1) 
    ! dx(2) = xi(2) - xj(2) 
    ! dx(3) = xi(3) - xj(3)     
    dx=x1-x2
    dy=y1-y2
    dz=z1-z2
    
    r12= sqrt((dx)**2 + (dy)**2 + (dz)**2 )
    
    d12 = r12 - b0
    
    ! Computa o potential
   EPOT= bk * (d12)**2
    
    ! Acumulador do potential
   EBOND = EBOND + EPOT

    ! Computa a forca
    force = 2 * bk * ( d12 )
    
    ! Decompoe para os eixos
    fx1 = - force * dx
    fy1 = - force * dy
    fz1 = - force * dz

    ! Transfere para o atomo "j", com o sinal inverso.
    fx2 = - fx1
    fy2 = - fy1
    fz2 = - fz1

    ! Acumula a força para os atomos "i" e "j"
    fx(b1) = fx(b1) + fx1
    fy(b1) = fy(b1) + fy1
    fz(b1) = fz(b1) + fz1

    fx(b2) = fx(b2) + fx2
    fy(b2) = fy(b2) + fy2
    fz(b2) = fz(b2) + fz2
    
enddo

end subroutine ForceBonds
