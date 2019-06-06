!------------------------------------------------------------------------------
!        MetroConf - Optimizer
!------------------------------------------------------------------------------
! TITLE         : Optimizer
! PROJECT       : MetroConf
! MODULE        : ForceAngles
! AFFILIATION   : Universidade Federal de Juiz de Fora
! DATE          : qui jun  6 17:32:46 -03 2019
! REVISION      : V 1.0
!> @author
!> Diego Enry Barreto Gomes
!> @brief
!> Calcula a energia, força devido as interacoes ligadas do tipo Angulo (ANGLES)
!
!> @param[in]  coord
!> @param[out] force
!------------------------------------------------------------------------------
!
!       O(1)            b
!      /  \            / \
!     /    \          /   \
!    H(2)   H(3)     a     c
!

subroutine ForceAngles
use mol

! Indice dos atomos e parametros de ligacoes
integer :: a1   ! Atom index i
integer :: a2   ! Atom index j
integer :: a3   ! Atom index k
real    :: ak   ! Constante de mola (Kb)
real    :: a0   ! Distancia de equilibrio (b0)

! Coordenadas temporarias
real :: x1(3)
real :: x2(3)
real :: x3(3)

! Valores para o calculo do angulo
real :: v1_norm
real :: v2_norm
real :: dot
real :: angle_rad
real :: angle_deg

! force
real :: angle_force
real :: angle_energy
real :: v3(3)
real :: ut(3)
real :: norm
real :: p(3)
real :: ab(3)
real :: ba(3)
real :: bc(3)
real :: cb(3)

! Calculo dos angulos --------------------------------------------------
do i=1,molecule%num_angles

	a1 = angle_list%angle_1(i)  ! Atom index i (a) (1)
	a2 = angle_list%angle_2(i)	! Atom index j (b) (2)
	a3 = angle_list%angle_3(i)	! Atom index k (c) (3)
    ak = angle_list%angle_k(i)	! Constante de mola (K)
	a0 = angle_list%angle_0(i)	! Distancia de equilibrio (b0)
	
! [ Passo 1] Copia as coordenadas
    x1(1) = molecule%x(a1)    
    x1(2) = molecule%y(a1)
    x1(3) = molecule%z(a1)

    x2(1) = molecule%x(a2)
    x2(2) = molecule%y(a2)
    x2(3) = molecule%z(a2)

    x3(1) = molecule%x(a3)
    x3(2) = molecule%y(a3)
    x3(3) = molecule%z(a3)

! [ Passo 2 ] Cria os vetores de distancia AB, BC
    ab = x1-x2
    ba = x2-x1
    bc = x2-x3
    cb = x3-x2

! [ Passo 3 ] normaliza o vetor ba, e o vetor cb
! da no mesmo fazer o ba e bc.
!    v1_norm = sqrt ( sum ( ( ab )**2 ) )
!    v2_norm = sqrt ( sum ( ( cb )**2 ) )
!	dot = sum ( ( ab ) * ( cb ) )	

    dot     = dot_product ( ab, cb )
    
    v1_norm = sqrt ( dot_product ( ab, ab ) )
    
    v2_norm = sqrt ( dot_product ( cb, cb ) )

    angle_rad = acos ( dot / ( v1_norm * v2_norm ) )
    
    angle_deg = angle_rad * 180 / 3.14159265359  !aproximated pi

!  write(*,*) dot, angle_rad,angle_deg  

! Calcula a forca ------------------------------------------------------

    ! Computa o potential
    angle_energy = ak * (angle_deg - a0 )**2
    
    EANGLE = EANGLE + angle_energy
    
    ! Computa a forca -2k(theta-theta_0) 
    angle_force = - 2 * ak * (angle_deg - a0) 

!   write(*,*) a1,a2,a3,angle_deg,(angle_deg - a0), a0,ak,angle_force
    
! ------------------------------------------------------------------
! Redistribute forces to each atom
!
! ---------------------------------
! Force on Atom A

    call cross_product_3d(ba, bc, ut)
    call cross_product_3d(ba, ut, v3)
    
    norm = sqrt(v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3))
    p(1) = v3(1)/norm
    p(2) = v3(2)/norm
    p(3) = v3(3)/norm    

    fx(a1) = - (angle_force / v1_norm) * p(1) + fx(a1)
    fy(a1) = - (angle_force / v1_norm) * p(2) + fy(a1)
    fz(a1) = - (angle_force / v1_norm) * p(3) + fz(a1)

! Force on Atom c

    call cross_product_3d(ba, bc, ut)
    call cross_product_3d(cb, ut, v3)
    
    norm = sqrt(v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3))
    p(1) = v3(1)/norm
    p(2) = v3(2)/norm
    p(3) = v3(3)/norm    

    fx(a3) = - (angle_force / v2_norm) * p(1) + fx(a3)
    fy(a3) = - (angle_force / v2_norm) * p(2) + fy(a3)
    fz(a3) = - (angle_force / v2_norm) * p(3) + fz(a3)

! Force on Atom B
    fx(a2) = - fx(a1) - fx(a3) + fx(a2)
    fy(a2) = - fy(a1) - fy(a3) + fy(a2)
    fz(a2) = - fz(a1) - fz(a3) + fz(a2)

!   write(*,*) fx(a1),fy(a1), fz(a1)
!   write(*,*) fx(a2),fy(a2), fz(a2)
!   write(*,*) fx(a3),fy(a3), fz(a3)
!   write(*,*) 'EANGLE = ',EANGLE

enddo



end subroutine ForceAngles


! Math routines
subroutine cross_product_3d(v1,v2,v3)
!    O produto cruzado é o determinante da matriz
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
    implicit none
    real, intent(in)    :: v1(3)
    real, intent(in)    :: v2(3)
    real, intent(out)   :: v3(3)

    v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
    v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
    v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

    return
end subroutine cross_product_3d
