!------------------------------------------------------------------------------
!        MetroConf - Water Optimizer
!------------------------------------------------------------------------------
! TITLE         : Water Optimizer
! PROJECT       : MetroConf
! MODULE        : aForce
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

subroutine AngleForce(x,force)
implicit none
    real, intent(in)    :: x(3,3)
    real, intent(inout) :: force(3,3)
    real :: xi(3),xj(3),xk(3)
    real :: v1_norm,v2_norm
    real :: dot,angle_rad,angle_deg
! force
    real :: k,theta_0
    real :: angle_force,angle_energy
    real :: v3(3),ut(3)
    real :: norm
    real :: p(3)
    real :: ab(3),ba(3),bc(3),cb(3)
    

    ! Valores gerados de acordo com o General Amber Force Field v2.0
    ! Nao sao os ideais para a molecula de agua.
    k=41.600
    theta_0=106.490


! Calcula o angulo
    xi=x(2,:) 
    xj=x(1,:) 
    xk=x(3,:)

    ab = xi-xj
    ba = xj-xi
    bc = xj-xk
    cb = xk-xj
    
    v1_norm = sqrt ( sum ( ( ab )**2 ) )
    v2_norm = sqrt ( sum ( ( cb )**2 ) )

    dot = sum ( ( ab ) * ( cb ) )

    angle_rad = acos ( dot / ( v1_norm * v2_norm ) )

    angle_deg = angle_rad*180/3.14159265359  !aproximated pi
            write(*,*) dot, angle_rad,angle_deg

! Calcula a forca

    ! Computa o potential
    angle_energy= k * (angle_deg - theta_0 )**2
    
    ! Computa a forca -2k(theta-theta_0) 
    angle_force =  2 * k * (angle_deg - theta_0) 

    !write(*,*) angle_rad,angle_deg, angle_force


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

    force(2,1) = -(angle_force / v1_norm) * p(1) + force(2,1)
    force(2,2) = -(angle_force / v1_norm) * p(2) + force(2,2)
    force(2,3) = -(angle_force / v1_norm) * p(3) + force(2,3)

! Force on Atom c

    call cross_product_3d(ba, bc, ut)
    call cross_product_3d(cb, ut, v3)
    
    norm = sqrt(v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3))
    p(1) = v3(1)/norm
    p(2) = v3(2)/norm
    p(3) = v3(3)/norm    

    force(3,1) = -(angle_force / v2_norm) * p(1) + force(3,1)
    force(3,2) = -(angle_force / v2_norm) * p(2) + force(3,2)
    force(3,3) = -(angle_force / v2_norm) * p(3) + force(3,3)

! Force on Atom B
    force(1,1) = - force(2,1) - force(3,1) + force(1,1)
    force(1,2) = - force(2,2) - force(3,2) + force(1,2)
    force(1,3) = - force(2,3) - force(3,3) + force(1,3)

end subroutine AngleForce
