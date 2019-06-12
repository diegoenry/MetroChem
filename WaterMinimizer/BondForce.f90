!------------------------------------------------------------------------------
!        MetroConf - Water Optimizer
!------------------------------------------------------------------------------
! TITLE         : Water Optimizer
! PROJECT       : MetroConf
! MODULE        : bforce
! AFFILIATION   : Universidade Federal de Juiz de Fora
! DATE          : qui jun  6 17:32:46 -03 2019
! REVISION      : V 1.0
!> @author
!> Diego Enry Barreto Gomes
!
!> @brief
!> Calcula a energia, força devido as interacoes ligadas (BONDS)
!
!> @param[out] coord
!> @param[out] pdb
!------------------------------------------------------------------------------

subroutine BondForce(coord,force)
    real,intent(inout)  :: force(3,3)
    real,intent(inout)  :: coord(3,3)
    real :: r
    real :: energy
    real :: k
    real :: xi,yi,zi
    real :: xj,yj,zj
    real :: dx,dy,dz
    real :: fxi,fyi,fzi
    real :: fxj,fyj,fzj
    real :: rij
    real :: dij
    real :: bond_force
    
    ! Valores gerados de acordo com o General Amber Force Field v2.0
    ! Note que nao sao os ideais para a molecula de agua.
    k=371.40  ! Potencial de mola
    r_0=0.973 ! Distancia de equilibrio
    
    do i=2,3
    
    ! Oxygen
        xi=coord(1,1)
        yi=coord(1,2)
        zi=coord(1,3)

    ! Each hydrogen
        xj=coord(i,1)
        yj=coord(i,2)
        zj=coord(i,3)
        
        ! Calcula as distancias
        dx=xi-xj
        dy=yi-yj
        dz=zi-zj
        
        ! Calcula o raio
        rij= sqrt( (dx)**2 + (dy)**2 + (dz)**2 )
    
        ! Calcula a diferenca em relacao a r_0
        dij = rij - r_0

        ! Computa a forca
        bond_force = 2 * k * (dij )
        !write(*,*) 1,i,rij,dij,bond_force

        ! Decompoe para os eixos
        fxi = bond_force * dx
        fyi = bond_force * dy
        fzi = bond_force * dz

        ! Transfere para o atomo "j", com o sinal inverso.
        fxj = -fxi
        fyj = -fyi
        fzj = -fzi
        
        ! Acumula a força para os atomos "i" e "j"
        force(1,1) = force(1,1) + fxi
        force(1,2) = force(1,2) + fyi
        force(1,3) = force(1,3) + fzi
        
        force(i,1) = force(i,1) + fxj
        force(i,2) = force(i,2) + fyj
        force(i,3) = force(i,3) + fzj
        
    enddo
    
end subroutine BondForce
