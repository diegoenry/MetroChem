program test
real                :: force(3,3) 
real                :: coord(3,3)
character(len=30)   :: line
character(len=30)   :: pdb(3)

! Read input file
call readpdb(coord,pdb)

! Open output file
open(2,file='out.pdb')


! BIG STEP LOOP
DO istep=1,100

force=0.0d0

call bforce(coord,force)

call aforce(coord,force)

call Minimize(coord,force)

! ---------------------------------
! Write out result
!write(*,'(3f15.3)') angle_deg,energy,angle_force


write(2,'("MODEL",i5)') istep
do i=1,3
    write(2,'(a30,3f8.3)') pdb(i),coord(i,1),coord(i,2),coord(i,3)
enddo
write(2,'("ENDMDL")')

ENDDO

close(2)

end program test

subroutine readpdb(coord,pdb)
    character(len=30) :: line
    real :: coord(3,3)
    character(len=30) :: pdb(3)

    ! Open input file
    open(1,file='water.pdb')
    ! Go to ATOM record
    do while (index(line,'HETATM')==0)
        read(1,*) line
    end do
    backspace(1)

    ! Read Water
    do i=1,3
        read(1,'(a30,3f8.3)')  pdb(i),coord(i,1),coord(i,2),coord(i,3)
    enddo
    
    close(1)

end subroutine readpdb


subroutine Minimize(coord,force)
    real,intent(in)     :: force(3,3)
    real,intent(inout)  :: coord(3,3)
    do i=1,3
        coord(i,1) = coord(i,1)-force(i,1)*0.00001
        coord(i,2) = coord(i,2)-force(i,2)*0.00001
        coord(i,3) = coord(i,3)-force(i,3)*0.00001
    enddo
end subroutine Minimize



subroutine aforce(x,force)
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
    

    
!
!       O(1)            b
!      /  \            / \
!     /    \          /   \
!    H(2)   H(3)     a     c
!   



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

    dot = sum ( ( ab ) &
              * ( cb ) )

    angle_rad = acos ( dot / ( v1_norm * v2_norm ) )

    angle_deg = angle_rad*180/3.1415
    
! Calcula a forca

    k=41.600
    theta_0=106.490

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

end subroutine aforce


subroutine double_cross_product(v1,v2,cp)
    implicit none
    real, intent(in)  :: v1(3) ! input
    real, intent(in)  :: v2(3) ! input
    real, intent(out) :: cp(3) ! output
    real :: ut(3) 

    ! cross product  ( v1 x v2 ) 
    !ut = cross_product(v1,v2)
    ut(1) = v1(2)*v2(3) - v1(3)*v2(2)
    ut(2) = v1(3)*v2(1) - v1(1)*v2(3)
    ut(3) = v1(1)*v2(2) - v1(2)*v2(1)
    
    ! second cross product
    ! cross product: v1 x (v1 x v2) = v1 x ut
    !tmp = cross_product(v1,ut)
    cp(1) = v1(2)*ut(3) - v1(3)*ut(2)
    cp(2) = v1(3)*ut(1) - v1(1)*ut(3)
    cp(3)=  v1(1)*ut(2) - v1(2)*ut(1)
    
    return
end subroutine double_cross_product



subroutine cross_product_3d(v1,v2,v3)
!    O produtro cruzado é o determinante da matriz
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


subroutine bforce(coord,force)
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
    
    k=371.40
    r_0=0.973
    
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
    
end subroutine bforce
