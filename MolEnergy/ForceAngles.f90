subroutine ForceAngles
use mol

! Particle defining the angle
real :: p1(3)
real :: p2(3)
real :: p3(3)

real :: uv_dot

real :: v1(3)
real :: v1_norm

real :: v2(3)
real :: v2_norm

real :: angle
real :: angle_cos

real :: dangle

real :: angle_potential
real :: angle_force

EANGLE=0.0d0

do i=1,molecule%num_angles
    
    p1(1)=molecule%x(angle_list%angle_1(i))    
    p1(2)=molecule%y(angle_list%angle_1(i))
    p1(3)=molecule%z(angle_list%angle_1(i))

    p2(1)=molecule%x(angle_list%angle_2(i))
    p2(2)=molecule%y(angle_list%angle_2(i))
    p2(3)=molecule%z(angle_list%angle_2(i))

    p3(1)=molecule%x(angle_list%angle_3(i))
    p3(2)=molecule%y(angle_list%angle_3(i))
    p3(3)=molecule%z(angle_list%angle_3(i))

    ! 3D vectors j->i, j->k
    v1=p2-p1
    v2=p2-p3
    
    uv_dot = dot_product ( v1, v2 )

    v1_norm = sqrt ( dot_product ( v1, v1 ) )

    v2_norm = sqrt ( dot_product ( v2, v2 ) )

    angle_cos = uv_dot / v1_norm / v2_norm

    angle = acos ( angle_cos )    
    
    dangle = angle - angle_list%angle_0(i)
    
    ! Computa o potential
    angle_potential=-angle_list%angle_k(i) * (dangle )**2
    
    ! Acumulador do potential
    EANGLE=EANGLE+angle_potential

    ! Computa a forca
    angle_force = -angle_list%angle_k(i) * (dangle ) / 2    
    
    
    write(*,*)  angle_list%angle_1(i),&
                angle_list%angle_2(i),&
                angle_list%angle_3(i),&
                angle * 180/3.1415,angle_potential,angle_force
    ! Agora complica.
    ! Decomposicao para cada atomo e cada eixo.
    

enddo 

!write(*,*) 'Eangle = ',Eangle
! Don't forget to Setup set neighbourlist

end subroutine ForceAngles
