subroutine EnergyAngles
use mol

! Particle defining the angle
real :: p1(3)
real :: p2(3)
real :: p3(3)

real :: angle_rad
real :: dot

real :: angle_potential

Eangle=0.0d0

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

    ! 3D vectors
    v1norm = sqrt ( sum ( ( p1(1:3) - p2(1:3) )**2 ) )
    if ( v1norm == 0.0D+00 ) then
        angle_rad = 0.0D+00
        cycle
    end if

    v2norm = sqrt ( sum ( ( p2(1:3) - p3(1:3) )**2 ) )
    if ( v2norm == 0.0D+00 ) then
        angle_rad = 0.0D+00
        cycle
    end if
  
    dot = dot_product(v1norm,v2norm)

    angle_cos = dot / v1norm / v2norm

    angle = cos ( angle_cos )
    
    write(*,*) angle_rad * 180/3.1415
    
    !angle_potential=
    !Eangle=Eangle+angle_potential
enddo 

!write(*,*) 'Eangle = ',Eangle
! Don't forget to Setup set neighbourlist

end subroutine EnergyAngles
