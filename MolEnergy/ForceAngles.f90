subroutine ForceAngles
use mol

! Particle defining the angle
real :: p1(3)
real :: p2(3)
real :: p3(3)

real :: fa(3)
real :: fb(3)
real :: fc(3)


real :: uv_dot

real :: ut(3) ! temporary array to store cross product

real :: v1(3)
real :: v1_norm ! v1_norm = |ab|


real :: v2(3)
real :: v2_norm ! v2_norm = |bc|


real :: v3(3)
real :: v3_norm

real :: v1v2(3) ! cross product v1,v2
real :: v2v1(3) ! cross product v2,v1

real :: angle
real :: angle_cos

real :: dangle

real :: angle_potential
real :: angle_force


! pa-> is the normalized vector in plane abc, orthogonal do ba->
real :: pa(3)

! pc-> is the normalized vector in plane abc, orthogonal do cb->
real :: pc(3

)
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
    v1=p2-p1   ! ba
    v2=p2-p3   ! bc
    v3=p3-p2   ! cb
    
    
    uv_dot  = dot_product ( v1, v2 )

    v1_norm = sqrt ( dot_product ( v1, v1 ) )

    v2_norm = sqrt ( dot_product ( v2, v2 ) )

    angle_cos = uv_dot / v1_norm / v2_norm

    angle = acos ( angle_cos )    
    
    dangle = angle - angle_list%angle_0(i)
    
    ! Computa o potential
    angle_potential=-angle_list%angle_k(i) * (dangle )**2
    
    ! Acumulador do potential
    EANGLE=EANGLE+angle_potential

    ! Computa a forca -2k(theta-theta_0) 
    angle_force = -2 * angle_list%angle_k(i) * (dangle )    
    
!---------------------
! Forces on atom 1 (a) 

! 
    cp = double_cross_product(v1,v2)
    
    ! norm(tmp)
    norma = sqrt(cp(1)*cp(1) + cp(2)*cp(2) + cp(3)*cp(3))
    pa(1) = cp(1)/norma
    pa(2) = cp(2)/norma
    pa(3) = cp(3)/norma    
    
    !r = normalize(tmp)
    
    ! check
    !rnorm = sqrt(dot_product(r,r))
    !rnorm = 1
    
    fa(1) = (angle_force /v1_norm)*pa(1)
    fa(2) = (angle_force /v1_norm)*pa(2)
    fa(3) = (angle_force /v1_norm)*pa(3)


    
!---------------------
! Forces on atom 3 (a) 

    ! cross product  v1 x v2
    tmp1 = v1(2)*v2(3) - v1(3)*v2(2)
    tmp2 = v1(3)*v2(1) - v1(1)*v2(3)
    tmp3 = v1(1)*v2(2) - v1(2)*v2(1)
    
    !ut = cross_product(v1,v2)
    
    ut(1) = tmp1
    ut(2) = tmp2
    ut(3) = tmp3
    
    ! cross product: v2 x (v1 x v2) = v1 x ut
    tmp1 = v2(2)*ut(3) - v2(3)*ut(2)
    tmp2 = v2(3)*ut(1) - v2(1)*ut(3)
    tmp3 = v2(1)*ut(2) - v2(2)*ut(1)
    
    !tmp = cross_product(v1,ut)
    
    ! norm(tmp)
    norma = sqrt(tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3)
    tmp1 = tmp1/norma
    tmp2 = tmp2/norma
    tmp3 = tmp3/norma
    
    !r = normalize(tmp)
    r(1) = tmp1
    r(2) = tmp2
    r(3) = tmp3 
    
    ! check
    !rnorm = sqrt(dot_product(r,r))
    !rnorm = 1
    
    fc(1) = (angle_force /v2_norm)*r(1)
    fc(2) = (angle_force /v2_norm)*r(2)
    fc(3) = (angle_force /v2_norm)*r(3)

    
    write(*,*)  angle_list%angle_1(i),&
                angle_list%angle_2(i),&
                angle_list%angle_3(i),&
                angle * 180/3.1415,angle_potential,angle_force
    ! Agora complica.
    ! Decomposicao para cada atomo e cada eixo.
    

enddo 

end subroutine ForceAngles


!function cross_product(va,vb)
!
!return
!end function

