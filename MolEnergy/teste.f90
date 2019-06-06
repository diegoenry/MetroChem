program test
real :: a(3), b(3) ,c(3)
real :: ba(3),bc(3)
real :: uv_dot
real :: v1_norm
real :: v2_norm
real :: angle_cos
real :: angle
real :: angle_deg
real :: angle_force

real, parameter :: pi= 4.0d0*atan2(1.0d0,1.0d0)
real :: dcp(3)
real :: pa(3), pc(3)
real :: fa(3), fc(3), fb(3) 

character(len=30) :: line
real :: coord(3,3)
character(len=30) :: pdb(3)

! Read input file
call readpdb(coord,pdb)

! Open output file
open(2,file='out.pdb')

! Write 1st frame
write(2,'("MODEL",i5)') 1
do i=1,3
    write(2,'(a30,3f8.3)') pdb(i),coord(:,i)
enddo
write(2,'("ENDMDL")')


! BIG STEP LOOP
DO istep=1,10

! this is because of the BOND list. 
! H - O - H  ->  (2 - 1 - 3)
a=coord(:,2)
b=coord(:,1) 
c=coord(:,3)







! ANGLES ---------------------------------------------------------------
ba=b-a
bc=b-c

uv_dot  = dot_product ( ba, bc )

v1_norm = sqrt ( dot_product ( ba, ba ) )
v2_norm = sqrt ( dot_product ( bc, bc ) )

angle_cos = uv_dot / v1_norm / v2_norm
angle = acos ( angle_cos )    

angle_deg = angle * 180/pi

call aforce(angle_deg,energy,angle_force)

! ----------------------------------------------------------------------
! Redistribute forces to each atom
!
! ---------------------------------
! Force on Atom A
call double_cross_product(ba, bc, dcp)

! norm(tmp)
norma = sqrt(dcp(1)*dcp(1) + dcp(2)*dcp(2) + dcp(3)*dcp(3))
pa(1) = dcp(1)/norma
pa(2) = dcp(2)/norma
pa(3) = dcp(3)/norma    

force(1,1) = (angle_force / v1_norm)*pa(1) + fa(1,1)
force(1,2) = (angle_force / v1_norm)*pa(2) + fa(1,2)
force(1,3) = (angle_force / v1_norm)*pa(3) + fa(1,3)

! ---------------------------------
! Force on Atom C
call double_cross_product(bc, ba, dcp)

! norm(tmp)
norma = sqrt(dcp(1)*dcp(1) + dcp(2)*dcp(2) + dcp(3)*dcp(3))
pc(1) = dcp(1)/norma
pc(2) = dcp(2)/norma
pc(3) = dcp(3)/norma    

fc(1) = (angle_force / v2_norm)*pc(1)
fc(2) = (angle_force / v2_norm)*pc(2)
fc(3) = (angle_force / v2_norm)*pc(3)

! ---------------------------------
! Force on Atom B
fb=-fa-fc



! ---------------------------------
! Write out result
write(*,'(3f15.3)') angle_deg,energy,angle_force
! write(*,*) fa
! write(*,*) fc
! write(*,*) fb

call Steepest(coord,fa,fb,fc)

write(2,'("MODEL",i5)') istep+1
do i=1,3
    write(2,'(a30,3f8.3)') pdb(i),coord(:,i)
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
        read(1,'(a30,3f8.3)')  pdb(i),coord(:,i)
    enddo
    
    close(1)

end subroutine readpdb


subroutine Steepest(coord,fa,fb,fc)
    real,intent(in)     :: fa(3), fc(3), fb(3) 
    real,intent(inout)  :: coord(3,3)
    do i=1,3
        coord(:,1) = coord(:,1)-fb*0.001
        coord(:,2) = coord(:,2)-fa*0.001
        coord(:,3) = coord(:,3)-fc*0.001
    enddo
end subroutine Steepest


subroutine bforce(r,energy,force)
    real, intent(in)  :: r
    real, intent(out) :: energy
    real, intent(out) :: force
    real :: k
    
    k=371.40
    r_0=0.973
    ! Computa o potential
    energy= - k * (r - r_0 )**2
    
    ! Computa a forca -2k(theta-theta_0) 
    force = - 2 * k * (r - r_0 )  

end subroutine bforce



subroutine aforce(angle,energy,force)
    real, intent(in)  :: angle
    real, intent(out) :: energy
    real, intent(out) :: force
    real :: k
    
    k=41.600
    theta_0=106.490
    ! Computa o potential
    energy= - k * (angle - theta_0 )**2
    
    ! Computa a forca -2k(theta-theta_0) 
    force = - 2 * k * (angle - theta_0)   

end subroutine aforce


subroutine double_cross_product(v1,v2,cp)
    implicit none
    real, intent(in)  :: v1(3) ! input
    real, intent(in)  :: v2(3) ! input
    real, intent(out) :: cp(3) ! output
    real :: ut(3) ! output

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

