subroutine GenAngleList
use mol

character(len=3) :: a1,a2,a3

! We must create a angle list to optimize
!
allocate(angle_list%angle_1 (molecule%num_angles))
allocate(angle_list%angle_2 (molecule%num_angles))
allocate(angle_list%angle_3 (molecule%num_angles))
allocate(angle_list%angle_k (molecule%num_angles))
allocate(angle_list%angle_0 (molecule%num_angles))

k=0 ! angle counter
do i=1,molecule%num_angles*3,3
   
   a1=trim(molecule%type( ( molecule%angles(i) )))
   a2=trim(molecule%type( ( molecule%angles(i+1) )))
   a3=trim(molecule%type( ( molecule%angles(i+2) )))
      
    do j=1,forcefield%num_angle_types

        if  (a1 == forcefield%angle_types(j)(1:2) &
        .and. &
             a2 == forcefield%angle_types(j)(4:5) &
        .and. &
             a3 == forcefield%angle_types(j)(7:8) &
        .or. &
             a3 == forcefield%angle_types(j)(1:2) &
        .and. &
             a2 == forcefield%angle_types(j)(4:5) &
        .and. &
             a1 == forcefield%angle_types(j)(7:8) ) &
        then
                ! Assign angle parameters
                k=k+1
                angle_list%angle_1(k)=molecule%angles(i)
                angle_list%angle_2(k)=molecule%angles(i+1)
                angle_list%angle_3(k)=molecule%angles(i+2)
                angle_list%angle_k(k)=forcefield%angle_k(j)
                angle_list%angle_0(k)=forcefield%angle_0(j)
    
        end if
    enddo
enddo
! write it here
! 



write(*,*) "Verify angles"

if (verbose) then
    write(*,'(i5,1x,"!NANGLES")') molecule%num_angles
    do k=1,molecule%num_angles

        write(*,'(3i3,2f8.3)')  &
                    angle_list%angle_1(k),&
                    angle_list%angle_2(k),&
                    angle_list%angle_3(k),&
                    angle_list%angle_k(k),&
                    angle_list%angle_0(k)
    enddo
end if


end subroutine GenAngleList
