subroutine GenDihedralList
use mol
integer	:: trash

!
!V(dihedral) = Kchi(1 + cos(n(chi) - delta))
!
!Kchi: kcal/mole
!n: multiplicity
!delta: degrees

! Atomos 1 até 4
character(len=3) :: a1,a2,a3,a4

! We must create a dihedral list to optimize
allocate( dihedral_list%dihedral_1 (molecule%num_dihedrals))
allocate( dihedral_list%dihedral_2 (molecule%num_dihedrals))
allocate( dihedral_list%dihedral_3 (molecule%num_dihedrals))
allocate( dihedral_list%dihedral_k (molecule%num_dihedrals))
allocate( dihedral_list%dihedral_0 (molecule%num_dihedrals))

k=0 ! dihedral counter
do i=1,molecule%num_dihedrals*4,4
   
   a1=trim(molecule%type( ( molecule%dihedrals(i) )))
   a2=trim(molecule%type( ( molecule%dihedrals(i+1) )))
   a3=trim(molecule%type( ( molecule%dihedrals(i+2) )))
   a3=trim(molecule%type( ( molecule%dihedrals(i+3) )))
      
    do j=1,forcefield%num_dihedral_types

        if  (a1 == forcefield%dihedral_types(j)(1:2) &
        .and. &
             a2 == forcefield%dihedral_types(j)(4:5) &
        .and. &
             a3 == forcefield%dihedral_types(j)(7:8) &
        .or. &
             a3 == forcefield%dihedral_types(j)(1:2) &
        .and. &
             a2 == forcefield%dihedral_types(j)(4:5) &
        .and. &
             a1 == forcefield%dihedral_types(j)(7:8) ) &
        then
                ! Assign dihedral parameters
                k=k+1
                dihedral_list%dihedral_1(k)=molecule%dihedrals(i)
                dihedral_list%dihedral_2(k)=molecule%dihedrals(i+1)
                dihedral_list%dihedral_3(k)=molecule%dihedrals(i+2)
				dihedral_list%dihedral_4(k)=molecule%dihedrals(i+3)
                dihedral_list%dihedral_k(k)=forcefield%dihedral_k(j)
                dihedral_list%dihedral_0(k)=forcefield%dihedral_0(j)
    
        end if
    enddo
enddo
! write it here
! 


write(*,*) "Verify dihedrals"

if (verbose) then
    write(*,'(i5,1x,"!NdihedralS")') molecule%num_dihedrals
    do k=1,molecule%num_dihedrals

        write(*,'(3i3,2f8.3)')  &
                    dihedral_list%dihedral_1(k),&
                    dihedral_list%dihedral_2(k),&
                    dihedral_list%dihedral_3(k),&
                    dihedral_list%dihedral_k(k),&
                    dihedral_list%dihedral_0(k)
    enddo
end if


end subroutine GenDihedralList
