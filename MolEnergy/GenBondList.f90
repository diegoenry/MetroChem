subroutine GenBondList
use mol

character(len=2) :: ai,aj

! We must create a bond list to optimize
!
allocate(bond_list%bond_i (molecule%num_bonds))
allocate(bond_list%bond_j (molecule%num_bonds))
allocate(bond_list%bond_k (molecule%num_bonds))
allocate(bond_list%bond_0 (molecule%num_bonds))


k=0 ! Bond counter
do i=1,molecule%num_bonds*2,2
   
   ai=trim(molecule%type( ( molecule%bonds(i) )))
   aj=trim(molecule%type( ( molecule%bonds(i+1) )))
      
    do j=1,forcefield%num_bond_types

        if  (ai == forcefield%bond_types(j)(1:2) &
        .and. &
             aj == forcefield%bond_types(j)(4:5)) &
        then
                ! Assign bond parameters
                k=k+1
                bond_list%bond_i(k)=molecule%bonds(i)
                bond_list%bond_j(k)=molecule%bonds(i+1)
                bond_list%bond_k(k)=forcefield%bond_k(j)
                bond_list%bond_0(k)=forcefield%bond_0(j)

        end if
    enddo
enddo
! write it here


if (verbose) then
    write(*,'(i5,1x,"!NBONDS")') molecule%num_bonds
    do k=1,molecule%num_bonds

        write(*,'(2i3,2f8.3)')  &
                    bond_list%bond_i(k),&
                    bond_list%bond_j(k),&
                    bond_list%bond_k(k),&
                    bond_list%bond_0(k)
    enddo
end if

! 
end subroutine GenBondList
