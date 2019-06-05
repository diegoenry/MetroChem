subroutine EnergyBonds
use mol
!integer :: ai
!integer :: aj
character(len=2) :: ai,aj
real :: xi
real :: yi
real :: zi
real :: xj
real :: yj
real :: zj
real :: r2
real :: dij
real :: rij
real :: bond_potential

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
!                write(*,*) "[BOND]", &
!                    molecule%bonds(i),&
!                    molecule%bonds(i+1), &
!                    forcefield%bond_types(j),&
!                    forcefield%bond_k(j),&
!                    forcefield%bond_0(j)
                
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
! 


do i=1,molecule%num_bonds
    
    xi=molecule%x(bond_list%bond_i(i))    
    yi=molecule%y(bond_list%bond_i(i))
    zi=molecule%z(bond_list%bond_i(i))
    
    xj=molecule%x(bond_list%bond_j(i))
    yj=molecule%y(bond_list%bond_j(i))
    zj=molecule%z(bond_list%bond_j(i))
    
    r2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
    dij = 1.d0/sqrt(r2)
    rij = r2*dij
    
    bond_potential=bond_list%bond_k(i) * (rij * rij)
    write(*,*) bond_list%bond_i(i),bond_list%bond_j(i),rij,bond_potential
enddo 


! Don't forget to Setup set neighbourlist


end subroutine EnergyBonds
