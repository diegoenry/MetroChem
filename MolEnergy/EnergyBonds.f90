subroutine EnergyBonds
use mol
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
real :: EBOND=0.0d0

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
    
    bond_potential=-bond_list%bond_k(i) * (rij * rij)
    EBOND=EBOND+bond_potential
    write(*,*) bond_list%bond_i(i),bond_list%bond_j(i),rij,bond_potential
enddo 

write(*,*) 'EBOND = ',EBOND
! Don't forget to Setup set neighbourlist

end subroutine EnergyBonds
