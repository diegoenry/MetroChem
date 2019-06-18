program readihe
type dihedrals
    character(len=11),allocatable 	:: dihedral_types(:)
    integer,allocatable 			:: dihedral_1(:)
    integer,allocatable 			:: dihedral_2(:)
    integer,allocatable 			:: dihedral_3(:)
    integer,allocatable 			:: dihedral_4(:)
    ! triple cosine dihedrals
    real*8,allocatable    			:: dihedral_k1(:)
    real*8,allocatable    			:: dihedral_k2(:)
    real*8,allocatable    			:: dihedral_k3(:)
    real*8,allocatable   			:: dihedral_0(:)
end type dihedrals

type(dihedrals) :: dihedral_list

type parameters
    integer :: num_dihedral_types
    character(len=11),allocatable :: dihedral_types(:)
    real*8,allocatable :: dihedral_k(:)
    real*8,allocatable :: dihedral_k(:)
    real*8,allocatable :: dihedral_k(:)
    real*8,allocatable :: dihedral_0(:)
    real*8,allocatable :: dihedral_y(:)
end type parameters

type(parameters) :: forcefield

type dihedral_vector

end type dihedral_vector



character(len=32) 	:: line
integer 			:: trash 

forcefield%num_dihedral_types=5

! We must create a dihedral list to optimize
allocate(forcefield%dihedral_types(forcefield%num_dihedral_types))
allocate(forcefield%dihedral_k    (forcefield%num_dihedral_types))
allocate(forcefield%dihedral_0    (forcefield%num_dihedral_types))
allocate(forcefield%dihedral_y    (forcefield%num_dihedral_types))

allocate( dihedral_list%dihedral_1 (forcefield%num_dihedral_types))
allocate( dihedral_list%dihedral_2 (forcefield%num_dihedral_types))
allocate( dihedral_list%dihedral_3 (forcefield%num_dihedral_types))
allocate( dihedral_list%dihedral_k (forcefield%num_dihedral_types))
allocate( dihedral_list%dihedral_0 (forcefield%num_dihedral_types))

open(1,file="example/butane/butane.frcmod")


! Read Dihedral ANGLE parameters -------------------------------------
do while (index(line,'DIHE')==0)
   read(1,'(A32)',end=100) line
end do

write(*,*) line, forcefield%num_dihedral_types

do i=1,forcefield%num_dihedral_types
read(1,*)   forcefield%dihedral_types(i)
write(*,*) forcefield%dihedral_types(i)
STOP
enddo


do i=1,forcefield%num_dihedral_types
read(1,*)   forcefield%dihedral_types(i), &
			trash, &
            forcefield%dihedral_k(i),     &
            forcefield%dihedral_0(i),     &
            forcefield%dihedral_y(i)
enddo

100 continue
close(1)

end program readihe
