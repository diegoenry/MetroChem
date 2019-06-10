module mol
implicit none


integer :: step   ! Simulation step
logical :: verbose

type topology

    !Read from .PSF
    integer                      :: num_atoms
    integer	                     :: num_bonds
    integer	                     :: num_angles
    integer	                     :: num_dihedrals

    integer,allocatable          :: bonds(:)
    real*8,allocatable             :: bond_k(:)
    real*8,allocatable             :: bond_0(:)
    
    
    integer,allocatable          :: angles(:)
    integer,allocatable          :: dihedrals(:)
    
    integer,allocatable          :: id(:)
    character(len=4),allocatable :: name(:)
    character(len=6),allocatable :: type(:)
    real*8,allocatable             :: mass(:)
    real*8,allocatable             :: charge(:)    
    character(len=6),allocatable :: segid(:)
    integer,allocatable          :: resid(:)
    character(len=6),allocatable :: resname(:)
    
    ! Read from PDB file
    character(len=30),allocatable:: pdb(:)
    real*8,allocatable             :: x(:)
    real*8,allocatable             :: y(:)
    real*8,allocatable             :: z(:)

end type topology

type(topology) :: molecule



type parameters
    integer :: num_atom_types=2
    integer :: num_bond_types=2
    integer :: num_angle_types=3
    integer :: num_dihedral_types=4

    character(len=2),allocatable :: atom_types(:)
    real*8,allocatable :: atom_epsilon(:)
    real*8,allocatable :: atom_sigma(:)
    
    character(len=5),allocatable :: bond_types(:)
    real*8,allocatable :: bond_k(:)
    real*8,allocatable :: bond_0(:)

    character(len=8),allocatable :: angle_types(:)
    real*8,allocatable :: angle_k(:)
    real*8,allocatable :: angle_0(:)
    
    character(len=11),allocatable :: dihedral_types(:)
    real*8,allocatable :: dihedral_k(:)
    real*8,allocatable :: dihedral_0(:)
    real*8,allocatable :: dihedral_y(:)

end type parameters

type(parameters) :: forcefield


type bonds
    character(len=5),allocatable 	:: bond_types(:)
    integer,allocatable 			:: bond_i(:)
    integer,allocatable 			:: bond_j(:)
    real*8,allocatable    			:: bond_k(:)
    real*8,allocatable    			:: bond_0(:)
end type bonds

type(bonds) :: bond_list


type angles
    character(len=5),allocatable 	:: angle_types(:)
    integer,allocatable 			:: angle_1(:)
    integer,allocatable 			:: angle_2(:)
    integer,allocatable 			:: angle_3(:)
    real*8,allocatable    			:: angle_k(:)
    real*8,allocatable   				:: angle_0(:)
end type angles

type(angles) :: angle_list


! FORCE
real*8, allocatable :: fx(:)
real*8, allocatable :: fy(:)
real*8, allocatable :: fz(:)

! Potential
real*8 :: EBOND
real*8 :: EANGLE


end module mol
