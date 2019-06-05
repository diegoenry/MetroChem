module mol
implicit none


type topology

    !Read from .PSF
    integer                      :: num_atoms
    integer	                     :: num_bonds
    integer	                     :: num_angles
    integer	                     :: num_dihedrals

    integer,allocatable          :: bonds(:)
    real,allocatable             :: bond_k(:)
    real,allocatable             :: bond_0(:)
    
    
    integer,allocatable          :: angles(:)
    integer,allocatable          :: dihedrals(:)
    
    integer,allocatable          :: id(:)
    character(len=4),allocatable :: name(:)
    character(len=6),allocatable :: type(:)
    real,allocatable             :: mass(:)
    real,allocatable             :: charge(:)    
    character(len=6),allocatable :: segid(:)
    integer,allocatable          :: resid(:)
    character(len=6),allocatable :: resname(:)
    
    ! Read from PDB file
    character(len=30),allocatable :: pdb(:)
    real,allocatable             :: x(:)
    real,allocatable             :: y(:)
    real,allocatable             :: z(:)

end type topology

type(topology) :: molecule



type parameters
    integer :: num_atom_types=2
    integer :: num_bond_types=2
    integer :: num_angle_types=3
    integer :: num_dihedral_types=4

    character(len=2),allocatable :: atom_types(:)
    real,allocatable :: atom_epsilon(:)
    real,allocatable :: atom_sigma(:)
    
    character(len=5),allocatable :: bond_types(:)
    real,allocatable :: bond_k(:)
    real,allocatable :: bond_0(:)

    character(len=8),allocatable :: angle_types(:)
    real,allocatable :: angle_k(:)
    real,allocatable :: angle_0(:)
    
    character(len=11),allocatable :: dihedral_types(:)
    real,allocatable :: dihedral_k(:)
    real,allocatable :: dihedral_0(:)
    real,allocatable :: dihedral_y(:)

end type parameters

type(parameters) :: forcefield


type bonds
    character(len=5),allocatable :: bond_types(:)
    integer,allocatable :: bond_i(:)
    integer,allocatable :: bond_j(:)
    real,allocatable :: bond_k(:)
    real,allocatable :: bond_0(:)
end type bonds

type(bonds) :: bond_list



! FORCE
real, allocatable :: fx(:)
real, allocatable :: fy(:)
real, allocatable :: fz(:)



end module mol
