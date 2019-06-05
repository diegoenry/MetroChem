module mol
implicit none
! atom
! atomtypes
! natom
! nbonds
! rotatable bonds

! From Tripos MOL2 documentation
! atom_id atom_name x y z atom_type [subst_id
!        [subst_name [charge [status_bit]]]]

! bonds
! bond_id origin_atom_id target_atom_id bond_type

! 1 = single
! 2 = double
! 3 = triple
! am = amide
! ar = aromatic
! du = dummy
! un = unknown (cannot be determined from the parameter tables)
! nc = not connected

type molecule
    character(len=128)           :: mol_name
    integer                      :: num_atoms
    integer	                     :: num_bonds
    integer	                     :: num_feat
    integer	                     :: num_sets
    character(len=20)            :: mol_type
    character(len=20)            :: charge_type
end type molecule

type atom
    integer,allocatable          :: atom_id(:)
    character(len=4),allocatable :: atom_name(:)
    character(len=4),allocatable :: atom_type(:)
    real,allocatable             :: x(:)
    real,allocatable             :: y(:)
    real,allocatable             :: z(:)
    integer,allocatable          :: subst_id(:)
    character(len=6),allocatable :: subst_name(:)
    real,allocatable             :: charge(:)
end type atom

type bond
    integer,allocatable :: a1(:)
    integer,allocatable :: a2(:)
    integer,allocatable :: type(:)
end type bond

type(molecule) :: mol_info
type(atom)     :: coordinates
type(bond)     :: bonds

end module mol
