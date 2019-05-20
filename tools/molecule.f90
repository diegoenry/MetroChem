module mol

  type atom
     integer :: atmtype
     real    :: x
     real    :: y
     real    :: z
     real    :: charge
  end type atom

  type molecule
     character(len=14) :: mol_id
     integer :: natm
     real    :: charge
  end type molecule

end module mol
