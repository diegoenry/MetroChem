!http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
!Mol2 files are written out in a free format to avoid the restrictions created by fixed format files.

! Heavly inspired by the python reader from
! Trent Balius and Sudipto Mukherjee in
! the Rizzo Research Group at Stony Brook University released in 2012
! http://docking.org/~tbalius/code/for_dock_3.7/mol2.py

program readmol2
use mol

call read_mol2()
call write_mol2()

end program readmol2

