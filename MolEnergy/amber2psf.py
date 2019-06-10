import parmed as pmd

mol = pmd.load_file('ethane.prmtop','ethane.rst7')
mol.write_psf('ethane.psf',vmd=True)
quit()

