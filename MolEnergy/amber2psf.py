import parmed as pmd

mol = pmd.load_file('butane.prmtop','butane.rst7')
mol.write_psf('butane.psf',vmd=True)
quit()

