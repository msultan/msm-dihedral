place metrics_degree.py in your MSM library directory, mine is:
/home/mlawrenz/Programs/epd-7.1-1-rh5-x86_64/lib/python2.7/site-packages/msmbuilder

place diehdrals_chi.py in your MSM library geomtery folder:
/home/mlawrenz/Programs/epd-7.1-1-rh5-x86_64/lib/python2.7/site-packages/msmbuilder/geometry/

TestDihedrals.py key lines are:
mymetric=metrics.Dihedral(metric='euclidean', p=2, angles=['chi'])
ptraj=mymetric.prepare_trajectory(t, degrees=True)

this creates the metric object for Dihedrals, with only phi, psi, chi, and
combos of these as options. you pass in combos using 
'phi/psi'
'phi psi'
'phi','psi'


the script prints out *.dat files that are compared with g_chi output in this
directory, and a matrix of size len(residues)xlen(residues)xtrajlength filled
with the dihedral values in degrees. Note that the residue ID's here differ from gromacs after
GLU195-LYS237 due to the missing loop. gromacs continues the numbers with LYS197.

if you want you can modify dihedrals_chi.py to include chi2/3, just modify 
_get_indices_chi1 to include the correct atom indices.
