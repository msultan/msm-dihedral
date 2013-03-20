import numpy
import os
from msmbuilder import metrics_test as metrics
from msmbuilder import Trajectory, Conformation

conf=Conformation.LoadFromPDB('protein_common_residues.pdb')
frames=numpy.where((conf['ResidueNames']!='ALA')&(conf['ResidueNames']!='GLY'))
resids=dict()
resids['phi']=set(conf['ResidueID'])
resids['psi']=set(conf['ResidueID'])
resids['chi1']=set(conf['ResidueID'][frames])


#Example 1: chi1 only
mymetric=metrics.Dihedral(metric='euclidean', p=2, angles='phi')
t=Trajectory.LoadFromLHDF('/home/shukla/GPCR_Mar27/b2ar3p0g2rh1_car/Trajectories/trj1.lh5')
phi_ptraj=mymetric.prepare_trajectory(t, degrees=True)

mymetric=metrics.Dihedral(metric='euclidean', p=2, angles='psi')
t=Trajectory.LoadFromLHDF('/home/shukla/GPCR_Mar27/b2ar3p0g2rh1_car/Trajectories/trj1.lh5')
psi_ptraj=mymetric.prepare_trajectory(t, degrees=True)

mymetric=metrics.Dihedral(metric='euclidean', p=2, angles='chi1')
t=Trajectory.LoadFromLHDF('/home/shukla/GPCR_Mar27/b2ar3p0g2rh1_car/Trajectories/trj1.lh5')
chi_ptraj=mymetric.prepare_trajectory(t, degrees=True)

mymetric=metrics.Dihedral(metric='euclidean', p=2, angles='phi psi chi1')
ptraj=mymetric.prepare_trajectory(t, degrees=True)


print chi_ptraj.shape
print phi_ptraj.shape
print psi_ptraj.shape
import pdb
pdb.set_trace()
