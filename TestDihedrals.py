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
mymetric=metrics.Dihedral(metric='euclidean', p=2, angles='chi1')
t=Trajectory.LoadFromLHDF('/home/shukla/GPCR_Mar27/b2ar3p0g2rh1_car/Trajectories/trj1.lh5')
size=t['XYZList'].shape[0]
residues=numpy.zeros((len(resids['chi1']), size))
#t.SaveToXTC('trj1.xtc')
ptraj=mymetric.prepare_trajectory(t, degrees=True)
dihedral_count=0
for (n, id) in enumerate(resids['chi1']):
    frames=numpy.where(t['ResidueID']==id)[0]
    name=t['ResidueNames'][frames[0]]
    numpy.savetxt('ptraj-files/chi1%s%i.dat' % (name, id), ptraj[:,n])
    init=numpy.where(residues[n]==0)[0][0] # this more useful when looping over trajectories, to extend the matrix 
    final=init+size
    residues[n][init:final]=ptraj[:,n]
numpy.savetxt('chi1_residue_matrix.txt', residues)

#Example 2: getting phi, psi, chi1 or combo
mymetric=metrics.Dihedral(metric='euclidean', p=2, angles='phi psi chi1')
ptraj=mymetric.prepare_trajectory(t, degrees=True)
n=0
for angle in ['phi', 'psi', 'chi1']:
    m=0
    print "on angle %s" % angle
    residues=numpy.zeros((len(resids[angle]), size))
    for id in resids[angle]:
        frames=numpy.where(t['ResidueID']==id)[0]
        name=t['ResidueNames'][frames[0]]
        numpy.savetxt('ptraj-combo-files/%s%s%i.dat' % (angle, name, id), ptraj[:,n])
        init=numpy.where(residues[m]==0)[0][0] # this more useful when looping over trajectories 
        final=init+size
        residues[m][init:final]=ptraj[:,n] #
        m+=1
        n+=1
    numpy.savetxt('%s_residue_matrix.txt' % angle,  residues)
