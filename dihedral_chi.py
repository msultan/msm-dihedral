'''
Compute dihedral angles for each frame in a trajectory

Implemented in C behind the scenes, so it's pretty fast I think
'''


from Emsmbuilder import _dihedral_wrap
import numpy as np

def get_indices(trajectory_or_conformation, angles='phi/psi', residues=None):
    '''Get the atom indices of the quartets of atoms involved in
    the dihedral angle interactions
    
    trajectory_or_conformation can be a Conformation or a Trajectory object
    
    angles can be a single string or a list of strings, that can be parsed
    as any combination of phi, psi and or chi
        ['phi','psi']
        'phi/psi'
        'phi psi'
    are all valid.
    
    Returns a N x 4 array
    '''
    
    if isinstance(angles, basestring):
        angles = angles.lower()
        angles = angles.replace('/', ' ').replace('-', ' ').split()
    try:
        angles = [angle.lower() for angle in angles]
    except:
        raise ValueError("I can't parse %s. Please supply a string like phi/psi" % str(angles))
    indices = np.zeros((0,4), dtype=int)  
    if residues is None:
        print "Using the following residues for dihedral metric: all" 
    else:
        print "Using the following residues for dihedral metric: ", residues
      
    for angle in angles:
        if angle == 'phi':
            indices = np.vstack((indices, _get_indices_phi(trajectory_or_conformation, residues)))
        elif angle == 'psi':
            indices = np.vstack((indices, _get_indices_psi(trajectory_or_conformation, residues)))
        elif angle == 'chi1':
            indices = np.vstack((indices, _get_indices_chi1(trajectory_or_conformation, residues)))
        else:
            raise ValueError("Uncregozied angle type: %s. Only phi, psi and chi are supported" % angle)
    return indices


def compute_dihedrals(trajectory_or_conformation, indices, degrees=True):
    '''Compute the dihedral angles in a XYZList
    
    xyzlist should be a 3D array (num_frames x num_atoms x num_dims) with
    the coordinates of the atoms.
    
    indices should be a 2D array (num_quartets x 4) with each row containing the
    indices of atoms you want to extract the dihedral angles for. This is the
    form returned by get_indices()
    
    degrees (boolean) controls whether the results are returned in degrees (if
    true) or radians (if false)
    
    Returns: a num_frames x num_quartets array of dihedral angles
    '''
    try:
        xyzlist = trajectory_or_conformation['XYZList']
    except KeyError:
        xyz = trajectory_or_conformation['XYZ']
        num_atoms, num_dims = xyz.shape
        xyzlist = np.reshape(xyz, (1, num_atoms, num_dims), order='C')
    
    num_frames, num_atoms, num_dims = xyzlist.shape
    if num_dims != 3:
        raise ValueError('num_dims must be three')
    num_quartets, four = indices.shape
    if four != 4:
        raise ValueError('indices must have 4 columns')
    
    indices = np.array(indices, dtype=np.int, order='C')
    
    if xyzlist.dtype == np.float32:
        xyzlist = np.array(xyzlist, dtype=np.float32, order='C')
        results = np.zeros((num_frames, num_quartets), dtype=np.float32)
        _dihedral_wrap.dihedrals_from_traj_float_wrap(results, xyzlist, indices)
    elif xyzlist.dtype == np.double:
        xyzlist = np.array(xyzlist, dtype=np.double, order='C')
        results = np.zeros((num_frames, num_quartets), dtype=np.double)
        _dihedral_wrap.dihedrals_from_traj_wrap(results, xyzlist, indices)
    else:
        raise ValueError('Unsupported type of xyzlist: %s' % xyzlist.dtype)
        
    if degrees:
        results *= 180.0 / np.pi
    return results


def _get_indices_phi(conformation, res):
    '''Get the atom indices of the quartets of atoms involved in
    each of the phi dihedral angles
    
    conformation can be a Conformation or a Trajectory object
    
    Returns a num_residues x 4 array
    '''
    
    NResi = conformation.GetNumberOfResidues()
    AID = conformation.GetEnumeratedAtomID()
    RID = conformation.GetEnumeratedResidueID()
    AName = conformation["AtomNames"]
    Indices = []
    if res==None:
        res=range(NResi)
    for i in res:
        try:
            a0 = np.where((AName == "C") & (RID == i))[0][0]
            a1 = np.where((AName == "N") & (RID == (i + 1)))[0][0]
            a2 = np.where((AName == "CA") & (RID == (i + 1)))[0][0]
            a3 = np.where((AName == "C") & (RID == (i + 1)))[0][0]
        except:
            pass
        Indices.append([a0, a1, a2, a3])
    return(np.array(Indices))


def _get_indices_psi(conformation, res):
    '''Get the atom indices of the quartets of atoms involved in
    each of the psi dihedral angles
    
    conformation can be a Conformation or a Trajectory object
    
    Returns a num_residues x 4 array
    '''
    
    NResi = conformation.GetNumberOfResidues()
    AID = conformation.GetEnumeratedAtomID()
    RID = conformation.GetEnumeratedResidueID()
    AName = conformation["AtomNames"]
    Indices = []
    if res==None:
        res=range(NResi)
    for i in res:
        try:
            a0 = np.where((AName == "N") & (RID == i))[0][0]
            a1 = np.where((AName == "CA") & (RID == i))[0][0]
            a2 = np.where((AName == "C") & (RID == i))[0][0]
            a3 = np.where((AName == "N") & (RID == (i + 1)))[0][0]
        except:
            pass
        Indices.append([a0, a1, a2, a3])
        
    return(np.array(Indices))

def _get_indices_chi1(conformation, res):
    '''Get the atom indices of the quartets of atoms involved in
    each of the chi dihedral angles
    
    conformation can be a Conformation or a Trajectory object
    MORGAN: have to make changes to get chi for residues:
        SER: CB=OG, CG=HG
        THR:  CG=OG1
        CYX: CG is SG
        GLY: none, ALA: none
    
    Returns a num_residues x 4 array
    '''
    
    NResi = conformation.GetNumberOfResidues()
    AID = conformation.GetEnumeratedAtomID()
    RID = conformation.GetEnumeratedResidueID()
    AName = conformation["AtomNames"]
    RName = conformation["ResidueNames"]
    Indices = []
    #modified 
    if res==None:
        res=range(NResi)
    for i in res:
        indices=np.where(RID==i)[0]
        if len(indices)<1:
            print "warning: residues %s does not exist" % i
            continue
        name=RName[indices]
        #skip capping fragments (have CA, N in them)
        if (name[0]=='ACE') or (name[0]=='NME') or (name[0]=='NMA'):
            print "Warning: skipping residue %s %s" % (i, name)
            continue
        testres=np.where((AName=="CA")&(RID==i))[0]
        if len(testres)==0: #skip non protein residues
            print "Warning: skipping residue %s %s" % (i, name)
            continue
        a0 = np.where((AName=="N")&(RID==i))[0][0]
        a1 = np.where((AName=="CA")&(RID==i))[0][0]
        if np.where((AName=="CB")&(RID==i))[0]:
            a2 = np.where((AName=="CB")&(RID==i))[0][0]
        #elif np.where((AName=="OG")&(RID==i)&(RName=='SER'))[0]:
        #    a2 = np.where((AName=="OG")&(RID==i)&(RName=='SER'))[0][0]
        else:
            a2 = None 
            print "Warning: skipping residue %s %s" % (i, name)
        if np.where((AName=="CG")&(RID==i))[0]:
            a3 = np.where((AName=="CG")&(RID==i))[0][0]
        elif np.where((AName=="CG1")&(RID==i))[0]:
            a3 = np.where((AName=="CG1")&(RID==i))[0][0]
        elif np.where((AName=="OG")&(RID==i)&(RName=='SER'))[0]:
            a3 = np.where((AName=="OG")&(RID==i)&(RName=='SER'))[0][0]
        elif np.where((AName=="OG1")&(RID==i)&(RName=='THR'))[0]:
            a3 = np.where((AName=="OG1")&(RID==i)&(RName=='THR'))[0][0]
        elif np.where((AName=="SG")&(RID==i)&(RName=='CYX'))[0]:
            a3 = np.where((AName=="SG")&(RID==i)&(RName=='CYX'))[0][0]
        elif np.where((AName=="SG")&(RID==i)&(RName=='CYS'))[0]:
            a3 = np.where((AName=="SG")&(RID==i)&(RName=='CYS'))[0][0]
        else:
            a3 = None
            print "Warning: skipping residue %s %s" % (i, name)
        if a3 == None or a2 == None:
            continue
        else:
            Indices.append([a0, a1, a2, a3])
    return(np.array(Indices))

