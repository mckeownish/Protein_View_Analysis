'''
------------------------------------------------------------------------------------------------
Weclome to Carbonara's Data Tools (CDT)!
------------------------------------------------------------------------------------------------

This package provides processing for Carbonara specific data input/outputs 
and generating useful inputs for external programs

------------------------------------------------------------------------------------------------
Package functionalities:
------------------------------------------------------------------------------------------------

- Coordinate, secondary structure and sequence extraction from xyz.dat and fingerprint.dat files

>>> extract_coords(coords_file)
>>> extract_sequence(fingerprint_file)

------------------------------------------------------------------------------------------------

- Generation of CA only PDB from xyz.dat file

>>> Carbonara_2_PDB(coords_file, fingerprint_file, output_file)

------------------------------------------------------------------------------------------------

- Finding geometric CB positions from CA only

>>> infer_CB_positions(CA_xyz)

------------------------------------------------------------------------------------------------

- Generation of CA + CB PDB (ready for reconstruction constraints)

>>> interlace_CA_CB_write(CA_xyz, CB_xyz, protein, output_pdb)
>>> CA_PDB_2_CB_PDB(CA_pdb, output_pdb)

------------------------------------------------------------------------------------------------

- (coming soon!) Reconstruction of full atomistic structure using MODELLER

>>> {CA_2_AA}
>>> {CACB_2_AA}

------------------------------------------------------------------------------------------------

- DSSP prediction of residue secondary structure

>>> 

------------------------------------------------------------------------------------------------

- (coming soon!) Identification of Varying Sections 


------------------------------------------------------------------------------------------------
'''

import pandas as pd
import numpy as np

from scipy.spatial.distance import cdist

import os
import subprocess
import shutil
from tqdm import tqdm

from Bio.PDB import PDBParser
import biobox as bb
from Bio.PDB.DSSP import DSSP

import shutil

import plotly.graph_objects as go



def extract_coords(coords_file):
    
    coords = np.genfromtxt(coords_file)
    coords = coords[~np.isnan(coords).any(axis=1)]
    
    return coords


def extract_sequence(fingerprint_file):
    
    seq = open(fingerprint_file, 'r').readlines()[2][:-1]
    
    return seq


def Carbonara_2_PDB(coords_file, fingerprint_file, output_file):
    
    '''
    Writes alpha carbon PDBs from Carbonara output
    
    Input
        coords_file      : coordinates of the carbon alpha chain
        fingerprint_file : Carbonara specific format containing secondary structure and sequence
        output_file      : define name of write output
    '''
    
    # read in coordinates and fingerprint 
    coords = extract_coords(coords_file)
    size = coords.shape[0]
    seq = extract_sequence(fingerprint_file)
    
    # map the AA shorthand to 3 letter abr.     
    aa_map = {
            'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
            'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
            'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
            'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
            'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
        }
    
    # map sequence to 3 character abr.
    seq_3 = []
    for a in list(seq):
        seq_3.append(aa_map[a])
    
    # create dataframe for biobox molecule type
    df = pd.DataFrame({'atom':['ATOM']*size, 'index':np.arange(size), 'name':['CA']*size, 
                       'resname':seq_3, 'chain':['A']*size, 'resid':np.arange(size),
                       'occupancy':[1]*size, 'beta':[50]*size, 'atomtype':['C']*size, 
                       'radius':[1.7]*size, 'charge':[0]*size})
    
    # take full advantage of Matteo's lovely biobox library - manually 'create' a molecule 
    molecule = bb.Molecule()
    molecule.data = df
    molecule.coordinates = np.expand_dims(coords, axis=0)
    
    # write out!
    molecule.write_pdb(output_file)
    

def infer_CB_positions(CA_xyz):
    
    '''
    Returns CB positions (excluding end of chains)
    '''
    
    CA_vecs = np.diff(CA_xyz, axis=0)
    normals = np.diff(CA_vecs, axis=0)
    normals = normals/np.linalg.norm(normals, axis=1)[:,None]

    av_bond_len = 3.8

    normals = normals*av_bond_len

    CB_xyz = CA_xyz[1:-1] - normals # minus as to face outwards
    
    return CB_xyz


    

    
def interlace_CA_CB_write(CA_xyz, CB_xyz, protein, output_pdb):
    
    # positions of CB - wont add CB to the residues
    gly_idx = np.where(protein.data['resname'].values== 'GLY')[0]
  
    # atom names
    atom_name_lst = []

    # interlaced CA + CB coordinates
    coordinates = CA_xyz[0]

    # number of CA
    size = CA_xyz.shape[0]

    # creating index for resid
    idx_counter = 0
    idx_lst = []

    # new interlaced CA CB sequence
    new_seq = []

    # extract CA sequence
    ca_seq = protein.data['resname'].values 
    
    for i in range(size):

        # ...CA bit...
        idx_lst.append(idx_counter)
        atom_name_lst.append('CA')
        new_seq.append(ca_seq[i])

        if i > 0:
                coordinates = np.vstack([coordinates, CA_xyz[i]])

        # ...CB bit...
        if i not in gly_idx:

            if (i > 0) & (i < size-1):
                idx_lst.append(idx_counter)
                atom_name_lst.append('CB')
                coordinates = np.vstack([coordinates, CB_xyz[i-1]])
                new_seq.append(ca_seq[i])

        idx_counter += 1

        
    tot_size = int(CA_xyz.shape[0]+CB_xyz.shape[0]-gly_idx.shape[0])

    if coordinates.shape[0] == tot_size:

        # create dataframe for biobox molecule type
        df = pd.DataFrame({'atom':['ATOM']*tot_size, 'index':np.arange(tot_size),
                           'name':atom_name_lst, 'resname':new_seq,
                           'chain':['A']*tot_size, 'resid':idx_lst,
                           'occupancy':[1]*tot_size, 'beta':[50]*tot_size,
                           'atomtype':['C']*tot_size, 'radius':[1.7]*tot_size, 
                           'charge':[0]*tot_size})

    else:
        raise ValueError('Total number of CA + CB - no. GLY res does not equal the coordinate size!')
        
        
    molecule = bb.Molecule()
    molecule.data = df
    molecule.coordinates = np.expand_dims(coordinates, axis=0)
    
    molecule.write_pdb(output_pdb)
    
    
    
def CA_PDB_2_CB_PDB(CA_pdb, output_pdb):
    
    # Load protein into BioBox object
    protein = bb.Molecule(CA_pdb)
     
    # Get CA coordinates
    CA_xyz = protein.coordinates[0]
    
    # infer the CB positions
    CB_xyz = infer_CB_positions(CA_xyz)
    
    interlace_CA_CB_write(CA_xyz, CB_xyz, protein, output_pdb)
    
    
    
def DSSP_structure_extractor(pdb_file):
    
    '''
    Use DSSP for predict secondary structure from a PDB file
    
    returns a numpy array of residue secondary structure labels
    '''
    
    p = PDBParser()
    structure = p.get_structure("PDB_file", pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)

    simplify_dict = {'H': 'H', 'B': 'S', 'E': 'S', 'G': 'H', 'I': 'H', 'T': '-', 'S': '-', '-': '-', ' ': '-'}
    secondary_struct = []

    for key in list(dssp.keys()):
        secondary_struct.append( simplify_dict[ dssp[key][2] ] )
        
    return np.asarray(secondary_struct)    
        
    
    
    
def get_secondary(fingerprint_file):
    return list(np.loadtxt(fingerprint_file, str)[2])


def read_coords(coords_file):
    
    coords = np.genfromtxt(coords_file)
    coords = coords[~np.isnan(coords).any(axis=1)]
    
    return coords


def section_finder(ss):
    
    '''Find protein sub-unit sections from the full secondary structure'''
    
    sections = []
    structure_change = np.diff(np.unique(ss, return_inverse=True)[1])

    for i, c in enumerate( structure_change ):

        if c!=0:
            sections.append(ss[i])

        if i==structure_change.shape[0]-1:
            sections.append(ss[i])
            
    sections = np.array(sections)
    
    return sections #, linker_indices #, structure_change


def find_sheet_indices(sections):
    
    '''Find sheet sub-unit section indices'''

    sheet_indices = np.where(sections=='S')[0]
    return sheet_indices


def find_linker_indices(sections):
    
    '''Find linker sub-unit section indices'''
    
    linker_indices = np.where(sections=='-')[0]
    return linker_indices


def generate_random_structures(coords_file, fingerprint_file):
    
    '''Generate random structures changing one linker section at a time
    
    Parameters
    coords_file:       /path/ to CA coordinates.dat file
    fingerprint_file:  /path/ to fingerprint.dat file
    
    Return
    Generated structures are written to ~/rand_structures/.. section_*LINKERINDEX*.dat as xyz
    Linker Indices
    ''' 
    
    linker_indices = find_linker_indices( section_finder( get_secondary(fingerprint_file) ) ) 
    
    current = os.getcwd()
    random = 'rand_structures'
    random_working = os.path.join(current, random)

    if os.path.exists(random_working) and os.path.isdir(random_working):
        shutil.rmtree(random_working)

    os.mkdir(random_working)

    # try:
        
    # except OSError as error:
    #     print(str(error)[11:])

    # print('Beginning random structures generation \n')
    
    rand_file_dict = {}
    for l in tqdm(linker_indices):
        
        outputname = random_working+'/section_'+str(l)
        
#         !./generate_structure {fingerprint_file} {coords_file} {outputname} {l}
        result = subprocess.run(['./generate_structure', fingerprint_file, coords_file, outputname, str(l)], capture_output=True, text=True)

    # print('')
    # print('Finished generating random structures')
    
    return linker_indices


def sheet_group_mask(ss):
     
    '''Groups adjacent sheets in secondary structure file and returns a grouping mask ( 0 : not a sheet;  1+: sheet )
    
    Parameters
    ss (numpy array):            Secondary structure labels (array of strings)
    
    Returns
    sheet_groups (numpy array):  Mask of grouped sheet sections
    '''
    
    sheet_mask = (ss == 'S')*1
    sheet_groups = np.zeros(ss.shape[0])
    group = 1
    
    if sheet_mask[0] == 1:
        label = True
    else:
        label = False

    for i, c in enumerate(np.diff(sheet_mask)):
        
        
        if c == 1:
            label = True

        elif c==-1:
            label=False
            group += 1

        else:
            pass 

        if label == True:
            if ss[i+1] == 'S':
                sheet_groups[i+1] = group
                
    return sheet_groups


def linker_group_mask(ss):
    
    '''Groups adjacent linkers in secondary structure file and returns a grouping mask ( 0 : not a linker;  1+: linker )
    
    Parameters
    ss (numpy array):             Secondary structure labels (array of strings)
    
    Returns
    linker_groups (numpy array):  Mask of grouped linker sections
    '''
    
    linker_mask = (ss == '-')*1
    linker_groups = np.zeros(ss.shape[0])
    group = 1
    
    # checking first index for linker 
    if linker_mask[0] == 1:
        label = True
        linker_groups[0] = group
    else:
        label = False

    for i, c in enumerate(np.diff(linker_mask)):
    
        if c == 1:
            label = True

        elif c==-1:
            label=False
            group += 1

        else:
            pass 

        if label == True:
            
            linker_groups[i+1] = group
                
    return linker_groups #, linker_mask


def get_sheet_coords(coords, sheet_groups):

    '''Finds CA coordinates of 
    
    Parameters
    coords (numpy array):        xyz coordinates of all protein CA atoms
    sheet_groups (numpy array):  Mask of grouped sheet sections
    
    Returns
    sheet_coords (numpy array):  xyz coordinates of CA atoms in each sheet structure [ [...sheet 1 coords...] [...sheet 2 coords...] ... ]
    '''
    
    sheet_coords = []

    for g in np.unique(sheet_groups):
        if g>0:
            sheet_coords.append(coords[sheet_groups==g])
            
    sheet_coords = np.asarray(sheet_coords)
    
    return sheet_coords



def get_section_groupings(ss, structure_change):
    
    group = 0
    structural_groups = np.zeros(ss.shape)
    structural_groups[0] = group

    for i, c in enumerate(structure_change):

        if c != 0:
            group += 1

        structural_groups[i+1] = group
    return structural_groups


def list_nohidden(path):
    lst = []
    for f in os.listdir(path):
        if not f.startswith('.'):
            lst.append(f)
    return lst


def sheet_pipe(coords_file, fingerprint_file):
   
    coords = read_coords(coords_file)
    ss = get_secondary(fingerprint_file)
    sheet_groups = sheet_group_mask( np.asarray(ss) )
    sheet_coords = get_sheet_coords(coords, sheet_groups)

    return sheet_coords


def sheet_pairwise_bond_number(sheet_coords, thr=5.5):
    
    '''Finds the number of pairs of CA atoms within some threshold between all sheet sections
    
    Parameters
    sheet_coords (numpy array): xyz coordinates of CA atoms in each sheet structure [ [...sheet 1 coords...] [...sheet 2 coords...] ... ]
    thr (float) {optional}:     Cutoff distance for inter-sheet bonding (default = 5.5 Å)
    
    Returns
    pairwise_bond_num (numpy array): Lower triangular array containing the number of individual CA bonds within threshold between each sheet pair
    
    '''
    
    number_bonds = 0

    pairwise_bond_num = np.zeros([len(sheet_coords), len(sheet_coords)])

    for i in range(1,len(sheet_coords)):

        for j in range(0,i):

            arr1, arr2 = sheet_coords[j], sheet_coords[i]
            dist_matrix = cdist(arr1, arr2)
            indices = np.where(dist_matrix < thr)

            pairwise_bond_num[i,j] = indices[0].shape[0]

            number_bonds += indices[0].shape[0]
    return pairwise_bond_num 


def random_bond_finder(rand_file_dir, fingerprint_file, linker_indices):
   
    # grouping all random structure for each linker together
    
    struture_lst = list_nohidden(rand_file_dir)
   
    linker_file_dict = {}
    for l in linker_indices:
        tmp = []
       
        for file in np.sort(struture_lst):
            if str(l) == file.split('_')[1]:
                tmp.append(file)

        linker_file_dict[l] = tmp
       
    # Pairwise sheet bonds for each random str for each linker
    linker_bond_dict = {}

    for l in linker_indices:
       
        tmp = []
       
        for file in linker_file_dict[l]:
            coords_file = rand_file_dir+file
            sheet_coords = sheet_pipe(coords_file, fingerprint_file)
            tmp.append( sheet_pairwise_bond_number(sheet_coords) )
   
        linker_bond_dict[l] = tmp
    
    return linker_bond_dict 


def find_non_varying_linkers():
    
    initial_coords_file = 'Fitting/coordinates1.dat'
    fingerprint_file = 'Fitting/fingerPrint1.dat'
    
    # Reference initial structure
    sheet_coords = sheet_pipe(initial_coords_file,
                              fingerprint_file)
    ref_bonds = sheet_pairwise_bond_number(sheet_coords, thr=5.5)

    # Generate the random structure changing each linker section
    linker_indices = generate_random_structures(initial_coords_file,
                                                        fingerprint_file)
    
    # Calculate the number of inter-sheet bonds for each rand struct
    linker_bond_arr_dict = random_bond_finder('rand_structures/', 
                                              fingerprint_file,
                                              linker_indices)
    
    # Find number of bond breaks relative to initial structure
    bond_breaks_dict = {}

    for l in linker_indices:

        bond_break_lst = []
        for bond_arr in linker_bond_arr_dict[l]:


            bond_break_lst.append( (ref_bonds > bond_arr).sum() )

        bond_breaks_dict[l] = sum(bond_break_lst)/(len(linker_bond_arr_dict[l])+1)
        
    
    # Linker indices that cause no bond breaks
    conds = np.asarray(list(bond_breaks_dict.values())) < 0.001
      
    allowed_linker = linker_indices[conds]

    if 0 in linker_indices:
        linker_indices = np.delete(linker_indices, np.where(linker_indices==0)[0].item())


    if 0 in allowed_linker:
        allowed_linker = np.delete(allowed_linker, np.where(allowed_linker==0)[0].item())

    return allowed_linker, linker_indices


# ------ Carbonara Setup Funcs ---------


def setup_working_directory():
    
    current = os.getcwd()
    working = 'Fitting'
    working_path = os.path.join(current, working)
    
    if os.path.exists(working_path):
        shutil.rmtree(working_path)
        print('Removing existing working directory')
        
    os.makedirs(working_path)
    os.mkdir(working_path+'/fitdata')

    print('Complete')
    return working_path


def pdb_2_biobox(pdb_file):
    M = bb.Molecule()
    M.import_pdb(pdb_file)
    return M


def extract_CA_coordinates(M):
    ca_idx = (M.data['name']=='CA').values
    ca_coords = M.coordinates[0][ca_idx]
    
    if ca_coords.shape[0] != M.data['resid'].nunique():
        raise Exception("You better check your PDB... The number of CA atoms does not equal the number of ResIDs in your PDB file!") 
    else:
        return ca_coords

    
def extract_sequence(M):
    
    
    aa_names = {
                'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
                'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
                'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
                'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
                }

    names_aa = {y: x for x, y in aa_names.items()}
    
    ca_idx = (M.data['name']=='CA').values
    resnames = M.data['resname'][ca_idx].map(names_aa).values
    
    if resnames.shape[0] != M.data['resid'].nunique():
        raise Exception("You better check your PDB... The number of CA atoms does not equal the number of ResIDs in your PDB file!") 
    else:
        return resnames


def write_fingerprint_file(number_chains, sequence, secondary_structure, working_path):
    
    assert isinstance(number_chains, int), 'Yikes... The number of chains is not int type!'
    
    if number_chains > 1:
        print('Are sure you have more than one chain - if not this will cause segmentation errors later! You have been warned...')
    
    seq_run = ''.join(list(sequence))
    ss_run = ''.join(list(secondary_structure))
    
    if len(seq_run) != len(ss_run):
        raise Exception("Uh Oh... The length of sequence and secondary structure is not equal!") 
    
    f = open(working_path+"/fingerPrint1.dat", "w")
    f.write(str(number_chains))
    f.write('\n \n')
    f.write(seq_run)
    f.write('\n \n')
    f.write(ss_run)
    f.close()
    
    
def write_coordinates_file(coords, working_path):
    
    assert type(coords).__module__ == np.__name__, 'Thats never good... the CA coordinates are not a numpy array'
    np.savetxt(working_path+'/coordinates1.dat', coords, delimiter=' ', fmt='%s',newline='\n', header='', footer='')
    
    
def write_mixture_file(working_path):
    # if default:
    f = open(working_path+"/mixtureFile.dat", "w")
    f.write(str(1))
        
#     else:
#          copy input file


def write_varysections_file(varying_sections, working_path):
    # auto: run beta sheet breaking code; write output sections to file
    f = open(working_path+"/varyingSectionSecondary1.dat", "w")
    for i, s in enumerate(varying_sections):
        f.write(str(s))
        
        if i < len(varying_sections)-1:
            f.write('\n')
    f.close()

    
def write_saxs(SAXS_file, working_path):
    
    saxs_arr = np.genfromtxt(SAXS_file)
    
    if saxs_arr.shape[1] == 3:
        saxs_arr = saxs_arr[:,:2]
        
    np.savetxt(working_path+'/Saxs.dat', saxs_arr, delimiter=' ', fmt='%s',newline='\n', header='', footer='')


def read_dssp_file(dssp_filename):
    
    simplify_dict = {'H': 'H', 'B': 'S', 'E': 'S', 'G': 'H', 'I': 'H', 'T': '-', 'S': '-', '-': '-', ' ': '-'}
    
    lines=[]
    with open(dssp_filename) as input_data:
        # Skips text before the beginning of the interesting block:
        for line in input_data:
            if line.strip() == '#  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA': 
                break
        # Reads text until the end of the block:
        for line in input_data:  # This keeps reading the file
            lines.append(simplify_dict[line[16]])
    return ''.join(lines)
    
    
def simplify_secondary(dssp_struct):
    
    simplify_dict = {'H': 'H', 'B': 'S', 'E': 'S', 'G': 'H', 'I': 'H', 'T': '-', 'S': '-', '-': '-', ' ': '-'}
    
    secondary_structure = []
    
    for s in dssp_struct:
        
        if s not in list(simplify_dict.keys()):
            print('>>> ', s, ' <<<')
            raise Exception('Secondary structure not recognised!')
            
        secondary_structure.append(simplify_dict[s])
        
    return secondary_structure


def write_sh_file(working_path, fit_n_times, min_q, max_q, max_fit_steps):
    
    curr = os.getcwd()
    run_file = curr + '/RunMe.sh'

    with open(run_file, 'w+') as fout:
        fout.write('#!/bin/bash')
        
        saxs_file = working_path+'/Saxs.dat'
        FP_file = working_path+"/fingerPrint1.dat"
        coords_file = working_path+'/coordinates1.dat'
        varying_file = working_path+"/varyingSectionSecondary1.dat"
        mixture_file = working_path+"/mixtureFile.dat"
        
        # Auto assign min / max q from SAXS profile
        # saxs_arr = np.genfromtxt(saxs_file)
        # min_q = np.round(saxs_arr[:,0].min(),2)
        # max_q = np.round(saxs_arr[:,0].max(),2)
        
        fout.write('\nfor i in {1..'+str(fit_n_times)+'}')

        fout.write('\n\ndo')
        fout.write('\n\n   echo " Run number : $i "')
        fout.write('\n\n   ./predictStructure ' + saxs_file + ' ' + working_path+'/' + ' ' + coords_file + ' ' + 'none' + ' ' + varying_file + ' ' + '1' + ' ' + 'none' + \
                   ' ' + 'none' + ' ' + str(min_q) + ' ' + str(max_q) + ' ' + str(max_fit_steps) + ' ' + working_path+'/fitdata/fitmolecule$i' + ' ' + working_path+'/fitdata/scatter$i.dat' + ' ' + mixture_file + ' ' +'1')
                   
        fout.write('\n\ndone')
        
    print('Successfully written bash script to: ', run_file) 




def SAXS_selection_plotter(SAXS_file, min_q, max_q):

    SAXS = np.genfromtxt(SAXS_file)

    q = SAXS[:,0]
    I = SAXS[:,1]

    q_selected = q[(q>=min_q)&(q<=max_q)]
    q_grey = q[(q<min_q) | (q>=max_q)]

    I_selected = I[(q>=min_q)&(q<=max_q)]
    I_grey = I[(q<min_q) | (q>=max_q)]

    fig = go.Figure( data=[go.Scatter(x=q_grey, y=I_grey, mode='markers', line=dict(color="grey"), opacity=0.7, name='Excluded')])
    fig.add_trace( go.Scatter(x=q_selected, y=I_selected, mode='markers', line=dict(color="crimson"), name='Selected') )
    fig.update_layout(title='Selected q range', yaxis_type = "log", template='plotly_white',
                    width=800, height=700, font_size=28)
    fig.update_xaxes(title='q')
    fig.update_yaxes(title='I')

    return fig


def get_minmax_q(SAXS_file):
    
    SAXS = np.genfromtxt(SAXS_file)

    q = SAXS[:,0]

    q_exp_min = float(np.round(q.min(),2))
    q_exp_max = float(np.round(q.max(),2))

    q_spread = q_exp_max - q_exp_min

    q_Q1 = float(np.round(0.00*q_spread,2))
    q_Q3 = float(np.round(0.45*q_spread,2))

    return q_exp_min, q_exp_max, q_Q1, q_Q3


def sort_by_creation(file_lst):

    files = list(filter(os.path.isfile, file_lst))
    files.sort(key=lambda x: os.path.getmtime(x))
    return files