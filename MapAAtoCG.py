#!/usr/bin/python3

import numpy as np
import MDAnalysis as md
import subprocess

gmx_loc= '/compilations/gromacs/gromacs-2021.3/bin/'
GRO = '/data/manel/cholesterol/aa/DPPC-CHOL/70-30/confout.gro'
XTC = '/data/manel/cholesterol/aa/DPPC-CHOL/70-30/traj_comp.xtc'
TPR = '/data/manel/cholesterol/aa/DPPC-CHOL/70-30/topol.tpr'
u = md.Universe(GRO)

lip_resnames            = []
lip_bead_assignments    = []
lip_bead_names          = []

### DPPC mapping
resname = 'DPPC'
NC3  = ['N', 'C13', 'H13A', 'H13B', 'H13C', 'C14', 'H14A', 'H14B', 'H14C', 'C15', 'H15A',
        'H15B', 'H15C', 'C12', 'H12A', 'H12B', 'C11', 'H11A', 'H11B']      
PO4  = ['P', 'O13', 'O14', 'O11', 'O12']
GL1  = ['C2', 'HS', 'O21', 'C21', 'O22', 'C1', 'HA', 'HB' ]
GL2  = ['C3', 'HX', 'HY', 'O31', 'O32', 'C31']
C1A  = ['C22', 'H2R', 'H2S', 'C23', 'H3R', 'H3S', 'C24', 'H4R', 'H4S', 'C25', 'H5R', 'H5S']
C2A  = ['C26', 'H6R', 'H6S', 'C27', 'H7R', 'H7S', 'C28', 'H8R', 'H8S', 'C29', 'H9R', 'H9S']
C3A  = ['C210', 'H10R', 'H10S', 'C211', 'H11R', 'H11S', 'C212', 'H12R', 'H12S', 'C213', 'H13R', 'H13S']
C4A  = ['C214', 'H14R', 'H14S', 'C215', 'H15R', 'H15S', 'C216', 'H16R', 'H16S']
C1B  = ['C32', 'H2Y', 'H2X', 'C33', 'H3Y', 'H3X', 'C34', 'H4Y', 'H4X', 'C35', 'H5Y', 'H5X']
C2B  = ['C36', 'H6Y', 'H6X', 'C37', 'H7Y', 'H7X', 'C38', 'H8Y', 'H8X', 'C39', 'H9Y', 'H9X']
C3B  = ['C310', 'H10Y', 'H10X', 'C311', 'H11Y', 'H11X', 'C312', 'H12Y', 'H12X', 'C313', 'H13Y', 'H13X']
C4B  = ['C314', 'H14Y', 'H14X', 'C315', 'H15Y', 'H15X', 'C316', 'H16Y', 'H16X']
bead_assignments = [NC3, PO4, GL1, GL2, C1A, C2A, C3A, C4A, C1B, C2B, C3B, C4B]
bead_names       = ['NC3', 'PO4', 'GL1', 'GL2', 'C1A', 'C2A', 'C3A', 'C4A',
                    'C1B', 'C2B', 'C3B', 'C4B']
lip_resnames.append(resname)
lip_bead_assignments.append(bead_assignments)
lip_bead_names.append(bead_names)

### Cholesterol mapping
resname = 'CHL1'
ROH  = ['C2', 'C2', 'C4', 'C4', 'O3', 'O3', 'O3', 'O3', 'O3', 'O3']   
R1   = ['C1']                                     
R2   = ['C6', 'C6', 'C6', 'C6', 'C5', 'C7']             
R3   = ['C11', 'C9']                              
R4   = ['C15', 'C15', 'C15', 'C15', 'C15', 'C14','C14','C16','C16']          
R5   = ['C12']                                            
R6   = ['C19', 'C19', 'C19', 'C19', 'C19', 'C19', 'C19', 'C19', 'C19', 'C10']          
R7   = ['C18', 'C18', 'C18', 'C18', 'C18', 'C18', 'C18', 'C18', 'C13']           
C1   = ['C20']
C2   = ['C24']
bead_assignments = [ROH, R1, R2, R3, R4, R5, R6, R7, C1, C2]
bead_names       = ['ROH', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'C1', 'C2']
bead_sizes       = ['S',    'T',  'S',  'S',  'S',  'T',  'T',  'T',  'S',  'R']
lip_resnames.append(resname)
lip_bead_assignments.append(bead_assignments)
lip_bead_names.append(bead_names)

bead_agg = []
for i, molecule in enumerate(lip_resnames):
    targets = u.select_atoms('resname {}'.format(molecule)).residues
    for j in range(len(targets)):
        target = targets[j].atoms
        for bead in lip_bead_assignments[i]:
            agg_whole = u.select_atoms('name empty')
            for atom in bead:
                agg = target.select_atoms(f'name {atom}')
                agg_whole += agg
            bead_agg.append(agg_whole)
with md.selections.gromacs.SelectionWriter('index.ndx', mode='w') as ndx:
    for idx, agg in enumerate(bead_agg):
        ndx.write(agg)
        
p = subprocess.Popen("{}gmx trjconv -f {} -s {} -pbc mol -o pbc.xtc".format(gmx_loc, XTC, TPR)
                    , stdin=subprocess.PIPE, shell=True, universal_newlines=True)
p.communicate('System')
p.wait()

for atom in u.atoms:
    atom.name = 'H'
u.atoms.write('all_hydrogen.gro')

p = subprocess.Popen("{}gmx traj -f {} -s all_hydrogen.gro -n index.ndx -com -ng {} -oxt cg.gro".format(gmx_loc, GRO, len(bead_agg))
                     , stdin=subprocess.PIPE, shell=True, universal_newlines=True)
p.communicate('\n'.join(map(str,list(range(len(bead_agg))))))
p.wait()

p = subprocess.Popen("{}gmx traj -f pbc.xtc -s all_hydrogen.gro -n index.ndx -com -ng {} -oxt cg.xtc".format(gmx_loc, len(bead_agg))
                     , stdin=subprocess.PIPE, shell=True, universal_newlines=True)
p.communicate('\n'.join(map(str,list(range(len(bead_agg))))))
p.wait()


### Fix names in GRO for easier parameterization
names = []
for i, molecule in enumerate(lip_resnames):
    targets = u.select_atoms('resname {}'.format(molecule)).residues
    gro_names    = lip_bead_names[i]*len(targets)
    names += gro_names
    
with open('cg.gro') as cggro_in:
    gro_file = cggro_in.readlines()
for idx, line in enumerate(gro_file[2:-1]):
    gro_file[2+idx]=line.replace(line[12:15],names[idx])
with open('cg_good.gro', 'w+') as file_out:
    for line in gro_file:
        file_out.write(line)