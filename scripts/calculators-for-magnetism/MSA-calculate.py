import numpy as np
from datetime import datetime
from ase import Atoms
from ase.io import read
from ase.neighborlist import NeighborList

def build_cell(element, position, cell, pbc):
    atom_object = Atoms(element, positions=position, cell=cell, pbc=pbc)
    return atom_object

def read_cell(file_path, element, pbc):
    atom_object = read(file_path)
    atom_object.pbc = pbc
    element_list = atom_object.get_chemical_symbols()
    indices = [index for index, value in enumerate(element_list) if value != element]
    del atom_object[indices]
    return atom_object

def get_neighborlist(atom_object, cutoff):
    atom_num = len(atom_object)
    cutoff = [cutoff]*atom_num
    nl = NeighborList(cutoff, self_interaction=False, bothways=True)
    nl.update(atom_object)
    atom_neigh = []
    for i in range(atom_num):
        indices, offsets = nl.get_neighbors(i)
        dxyz = atom_object.positions[indices] + np.dot(offsets, atom_object.cell) - atom_object.positions[i]
        dr = np.sqrt(np.sum(dxyz**2, axis=1))
        atomR = np.column_stack((indices,dxyz,dr))
        atomR[:, 0] = atomR[:, 0].astype(int)
        atom_neigh.append(atomR)
    return atom_neigh

def get_mm(atom_object, magnetic_moment):
    atom_num = len(atom_object)
    mmx = [[0, 0, 0]]*atom_num
    mmy = [[0, 0, 0]]*atom_num
    mmz = [[0, 0, 0]]*atom_num
    for i in range(atom_num):
        mmx[i][0] = magnetic_moment[i]
        mmy[i][1] = magnetic_moment[i]
        mmz[i][2] = magnetic_moment[i]
    return mmx, mmy, mmz

def calculate_MSA(atom_neigh, magnetic_moment=[[0,0,1],[0,0,2]]):
    dipole_lst = []
    
    for i in range(len(atom_neigh)):
        mmi = magnetic_moment[i]
        dipole = 0
        for j in range(len(atom_neigh[i])):
            indices = atom_neigh[i][:,0][j].astype(int)
            mmj = magnetic_moment[indices]
            dxyz = atom_neigh[i][:,1:4][j]
            dr = atom_neigh[i][:,4][j]

            temp = np.dot(mmi, mmj)-3*np.dot(mmi, dxyz)*np.dot(mmj, dxyz)/(dr**2)
            temp_A = temp*((9.274*1e-24) ** 2)
            constant = 1e-7
            J_dipole = constant*temp_A/(2*dr**3)*1e30
            miueV_dipole = J_dipole * 6.241506363 * 1e24

            dipole = dipole + miueV_dipole
        
        dipole_lst.append(dipole)
    return dipole_lst

    
if __name__ == '__main__':
    start_time = datetime.now()
    file_path = './POSCAR'
    pbc = [1, 1, 0]
    mag_element = 'Cr'
    mag_moment = [3, 3]
    rcut = 1000
    
    atom_object = read_cell(file_path, mag_element, pbc)
    atom_neigh = get_neighborlist(atom_object, rcut/2)
    mmx, mmy, mmz = get_mm(atom_object, magnetic_moment=mag_moment)
    dipole_z = calculate_MSA(atom_neigh, mmz)
    dipole_y = calculate_MSA(atom_neigh, mmy)
    dipole_x = calculate_MSA(atom_neigh, mmx)
    for i in range(len(dipole_x)):
        print('第{0}个原子:\n x方向: {1:10f} μeV, y方向: {2:.10f} μeV, z方向: {3:.10f} μeV'.format(i+1, dipole_x[i], dipole_y[i], dipole_z[i]))
        print(' out of plane - in plane: {0:.10f} μeV'.format(dipole_z[i]-min(dipole_x[i], dipole_y[i])))
    end_time = datetime.now()
    print(end_time-start_time)
    
