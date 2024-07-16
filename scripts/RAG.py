#!/usr/bin/env python
"""
  Young-Jae Choi, POSTECH, Korea, Rep. of.
  Inquiries: ssrokyz@gmail.com
"""

import numpy as np

def count2list(dict_in):
    list_out = []
    for i in range(len(dict_in.keys())):
        key = list(dict_in.keys())[i]
        num = 0
        itera = dict_in[key]
        for j in range(itera):
            list_out.append(key)
            num += 1
    return list_out

def list2count(list_inp):
    """
    input: list.
    output: dict of counts of elements.
    """
    keys = list(set(list_inp))
    dict_out = dict()
    for i in keys:
        dict_out[i] = 0
    for i in list_inp:
        dict_out[i] += 1
    return dict_out

def covalent_expect(input):
    """ Returns covalent bond expectation value of system 
    input : dict or list
    e.g. {'Si':1, 'H':2} or ['Si', 'H', 'H']
    """

    from ase.atoms import symbols2numbers as s2n
    from ase.data import covalent_radii
    
    if isinstance(input, list):
        num_spec_dict = list2count(input)
    elif isinstance(input, dict):
        num_spec_dict = input
    else:
        raise TypeError("input is not list nor dict")
    
    tot_num = np.sum(list(num_spec_dict.values()))
    r_sum = 0
    for key, value in num_spec_dict.items():
        r_sum += covalent_radii[s2n([key])[0]] * value
    expect_value = r_sum / tot_num
    return expect_value

def random_atoms_gen(
    backbone,
    num_spec_dict = None,
    fix_ind_dict  = None,
    pin_the_fixed = False,
    cutoff_radi   = None,
    cutoff_frac   = None,
    random_radi   = None,
    random_frac   = None,
    strain        = None,
    strain_ratio  = [1.,1.,1.],
    vacuum        = None,
    vacuum_ratio  = None,
    sort_atoms    = False,
    max_trial_sec = 5,
    log           = True,
    ):
    """ 
    RAG generates randomly positioned image based on backbone structure.

    backbone : An ASE atoms object
        All atomic positions of the backbone structure will be the lattice sites.

    num_spec_dict : dict or None
        Number of atoms for each species. If None, set to be identical to the backbone structure.
        "V" correspond to vacancy (Not placing any atom at the lattice site).
        Note that the following condition must be satisfied.
            --> np.sum(num_spec_dict.values()) == len(backbone)
        E.g.) {'Ge': 3, 'Te': 4, 'V': 1}
           --> Condition: the backbone structure must have 8 atoms, totally.

    fix_ind_dict : dict or list or None
        Dict of atomic indices in the backbone for each species whose sites will not be shuffled.
        It can be a list. The site will not be shuffled for the atoms with indices included in the list.
        Note that every fixed atom will also have positional deviation, unless 'pin_the_fixed' set to be 'True'.
        Set to "None" if all atoms' positions must be shuffled.
        E.g.) {'Te': [0, 2, 4, 6]}
            * I.e. The lattice sites at the positions, backbone.positions[[0,2,4,6]], will be filled by Te atoms, and others will be shuffled.

    pin_the_fixed : Boolean
        If true, atoms included in the fix_ind_dict will not be deviated from their lattice sites (even if the cutoff_radi and cutoff_frac are not zero).

    Note)
        If both cutoff_radi and cutoff_frac are provided simultaneously, occurs error.
        If both cutoff_radi and cutoff_frac are not provided, cutoff_radi is set to be zero by default. i.e. no constraint

        cutoff_radi : Float
            Cutoff radius for the minimum distance between two atoms.
            If there is a pair of atoms closer than the cutoff_radi, the position of the atom will be randomized again.

        cutoff_frac : Float
            Set cutoff_radi as {cutoff_frac *2 *(expected value of covalent radius of the material)}.
            i.e. When cutoff_frac == 1, all atoms are spaced with at least two times of expected covalent radius.

    Note)
        If both random_radi and random_frac are provided simultaneously, occurs error.
        If both random_radi and random_frac are not provided, random_radi is set to be zero by default. i.e. no deviation

        random_radi : float
            Maximum magnitude of the random deviation vector.
            Each atom will deviated from its lattice point by a random vector shorter than random_radi.

        random_frac : Float
            Set random_radi as {random_frac *0.5 *(position of RDF first peak)}
            RDF = radial distribution function or pair correlation function
            !!Caution!!
                When this option is switched on, RDF calculation will be carried out, which is a time consuming process for big systems.

    strain : List of three floats e.g. [0.,0.,5.]
        Values specify how much you magnify the provided backbone cell.
        Cell gets longer along lattice vectors.
        positions will stay scaled positions of backbone atoms.

    strain_ratio : List of three floats e.g. [1.,1.1,1.]
        Values specify how much you magnify the provided backbone cell.
        Cell gets longer along lattice vectors.
        positions will stay scaled positions of backbone atoms.

    vacuum : List of three floats e.g. [0.,0.,5.]
        Values specify how much you magnify the provided backbone cell with vacuum.
        Cell gets longer along lattice vectors.
        positions will stay absolute positions of backbone atoms.
        insert vacuum after strain (if provided)

    vacuum_ratio : List of three floats e.g. [1,1.1,1]
        Values specify how much you magnify the provided backbone cell with vacuum.
        Cell gets longer along lattice vectors.
        positions will stay absolute positions of backbone atoms.
        insert vacuum after strain (if provided)

    sort_atoms : True or False
        Sort atomic indices with increasing order in atomic numbers.

    max_trial_sec : float
        Maximum number of trials to get a new atomic position.
        If fails for more than "max_trial_sec" seconds, start to find new atoms from scratch.

    """

    ##
    backbone = backbone.copy()
    # Get spec_list and num_spec_dict
    if num_spec_dict is None:
        spec_list = backbone.get_chemical_symbols()
        num_spec_dict = list2count(spec_list)
        num_vacancy = 0

    # Get spec_list.
    else:
        if 'V' in num_spec_dict.keys():
            num_vacancy = num_spec_dict['V']
            del(num_spec_dict['V'])
        else:
            num_vacancy = 0

        spec_list = count2list(num_spec_dict)

    # Get num_fix_dict and num_shffl_spec_dict
    num_fix_dict = {}
    from copy import deepcopy
    num_shffl_spec_dict = deepcopy(num_spec_dict)
    if isinstance(fix_ind_dict, list) or isinstance(fix_ind_dict, np.ndarray):
        fix_ind_arr = np.array(deepcopy(fix_ind_dict))
        fixed_atoms = backbone.copy()[fix_ind_arr]
        fixed_spec = np.array(fixed_atoms.get_chemical_symbols())
        fix_ind_dict = {}
        for spec in np.unique(fixed_spec):
            fix_ind_dict[spec] = list(fix_ind_arr[fixed_spec==spec])
    if isinstance(fix_ind_dict, dict): 
        fix_ind_dict = deepcopy(fix_ind_dict)
        for key, value in fix_ind_dict.items():
            num_fix_dict[key] = len(value)
            num_shffl_spec_dict[key] -= len(value)
            fix_ind_dict[key] = np.array(value, dtype=np.int32).tolist()
            if key == 'V' and num_vacancy == 0:
                raise ValueError('The fix_ind_dict can not have "V", if num_spec_dict do not have "V"')
        if 'V' not in fix_ind_dict.keys():
            fix_ind_dict['V'] = []
            num_fix_dict['V'] = 0
            num_shffl_spec_dict['V'] = 0
    elif fix_ind_dict is None:
        fix_ind_dict = {}
        for key in num_spec_dict.keys():
            fix_ind_dict[key] = []
            num_fix_dict[key] = 0
        fix_ind_dict['V'] = []
        num_fix_dict['V'] = 0
    else:
        raise ValueError('Unknown type of fix_ind_dict.')

    # Get shffl_spec_list.
    shffl_spec_list = count2list(num_shffl_spec_dict)

    # Covalent bond length expectation value
    coval_expect = covalent_expect(spec_list)
                
    ## Cell strain adjustment.
    if strain_ratio is not None and strain is not None:
        raise ValueError("strain_ratio & strain parameters provided simultaneously. \
            Just provide one.")
    if strain is not None:
        strain = np.array(strain)
        if strain.shape != (3,):
            raise ValueError("Somethings wrong with strain parameter. Please check.")
        norm = np.linalg.norm(backbone.cell, axis=1)
        strain_ratio = strain / norm + 1
    if strain_ratio is not None:
        strain_ratio = np.array(strain_ratio)
        if strain_ratio.shape != (3,):
            raise ValueError("Somethings wrong with strain_ratio parameter. Please check.")
        backbone.set_cell(
            backbone.cell * np.expand_dims(strain_ratio, axis=1),
            scale_atoms = True,
            )
    if strain_ratio is None and strain is None:
        strain_ratio = [1.,1.,1.]
        backbone.set_cell(
            backbone.cell * np.expand_dims(strain_ratio, axis=1),
            scale_atoms = True,
            )

    ## Vacuum layer adjustment.
    if vacuum_ratio is not None and vacuum is not None:
        raise ValueError("vacuum_ratio & vacuum parameters provided simultaneously. \
            Just provide one.")
    if vacuum is not None:
        vacuum = np.array(vacuum)
        if vacuum.shape != (3,):
            raise ValueError("Somethings wrong with vacuum parameter. Please check.")
        norm = np.linalg.norm(backbone.cell, axis=1)
        vacuum_ratio = vacuum / norm + 1
    if vacuum_ratio is not None:
        vacuum_ratio = np.array(vacuum_ratio)
        if vacuum_ratio.shape != (3,):
            raise ValueError("Somethings wrong with vacuum_ratio parameter. Please check.")
        backbone.set_cell(
            backbone.cell * np.expand_dims(vacuum_ratio, axis=1),
            scale_atoms = False,
            )
    if vacuum_ratio is None and vacuum is None:
        vacuum_ratio = [1.,1.,1.]
        backbone.set_cell(
            backbone.cell * np.expand_dims(vacuum_ratio, axis=1),
            scale_atoms = True,
            )

    ## Determine cutoff radius.
    if cutoff_radi is not None and cutoff_frac is not None:
        raise ValueError("cutoff_radi & cutoff_frac parameters provided simultaneously. \
            Just provide one.")
    if cutoff_radi is not None:
        cutoff_r = cutoff_radi
    elif cutoff_frac is not None:
        cutoff_r = coval_expect * 2 * cutoff_frac
    else:
        cutoff_r = 0.

    ## Get random adjust radius
    from ase.build import make_supercell
    if random_frac is not None and random_radi is None:
        supercell = make_supercell(backbone,[[2,0,0],[0,2,0],[0,0,2]])
        from ase.optimize.precon.neighbors import estimate_nearest_neighbour_distance as rNN
        rdf_1st_peak = rNN(supercell)
        ran_radi = rdf_1st_peak / 2 * random_frac
        if log:
            print("")
            print("********* Please check carefully !!!! ***********".center(120))
            print("RDF 1st peak / 2 == {:.2f}".format(rdf_1st_peak/2).center(120))
            print("Positional deviation degree == {:.2f}".format(random_frac).center(120))
            print("==> Random deviation radius == {:.2f}".format(ran_radi).center(120))
            print("it is {:.2f} % of covalent-bond-length expectation value.".format(ran_radi / coval_expect * 100).center(120))
            print("")
            print("C.f. ) covalent-bond-length expectation value == {:.2f}".format(coval_expect).center(120))
            print("C.f. ) cutoff radius == {:.2f}".format(cutoff_r).center(120))
            print("C.f. ) cutoff radius / covalent bond expectation *2 == {:.2f} %".format(cutoff_r / coval_expect / 2 * 100).center(120))
            print("")
    elif random_radi is not None and random_frac is None:
        ran_radi = float(random_radi)
    else:
        raise ValueError('Check random_radi or random_frac parameters.')

    ## Main
    if num_vacancy != 0:
        # Choose vacancy indices.
        vacancy_ind = np.random.permutation(
            np.setdiff1d(
                range(len(backbone)),
                np.concatenate(list(fix_ind_dict.values())),
                True,
                ),
            )[:num_vacancy - num_fix_dict['V']]
        # Add fixed-vacancy indicies to the array.
        vacancy_ind = np.concatenate([vacancy_ind, fix_ind_dict['V']]).astype(int)

        # Remove vacancies from the backbone.
        vacancy_bool = np.array([True] *len(backbone))
        vacancy_bool[vacancy_ind] = False
        backbone = backbone[vacancy_bool]

        # Update fix_ind_dict.
        del(fix_ind_dict['V'])
        for key, value in fix_ind_dict.items():
            for i in range(len(value)):
                value[i] -= np.sum(vacancy_ind < value[i])
            fix_ind_dict[key] = value
    fix_ind_arr = np.concatenate(list(fix_ind_dict.values())).astype(np.int32)

    from time import time
    len_atoms = len(backbone)
    old_posi  = backbone.get_positions()
    cell      = backbone.get_cell()
    cell_inv  = np.linalg.inv(cell)

    new_posi = []
    while len(new_posi) < len_atoms:
        # Start time of this loop.
        time_i = time()
        while True:
            # Attach one more atom.
            new_posi.append(old_posi[len(new_posi)].copy())
            # Give positional deviation to the lastest atom.
            if not pin_the_fixed or len(new_posi)-1 not in fix_ind_arr:
                direc_vec = np.random.rand(3)-0.5
                direc_vec /= np.linalg.norm(direc_vec)
                new_posi[-1] += direc_vec * ran_radi
            # Get minimum distance from latest atom.
            rel_new_posi = np.array(new_posi) @ cell_inv
            if len(new_posi) != 1:
                min_dist = np.min(np.linalg.norm(((rel_new_posi[:-1] - rel_new_posi[-1] + np.array([0.5]*3)) % 1.0 - np.array([0.5]*3)) @ cell, axis=1))
            else:
                min_dist = cutoff_r + 1.
            # Get elapsed time of this loop.
            time_f = time()
            time_d = time_f - time_i
            # If the latest atom is properly positioned, break this loop.
            if min_dist > cutoff_r:
                if log:
                    if len(new_posi) % 100 == 0:
                        print("( {} th / {} ) new atom position found".format(len(new_posi), len_atoms))
                break
            # If the latest atom is too close to another atom, remove the latest one.
            elif time_d < max_trial_sec:
                new_posi.pop()
            # If failed for more than "max_trial_sec" seconds, restart from scratch.
            else:
                new_posi = []
                break
    new_atoms = backbone.copy()
    new_atoms.set_positions(new_posi, apply_constraint = False)
    # Shuffle positions
    shuffle_ind = np.setdiff1d(range(len(new_atoms)), np.concatenate(list(fix_ind_dict.values())), True)
    new_positions = new_atoms.get_positions().copy()
    new_positions[shuffle_ind] = np.random.permutation(new_positions[shuffle_ind])
    new_atoms.set_positions(new_positions, apply_constraint = False)
    
    # Correct chemical symbols
    new_species = np.array(['XX']*len(new_atoms))
    new_species[shuffle_ind] = shffl_spec_list
    if len(fix_ind_arr):
        new_species[fix_ind_arr] = count2list(num_fix_dict)

    ## Set the chemical symbols
    new_atoms.set_chemical_symbols(new_species)
    # Sort by chemical numbers
    if sort_atoms:
        new_atoms = new_atoms[np.argsort(new_atoms.get_atomic_numbers())]

    return new_atoms
