# Usage: python dp2ase.py train/ Ge Te
import os
import uuid
import shutil
import argparse

from dpdata import LabeledSystem
from ase.io import write

def to_ase(data):
	'''Convert System to ASE Atoms object.'''
	from ase import Atoms
	from ase.calculators.singlepoint import SinglePointCalculator

	structures = []
	species = [data['atom_names'][tt] for tt in data['atom_types']]

	for ii in range(data['coords'].shape[0]):
		structure = Atoms(
			symbols=species,
			positions=data['coords'][ii],
			pbc=not data.get('nopbc', False),
			cell=data['cells'][ii]
		)

		results = {
			'energy': float(data["energies"][ii]),
			'forces': data["forces"][ii]
		}
		if "virials" in data:
			# convert to GPa as this is ase convention
			# v_pref = 1 * 1e4 / 1.602176621e6
			vol = structure.get_volume()
			# results['stress'] = data["virials"][ii] / (v_pref * vol)
			results['stress'] = -data["virials"][ii] / vol

		structure.calc = SinglePointCalculator(structure, **results)
		structures.append(structure)

	return structures

if __name__ == '__main__':

	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('path', type=str,
		help='path to the folder containing all the data')

	parser.add_argument('type_map', nargs='+', type=str,
		help='the type_map of the data')

	args = parser.parse_args()

	print('\nReading data from: %s' % args.path)
	images = []
	for subdir, dirs, files in os.walk(args.path):
		for file in files:
			if file == 'type.raw':
				print(subdir)
				images += to_ase(LabeledSystem(
					subdir,
					fmt='deepmd/npy',
					type_map=args.type_map
				).data)

	print('\nWriting dataset to file: %s' % 'ase_dataset.db')
	filename = uuid.uuid4().hex
	write('/dev/shm/' + filename + '.db', images)
	shutil.move('/dev/shm/' + filename + '.db', 'ase_dataset.db')
