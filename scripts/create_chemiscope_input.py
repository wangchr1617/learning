# Usage: python create_chemiscope_input.py -a 8 ase_dataset.db descriptors.npy energies.npy
import argparse
import numpy as np
import chemiscope as cs

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from ase.io import read

if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('dataset',
		help='the ase dataset to use')

	parser.add_argument('descriptors',
		help='the descriptors for the dataset')

	parser.add_argument('energies',
		help='the energies for the dataset')

	parser.add_argument('-a', '--atomic', type=float,
		help='specify the cutoff for the local environment')

	args = parser.parse_args()

	print('Reading data!')
	frames = read(args.dataset, ':')
	descriptors = np.load(args.descriptors)
	energies = np.load(args.energies)

	print('Performing PCA!')
	scaler = StandardScaler()
	scaler.fit(descriptors)
	descriptors = scaler.transform(descriptors)

	pca = PCA(n_components=3, svd_solver='full')
	pca.fit(descriptors)
	descriptors_pc = pca.transform(descriptors)

	print('Creating chemiscope input file!')
	target = 'structure' if args.atomic is None else 'atom'
	properties = {
		'PCA': {
			'target': target,
			'values': descriptors_pc,
			'description': 'PCA of the descriptors',
		},
		'energies': {
			'target': target,
			'values': energies,
			'units': 'eV',
			'description': 'Energy per atom for the structures',
		},
	}

	## https://chemiscope.org/
	if target == 'structure':
		cs.write_input(
			path = 'chemiscope.json.gz',
			frames = frames,
			properties = properties
		)
	else:
		cs.write_input(
			path = 'chemiscope.json.gz',
			frames = frames,
			properties = properties,
			environments=cs.all_atomic_environments(frames, args.atomic)
		)

	print('DONE! upload chemiscope.json.gz to https://chemiscope.org/ to view')
