# Usage: python dp_descriptors.py frozen_model.pb ase_dataset.db Ge Te -a
import argparse
import numpy as np

from deepmd.infer import DeepPot
from ase.io import read
from tqdm import tqdm

if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('model',
		help='the deepmd model to use')

	parser.add_argument('dataset',
		help='the ase dataset to use')

	parser.add_argument('type_map', nargs='+', type=str,
		help='the type_map of the deepmd model')

	parser.add_argument('-a', '--atomic', action='store_true',
		help='store per atom descriptors and energies')

	args = parser.parse_args()

	dp = DeepPot(args.model)
	images = read(args.dataset, ':')

	descriptors = []
	energies = []

	for atoms in tqdm(images):
		coords = atoms.get_positions().reshape([1, -1, 3])
		cells = np.array(atoms.get_cell()).reshape([1, -1])
		symbols = atoms.get_chemical_symbols()
		atom_types = [args.type_map.index(s) for s in symbols]

		descriptor = dp.eval_descriptor(coords, cells, atom_types)

		if args.atomic:
			# save descriptors as float32 to reduce size
			descriptors.append(np.sum(descriptor, axis=0).astype(np.float32))
			energy = dp.eval(coords, cells, atom_types, atomic=True)[3].flatten()
			energies.append(energy)
		else:
			descriptors.append(np.sum(np.sum(descriptor, axis=0), axis=0))
			energy = atoms.get_potential_energy()
			energies.append(energy / len(atoms))

	if args.atomic:
		descriptors = np.concatenate(descriptors)
		energies = np.concatenate(energies)
	else:
		descriptors = np.array(descriptors)
		energies = np.array(energies)

	np.save('descriptors.npy', descriptors)
	np.save('energies.npy', energies)
