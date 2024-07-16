import sys
import time
start_time = time.time()

from ase.build import sort
from ase.io import read, write
from ase.io.vasp import read_vasp_xdatcar
from ase.neighborlist import neighbor_list
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
from scipy.stats import gaussian_kde

def configure_plot():
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

def slicing(file_input, pmax=10, output_file=False):
    atoms = read(file_input)  
    cell_center = atoms.get_cell().sum(axis=0) / 2
    new_positions = atoms.get_positions() - cell_center
    atoms.set_positions(new_positions)
    pos = pd.DataFrame(atoms.get_positions(), columns=["x", "y", "z"])
    filtered_indices = pos[(pos.abs() <= pmax).all(axis=1)].index
    filtered_atoms = atoms[filtered_indices]
    if output_file:
        write('./selected.vasp', sort(filtered_atoms))
    return filtered_atoms

def find_ABC(center_atom_index, neighbor_indices):
    return np.array([[first_neighbor, center_atom_index, second_neighbor] 
                     for first_neighbor in neighbor_indices 
                     for second_neighbor in neighbor_indices 
                     if first_neighbor != second_neighbor])

def analyze_file(filename, filetype='POSCAR', pmax=20, cutoff=4, theta_min=170, theta_max=180, mic=True, output_file=False, frame_interval=10):
    if filetype == 'POSCAR':
        atoms = slicing(filename, pmax, output_file)
        return compute_ALTBC(atoms, cutoff, theta_min, theta_max, mic)
    elif filetype == 'XDATCAR':
        structures = read_vasp_xdatcar(filename, index=slice(0, None, frame_interval))
        if len(structures[0]) > 1000:
            print("Number of atoms exceeds 1000. Please use POSCAR file type for analysis.")
            sys.exit(1)
        data_frames = [compute_ALTBC(structure, cutoff, theta_min, theta_max, mic) for structure in structures]
        return pd.concat(data_frames, ignore_index=True) if data_frames else pd.DataFrame()
    else:
        print("Unrecognized file type!")
        sys.exit(1)

def compute_ALTBC(atoms, cutoff, theta_min, theta_max, mic):
    i_list, j_list = neighbor_list('ij', atoms, cutoff)
    distances = atoms.get_all_distances(mic=mic)
    num_atoms = len(atoms)
    all_indices = []
    all_normalized_weights = []
    
    for center_atom_index in range(num_atoms):
        neighbor_indices = j_list[i_list == center_atom_index]
        if len(neighbor_indices) > 1:
            indices = find_ABC(center_atom_index, neighbor_indices)
            all_indices.append(indices)
            all_normalized_weights.append(np.ones((len(indices), 1)) / (len(neighbor_indices) - 1))
    
    if not all_indices:
        return pd.DataFrame()

    indices = np.vstack(all_indices)
    normalized_weights = np.vstack(all_normalized_weights)
    angles = atoms.get_angles(indices, mic=mic)
    df = pd.DataFrame(np.hstack([indices, normalized_weights, angles[:, None]]), columns=["A", "B", "C", "weight", "angle_ABC"])
    df = df[(df["angle_ABC"] >= theta_min) & (df["angle_ABC"] <= theta_max)].reset_index(drop=True)
    df["AB"] = distances[df["A"].astype(int), df["B"].astype(int)]
    df["BC"] = distances[df["B"].astype(int), df["C"].astype(int)]
    df = df[(df["AB"] <= cutoff) & (df["BC"] <= cutoff)].reset_index(drop=True)
    df['pair'] = df.groupby(['A', 'B'])['weight'].transform('sum')
    return df

def plot_common_settings():
    plt.xlabel("Distance r1 (Å)")
    plt.ylabel("Distance r2 (Å)")
    plt.xlim(2.5, 3.8)
    plt.ylim(2.5, 3.8)
    plt.xticks(np.arange(2.6, 3.9, 0.2))
    plt.yticks(np.arange(2.6, 3.9, 0.2))
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()

def plot_results(df, plot_type='scatter'):
    x, y, z = np.array(df["AB"]), np.array(df["BC"]), np.array(df["pair"])
    xy = np.vstack([x, y])
    kde = gaussian_kde(xy)
    
    if plot_type == 'scatter':
        z_kde = kde(xy) * z * z
        plot_scatter(x, y, z_kde)
    elif plot_type == 'kde':
        plot_kde(kde, x, y, z)

def plot_scatter(x, y, z_kde):
    cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', ['blue', 'cyan', 'yellow', 'red'])
    plt.figure(figsize=(6, 5))
    scatter = plt.scatter(x, y, c=z_kde, cmap=cmap, s=1, alpha=1)
    plt.colorbar(scatter, label='P(r1,r2)')
    plot_common_settings()
    plt.show()

def plot_kde(kde, x, y, z):
    m1 = np.max(kde(np.vstack([x, y])) * z * z)
    x1, y1 = np.arange(2.5, 3.8, 0.01), np.arange(2.5, 3.8, 0.01)
    X, Y = np.meshgrid(x1, y1)
    kde_values = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    k = kde_values * m1 / np.max(kde_values)
    plt.figure(figsize=(6, 5))
    plt.imshow(k, extent=[2.5, 3.8, 2.5, 3.8], origin='lower', cmap=mcolors.LinearSegmentedColormap.from_list('custom_cmap', ['blue', 'cyan', 'yellow', 'red']), aspect='auto')
    plt.colorbar(label='P(r1,r2)')
    plot_common_settings()
    plt.show()

def main():
    configure_plot()
    # df = analyze_file('./XDATCAR', filetype='XDATCAR', frame_interval=10)
    df = analyze_file('./POSCAR', filetype='POSCAR', pmax=20, output_file=False)
    if not df.empty:
        plot_results(df, plot_type='scatter')  # 'scatter' or 'kde'

if __name__ == "__main__":
    main()

end_time = time.time()
print(f"Execution time: {end_time - start_time} seconds")
