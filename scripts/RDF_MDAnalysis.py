#Calculate the rdf between A-B atoms, allowing for dynamic selection of regions
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF


u = mda.Universe('npt0.tpr', '10-20ns.xtc') #Topology and trajectory


select_A = u.select_atoms('name OW and (prop z > 30 and prop z < 36)', updating=True)  #Select A atoms
select_B = u.select_atoms('name OW and (prop z > 30 and prop z < 36)', updating=True)  #Select B atoms, Can be different from A


rdf = InterRDF(select_A, select_B, range=(0.0, 20.0), nbins=200, n_threads=12)
rdf.run()

rdf.results.rdf[0] = 0.0

np.savetxt("rdf.txt", np.column_stack((rdf.results.bins, rdf.results.rdf)), header="r (Angstrom)\t g(r)")