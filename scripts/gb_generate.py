# Usage: python gb_generate.py 0 0 1 5 1 2 0
from aimsgb import GrainBoundary, Grain
import sys

o1 = int(sys.argv[1])
o2 = int(sys.argv[2])
o3 = int(sys.argv[3])
sigma = int(sys.argv[4])
p1 = int(sys.argv[5])
p2 = int(sys.argv[6])
p3 = int(sys.argv[7])
filename = "GB_{}{}{}-{}-{}{}{}.vasp".format(o1, o2, o3, sigma, p1, p2, p3)
s_input = Grain.from_file("POSCAR")
gb = GrainBoundary([o1, o2, o3], sigma, [p1, p2, p3], s_input, uc_a=1, uc_b=1)
# gb.grain_a.translate_sites(range(len(gb.grain_a)), [0.2, 0, 0])
structure = Grain.stack_grains(gb.grain_a, gb.grain_b, direction=gb.direction)
structure.to(filename, "POSCAR")
