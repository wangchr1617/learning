# Usage: python capacity.py base/ Li Na 16
import os
import pandas as pd
from analyzer import BatteryAnalyzer
from pymatgen.io.vasp.inputs import Structure
from sys import argv
path = argv[1]
cations = argv[2:-1]
num_cations = eval(argv[-1])
slab_files = [os.path.join(path, i) for i in os.listdir(path)]
cap_dict = {}
for ion in cations:
    ion_dict = {}
    for file in slab_files:
        structure = Structure.from_file(file)
        slab_name = structure.formula.replace(' ', '')
        ele_list = structure.composition.elements
        # mass = structure.composition.weight
        s = [ele.name for ele in ele_list]
        
        oxi_dict = {}
        for ele in s:
            oxi_dict[ele] = +1
        oxi_dict["B"] = -2
        print(oxi_dict, ' -->')
        
        structure.add_oxidation_state_by_element(oxi_dict)
        battery = BatteryAnalyzer(struc_oxid=structure, cation=ion)
        capacity = battery.get_max_capgrav(num_cations) # return max grav capacity in mAh/g
        print('--------------------------------------')
        ion_dict[slab_name] = round(capacity, 2)
    cap_dict[ion] = ion_dict
df = pd.DataFrame.from_dict(cap_dict, orient='index').T
df.to_csv('capacity.csv')
