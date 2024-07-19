# Usage: python xyz2dp.py dataset.xyz
# sed 's/forces/force/g' train.xyz > dataset.xyz
import dpdata
import numpy as np
import sys

filename = sys.argv[1]
filetype = 'quip/gap/xyz'

ms = dpdata.MultiSystems.from_file(filename, fmt=filetype)
total_frames = sum(len(system) for system in ms)
print('# The dataset contains {} frames.'.format(total_frames))

# ms.to_deepmd_raw('./train')
ms.to_deepmd_npy('./train',set_size=total_frames)
