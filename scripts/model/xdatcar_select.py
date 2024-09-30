# Usage: python xdatcar_select.py ./XDATCAR 10000 15000 1
from ase.io import read, write
import sys
begin = sys.argv[2]
end = sys.argv[3]
step = sys.argv[4]
idx = '{}:{}:{}'.format(begin,end,step)
traj = read(sys.argv[1], index=idx)
write('./XDATCAR_new', traj)
