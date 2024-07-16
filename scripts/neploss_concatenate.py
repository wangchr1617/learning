import numpy as np
import sys
loss1 = np.loadtxt(sys.argv[1])
loss2 = np.loadtxt(sys.argv[2])
loss2[:, 0] = loss2[:, 0] + loss1[-1, 0]
loss = np.concatenate((loss1,loss2), axis=0)
np.savetxt('loss_all.out', loss)
