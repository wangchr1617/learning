from pynep.io import load_nep,dump_nep
import random
import sys

filename=sys.argv[1]
train_ratio = 1.0
rand = True

train_data = load_nep(filename, ftype="exyz")

nframes = [i for i in range(len(train_data))]
if rand: random.shuffle(nframes)

train_frame = nframes[:int(len(train_data)*train_ratio)]
dump_nep('./train.xyz'.format(filename), [train_data[i] for i in train_frame], ftype="exyz")

test_frame = nframes[int(len(train_data)*train_ratio):]
dump_nep('./test.xyz'.format(filename), [train_data[i] for i in test_frame], ftype="exyz")
