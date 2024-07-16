# Usage: python gb_info.py 1 1 0 30

from aimsgb import GBInformation
import sys

x1 = int(sys.argv[1])
x2 = int(sys.argv[2])
x3 = int(sys.argv[3])
sigma_max = int(sys.argv[4])
gb_info = GBInformation([x1, x2, x2], sigma_max, specific=False)

print(gb_info.__str__())
print(gb_info.axis)
