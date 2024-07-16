# Usage:
# Data from OSZICAR
'''
echo "lat, ene" > ev.csv
for i in */;
do
E=`tail -1 OSZICAR | awk '{print $5}'`
echo "$i,$E" >> ev.csv 
done
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def Birch_Murnaghan(p,x):
  E0, V0, B0, B1 = p
  return E0+(9*V0*B0/16)*((((V0/(x**3))**(2/3)-1)**3)*B1+\
         (((V0/(x**3))**(2/3)-1)**2)*(6-4*(V0/(x**3))**(2/3)))

def error(p,x,y):
  return Birch_Murnaghan(p,x) - y

f = "./ev.csv"
df = pd.read_csv(f)
x = np.array(df['lat'])
y = np.array(df['ene'])
p0 = [min(y), min(x)**3, 100, 1]
para = leastsq(error, p0, args=(x,y))
y_fitted = Birch_Murnaghan(para[0], x)
print("Emin=",para[0][0],"opt_lat=",para[0][1]**(1/3))
print("B0=",para[0][2],"B1=",para[0][3])
plt.figure(figsize=(8,6), dpi=300)
plt.scatter(x,y,20,c='r',label="DFT")
plt.plot(x,y_fitted,'-b',label="B-M Fitting")
plt.xlabel("Lattice constant", fontsize=20)
plt.ylabel("Energy(eV)", fontsize=20)
plt.legend()
plt.savefig('./ev.png', bbox_inches='tight')
