########################################################################################
# Python script to calculate RDF of specified atom and orbital                        ##
# Written by Yafei Jiang                                                              ##
# Email: jiangyafei730@163.com                                                        ##
# Usage: python RDF_analysis.py filename orbital                                               ## 
# Example: python RDF_analysis.py out.Os 5d                                                           ##
# Note: the file such as out.Os should be obtained from ADF calculation               ##                                                           ##
########################################################################################
import sys
import subprocess

## 1. Extract (basis "z", coefficient "c", mainquantumnumber "N" from output file of ADF
file_name=sys.argv[1]
orbital=sys.argv[2]
file_out=open(file_name,'r')
lines=file_out.readlines()

status1,line1 = subprocess.getstatusoutput("sed -n '/Valence Basis Sets/=' " + file_name)
# print(status1)
line_basis=int(line1.split('\n')[0])
Nbasis=int(lines[line_basis-1].split()[-1].split('\n')[0])  # Numbers of valence basis sets
basis=[]  # all zeta value of valence basis set 
mainquantumnumber=[]  # corresponding main quantum number
orbitals=[]  # corresponding orbital: S P D F
for i in range(Nbasis):
    basis.append(float(lines[line_basis+1+i].split()[2]))
    mainquantumnumber.append(int(lines[line_basis+1+i].split()[0]))
    orbitals.append(lines[line_basis+1+i].split()[1])

# status2,line2 = subprocess.getstatusoutput("sed -n '/=== S ===/=' " + file_name)
# status3,line3 = subprocess.getstatusoutput("sed -n '/=== P:y ===/=' " + file_name)
# status4,line4 = subprocess.getstatusoutput("sed -n '/=== D:xz ===/=' " + file_name)
# status5,line5 = subprocess.getstatusoutput("sed -n '/=== F:xyz ===/=' " + file_name)

orbital_dic={"S":0,"P":1,"D":2,"F":3}
# orbital="6p"
orbital_n = int(orbital[0])  # The main quantum number of orbital you input
orbital_ch = orbital[1].upper()  # The angular quantum number of orbital you input
orbital_l = orbital_dic[orbital_ch] # The angular quantum number of orbital you input

c=[]  #coefficient matrix
if orbital_ch == "F":
    status5,line5 = subprocess.getstatusoutput("sed -n '/=== F:xyz ===/=' " + file_name)
    if status5 == 0:
        column = orbital_n - orbital_l
        line_coeff=int(line5.split('\n')[1])
        f_index = [i for i in range(len(orbitals)) if orbitals[i] == orbital_ch]
        z=[basis[i] for i in f_index]
        N=[mainquantumnumber[i] for i in f_index]
        for i in range(len(z)):
            c.append(float(lines[line_coeff+6+i].split()[column]))

elif orbital_ch == "D":
    status4,line4 = subprocess.getstatusoutput("sed -n '/=== D:xz ===/=' " + file_name)
    if status4 == 0:
        column = orbital_n - orbital_l
        line_coeff=int(line4.split('\n')[1])
        d_index = [i for i in range(len(orbitals)) if orbitals[i] == orbital_ch]
        z=[basis[i] for i in d_index]
        N=[mainquantumnumber[i] for i in d_index]
        for i in range(len(z)):
            c.append(float(lines[line_coeff+6+i].split()[column]))

elif orbital_ch == "P":
    status3,line3 = subprocess.getstatusoutput("sed -n '/=== P:y ===/=' " + file_name)
    if status3 == 0:
        column = orbital_n - orbital_l
        line_coeff=int(line3.split('\n')[1])
        p_index = [i for i in range(len(orbitals)) if orbitals[i] == orbital_ch]
        z=[basis[i] for i in p_index]
        N=[mainquantumnumber[i] for i in p_index]
        if column < 5:
            for i in range(len(z)):
                c.append(float(lines[line_coeff+6+i].split()[column]))
        else:
            column = column - 4
            for i in range(len(z)):
                c.append(float(lines[line_coeff+6+i+len(z)+3].split()[column]))

elif orbital_ch == "S":
    status2,line2 = subprocess.getstatusoutput("sed -n '/=== S ===/=' " + file_name)
    if status2 == 0:
        column = orbital_n - orbital_l
        line_coeff=int(line2.split('\n')[1])
        s_index = [i for i in range(len(orbitals)) if orbitals[i] == orbital_ch]
        z=[basis[i] for i in s_index]
        N=[mainquantumnumber[i] for i in s_index]
        if column < 5:
            for i in range(len(z)):
                c.append(float(lines[line_coeff+6+i].split()[column]))
        else:
            column = column - 4
            for i in range(len(z)):
                c.append(float(lines[line_coeff+6+i+len(z)+3].split()[column]))

else:
    print("The orbital you input is not correct.")
    sys.exit()
	
####################################################################################
# 2. Calculate RDF of specified atom and orbital
####################################################################################
import numpy as np
import math


# z=np.loadtxt('basis.txt')  #sita
# c=np.loadtxt('coefficient.txt')  #coefficient matrix
# N=np.loadtxt('mainquantumnumber.txt')  #main quantum number of sita
M=len(z);

# AO nomalization coefficient
b=np.zeros(M)
for i in range(M):
    b[i]=(2 * z[i])**(N[i] + 0.5) / (math.factorial(2 * N[i]))**0.5
    
# MO nomalization coefficient
V=0
for i in range(M):
    for j in range(M):
        V=V + math.factorial(N[i]+N[j]) * b[i] * b[j] * c[i] * c[j] / (z[i] + z[j]) ** ( N[i] + N[j] + 1)

print(b)
print(V)
grid=800        
D_r=np.zeros(grid)  
r=np.zeros(grid)    
for m in range(grid):
    r[m]= m / 100  # radial variable 
    sumsum=0;
    for p in range(M):
        for q in range(M):
            sumpq=(1/V)*r[m]**2*c[p]*c[q]*b[p]*b[q]*(r[m]**(N[p]-1)* np.exp(-z[p]*r[m]))*(r[m]**(N[q]-1)* np.exp(-z[q]*r[m]));
        # AO:f1(r)=r^(N-1)*exp(-z(i,1)*r)  MO:f2(r)=E*sum(C(1,i)*f1)
        # F(r)=r^2*f2^2
            sumsum=sumsum+sumpq;
    D_r[m]=sumsum

bohr2A=0.529
with open(file_name.split(".")[1] + "_" + orbital + "-RDF.dat", 'w') as f:
    for i in range(len(r)):
        f.write("{:^10.6f}{:^10.6f}\n".format(r[i] * bohr2A, D_r[i]))

print("RDF calcuation is finished.")        	