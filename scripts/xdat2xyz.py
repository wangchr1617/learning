#lipai@mail.ustc.edu.cn
#convert XDATCAR to unwraped xyz file
import numpy as np
from copy import deepcopy

xdatcar = open('XDATCAR', 'r')
xyz = open('NEP-dataset.xyz', 'w')

system = xdatcar.readline()
scale = float(xdatcar.readline().rstrip('\n'))
print(scale)

#get lattice vectors
a1 = np.array([ float(s)*scale for s in xdatcar.readline().rstrip('\n').split() ])
a2 = np.array([ float(s)*scale for s in xdatcar.readline().rstrip('\n').split() ])
a3 = np.array([ float(s)*scale for s in xdatcar.readline().rstrip('\n').split() ])

comment='Lattice=\"'+str(a1[0])+' '+str(a1[1])+' '+str(a1[2])
comment=comment+str(a2[0])+' '+str(a2[1])+' '+str(a2[2])
comment=comment+str(a3[0])+' '+str(a3[1])+' '+str(a3[2])+'\"'

#Read xdatcar
element_names = xdatcar.readline().rstrip('\n').split()

element_dict = {}
element_numbers = xdatcar.readline().rstrip('\n').split()

Natom = 0
Ntype = len(element_names)
Nname=[]
for t in range(Ntype):
    Natom += int(element_numbers[t])
    for i in range(int(element_numbers[t])):
        Nname.append(element_names[t])

print(Ntype,Natom)

f_prev=np.zeros([Natom,3])
f_next=np.zeros([Natom,3])
while True:
    line = xdatcar.readline()
    if len(line) == 0:
        break
    xyz.write(str(Natom) + "\n"+comment+"\n")

    for atom in range(Natom):
        p = xdatcar.readline().rstrip('\n').split()
        f_next[atom,:] = np.array([ float(s) for s in p ])
        for x in range(3):
            if(f_next[atom,x]-f_prev[atom,x]<-0.5):
                f_next[atom,x]+=1
            elif(f_next[atom,x]-f_prev[atom,x]>0.5):
                f_next[atom,x]-=1

        c_coords=f_next[atom,0]*a1+f_next[atom,1]*a2+f_next[atom,2]*a3
        xyz.write(Nname[atom]+" "+str(c_coords[0])+" "+str(c_coords[1])+" "+str(c_coords[2])+"\n")

    f_prev=deepcopy(f_next)

xdatcar.close()
xyz.close()
