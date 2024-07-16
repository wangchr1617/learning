########################################################################################
# Python script to Get selected COHP and ICOHP from COHPCAR.lobster file              ##
# Written by Yafei Jiang                                                              ##
# Email: jiangyafei730@163.com                                                        ##
# Usage: python cohp.py                                                               ##
# update including spin unpolarization condition                                      ##
########################################################################################
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd 
from pandas import DataFrame


def process_file(file_name):
    '''Loading COHPCAR.lobster file and return numbers, titles and cohp data of interactions'''
    #file_name="COHPCAR.lobster"
    file_out=open(file_name,'r')
    lines=file_out.readlines()
    print("Reading COHPCAR.lobster...")
    num = int(lines[1].split()[0])
    spin = int(lines[1].split()[1])
    grid = int(lines[1].split()[2])
    titles=["Average"]
    titles.append(lines[3].split("(")[0].split(":")[1])
    for i in range(num-2):
        linesplit=lines[i+4].split("[")
        atom1=linesplit[1].split("]")[0]
        atom2=linesplit[2].split("]")[0]
        titles.append([atom1,atom2])
    
    data = []    
    for i in range(grid):
        data.append(list(map(float,lines[num+2+i].split())))
    
    cohp = np.array(data).T
    print(len(cohp))
    print("Including cohp of "+str(num))
    print(titles)
    # energy = cohp[0]
    print("COHP data is loaded.")
    return num,titles,cohp


def get_cohp(num,titles,cohp,orbital):
    '''readin interaction orbital paris and return cohp and icohp data of orbital paris'''
    if type(orbital[0]) is not list and type(orbital[1]) is not list:
        for i in range(num):
            if orbital[0] == "all":
                index0 = 1
            elif orbital[0] in titles[i][0] and orbital[1] in titles[i][1]:
                index0 = i
                break
            else:
                continue
        try:
            index0
        except NameError:
            print("Your input orbitals are not in COHPCAR.lobster!")
            return None,None,None,None
        else: 
            cohp_alpha = -cohp[index0*2+1]
            icohp_alpha = cohp[index0*2+2][int(np.argwhere(cohp[0]==0))]
            if len(cohp) == num * 4 + 1: # spin polarization
                cohp_beta = -cohp[num*2+index0*2+1]
                icohp_beta = cohp[num*2+index0*2+2][int(np.argwhere(cohp[0]==0))]  # icohp value at E-Ef=0eV
            else:  # spin unpolarization
                cohp_beta = np.zeros(len(cohp[0]))
                icohp_beta = 0
    elif type(orbital[0]) is list and type(orbital[1]) is not list:
        cohp_alpha = np.zeros(len(cohp[0]))
        cohp_beta = np.zeros(len(cohp[0]))
        icohp_alpha = 0
        icohp_beta = 0
        for j in range(len(orbital[0])):
            for i in range(num):
                if orbital[0][j] in titles[i][0] and orbital[1] in titles[i][1]:
                    index0 = i
            try:
                index0
            except NameError:
                print("Your input orbitals are not in COHPCAR.lobster!")
                return None,None,None,None
            else: 
                cohp_alpha_one = -cohp[index0*2+1]
                icohp_alpha_one = cohp[index0*2+2][int(np.argwhere(cohp[0]==0))]            
                cohp_alpha = cohp_alpha + cohp_alpha_one            
                icohp_alpha = icohp_alpha + icohp_alpha_one
                if len(cohp) == num * 4 + 1: # spin polarization
                    cohp_beta_one = -cohp[num*2+index0*2+1]
                    icohp_beta_one = cohp[num*2+index0*2+2][int(np.argwhere(cohp[0]==0))]
                    cohp_beta = cohp_beta + cohp_beta_one
                    icohp_beta = icohp_beta + icohp_beta_one
                else:  # spin unpolarization
                    cohp_beta = np.zeros(len(cohp[0]))
                    icohp_beta = 0
    elif type(orbital[0]) is not list and type(orbital[1]) is list:
        cohp_alpha = np.zeros(len(cohp[0]))
        cohp_beta = np.zeros(len(cohp[0]))
        icohp_alpha = 0
        icohp_beta = 0
        for j in range(len(orbital[1])):
            for i in range(num):
                if orbital[0] in titles[i][0] and orbital[1][j] in titles[i][1]:
                    index0 = i
            try:
                index0
            except NameError:
                print("Your input orbitals are not in COHPCAR.lobster!")
                return None,None,None,None
            else: 
                cohp_alpha_one = -cohp[index0*2+1]
                icohp_alpha_one = cohp[index0*2+2][int(np.argwhere(cohp[0]==0))]            
                cohp_alpha = cohp_alpha + cohp_alpha_one            
                icohp_alpha = icohp_alpha + icohp_alpha_one
                
                if len(cohp) == num * 4 + 1: # spin polarization
                    cohp_beta_one = -cohp[num*2+index0*2+1]
                    icohp_beta_one = cohp[num*2+index0*2+2][int(np.argwhere(cohp[0]==0))]
                    cohp_beta = cohp_beta + cohp_beta_one
                    icohp_beta = icohp_beta + icohp_beta_one
                else:  # spin unpolarization
                    cohp_beta = np.zeros(len(cohp[0]))
                    icohp_beta = 0
    else:
        cohp_alpha = np.zeros(len(cohp[0]))
        cohp_beta = np.zeros(len(cohp[0]))
        icohp_alpha = 0
        icohp_beta = 0
        for j in range(len(orbital[0])):
            for k in range(len(orbital[1])):
                for i in range(num):
                    if orbital[0][j] in titles[i][0] and orbital[1][k] in titles[i][1]:
                        index0 = i
                try:
                    index0
                except NameError:
                    print("Your input orbitals are not in COHPCAR.lobster!")
                    return None,None,None,None
                else: 
                    cohp_alpha_one = -cohp[index0*2+1]
                    icohp_alpha_one = cohp[index0*2+2][int(np.argwhere(cohp[0]==0))]                
                    cohp_alpha = cohp_alpha + cohp_alpha_one                
                    icohp_alpha = icohp_alpha + icohp_alpha_one
                    
                    if len(cohp) == num * 4 + 1: # spin polarization
                        cohp_beta_one = -cohp[num*2+index0*2+1]
                        icohp_beta_one = cohp[num*2+index0*2+2][int(np.argwhere(cohp[0]==0))]
                        cohp_beta = cohp_beta + cohp_beta_one
                        icohp_beta = icohp_beta + icohp_beta_one
                    else:  # spin unpolarization
                        cohp_beta = np.zeros(len(cohp[0]))
                        icohp_beta = 0
    return cohp_alpha,cohp_beta,icohp_alpha,icohp_beta


def get_orb():
    '''get orbital paris from user input and return corresponding orbital pairs in format of COHPCAR.lobster'''
    print("=="*30)
    orb1=["all","1s","2s","3s","4s","5s","6s","7s",\
    "2px","2py","2pz","3px","3py","3pz","4px","4py","4pz","5px","5py","5pz","6px","6py","6pz","7px","7py","7pz",\
    "3dxy","3dyz","3dxz","3dz2","3dx2-y2","4dxy","4dyz","4dxz","4dz2","4dx2-y2",\
    "5dxy","5dyz","5dxz","5dz2","5dx2-y2","6dxy","6dyz","6dxz","6dz2","6dx2-y2",\
    "4fy3x2","4fxyz","4fyz2","4fz3","4fxz2","4fzx2","4fx3",\
    "5fy3x2","5fxyz","5fyz2","5fz3","5fxz2","5fzx2","5fx3"]
    orb2=["all","1s","2s","3s","4s","5s","6s","7s",\
    "2p_x","2p_y","2p_z","3p_x","3p_y","3p_z","4p_x","4p_y","4p_z","5p_x","5p_y","5p_z","6p_x","6p_y","6p_z","7p_x","7p_y","7p_z",\
    "3d_xy","3d_yz","3d_xz","3d_z^2","3d_x^2-y^2","4d_xy","4d_yz","4d_xz","4d_z^2","4d_x^2-y^2",\
    "5d_xy","5d_yz","5d_xz","5d_z^2","5d_x^2-y^2","6d_xy","6d_yz","6d_xz","6d_z^2","6d_x^2-y^2",\
    "4f_y(3x^2-y^2)","4f_xyz","4f_yz^2","4f_z^3","4f_xz^2","4f_z(x^2-y^2)","4f_x(x^2-3y^2)",\
    "5f_y(3x^2-y^2)","5f_xyz","5f_yz^2","5f_z^3","5f_xz^2","5f_z(x^2-y^2)","5f_x(x^2-3y^2)"]
    orb_p=["2p","3p","4p","5p","6p","7p"]
    orb_d=["3d","4d","5d","6d"]
    orb_f=["4f","5f"]
    orb_p2=["_x","_y","_z"]
    orb_d2=["_xy","_yz","_xz","_z^2","_x^2-y^2"]
    orb_f2=["_y(3x^2-y^2)","_xyz","_yz^2","_z^3","_xz^2","_z(x^2-y^2)","_x(x^2-3y^2)"]
    # input orbital paris: such as "sum sum;4s 4s"
    inp = input("Please select interactions: \n all\n s px py pz\n dxy dyz dxz dz2 dx2-y2\n fy3x2 fxyz fyz2 fz3 fxz2 fzx2 fx3\n such as: 3dyz 2py \n such as: 3d 3d \n such as: all all;3d 3d;4s 4s \n").split(";")
    if inp[0] == "":
        status = 0
        print("End input!")
        return status,status,inp[0]
    else:
        status = 1
        inp2= [i.split() for i in inp]  # input orbital name for each orbial pairs
        orbitals=[]  # orbital name in COHPCAR.lobster
        for ss in inp2:  ## ss: input orbital name list for each orbial pairs
            orb_pairs=[]
            if ss[0] == 'all' and ss[1] == 'all':
                orb_pairs = ss
            elif ss[0] == 'all':  # such as all 2p, extend all into all orbitals list of atom 1
                orb_pairs.append(orbA)
                if ss[1] in orb_p:  # such as 2p
                    orb_temp=[ss[1]+i for i in orb_p2]  # 2p --> [2p_x, 2p_y, 2p_z]
                    orb_pairs.append(orb_temp)
                elif ss[1] in orb_d:  # such as 3d
                    orb_temp=[ss[1]+i for i in orb_d2]  # 3d --> ["3d_xy","3d_yz","3d_xz","3d_z^2","3d_x^2-y^2"]
                    orb_pairs.append(orb_temp)
                elif ss[1] in orb_f:  # such as 4f
                    orb_temp=[ss[1]+i for i in orb_f2]  # 4f --> [...]
                    orb_pairs.append(orb_temp)
                else:
                    orb_pairs.append(orb2[orb1.index(ss[1])])  # 2px --> 2p_x
            elif ss[1] == 'all':  # such as 2p all, extend all into all orbitals list of atom 2
                if ss[0] in orb_p:  # such as 2p
                    orb_temp=[ss[0]+i for i in orb_p2]  # 2p --> [2p_x, 2p_y, 2p_z]
                    orb_pairs.append(orb_temp)
                elif ss[0] in orb_d:  # such as 3d
                    orb_temp=[ss[0]+i for i in orb_d2]  # 3d --> ["3d_xy","3d_yz","3d_xz","3d_z^2","3d_x^2-y^2"]
                    orb_pairs.append(orb_temp)
                elif ss[0] in orb_f:  # such as 4f
                    orb_temp=[ss[0]+i for i in orb_f2]  # 4f --> [...]
                    orb_pairs.append(orb_temp)
                else:
                    orb_pairs.append(orb2[orb1.index(ss[0])])  # 2px --> 2p_x 
                orb_pairs.append(orbB)
            else:
                for sss in ss:  ## sss:input orbital name in input orbital name list for each orbial pairs
                    
                    if sss in orb_p:  # such as 2p
                        orb_temp=[sss+i for i in orb_p2]  # 2p --> [2p_x, 2p_y, 2p_z]
                        orb_pairs.append(orb_temp)
                    elif sss in orb_d:  # such as 3d
                        orb_temp=[sss+i for i in orb_d2]  # 3d --> ["3d_xy","3d_yz","3d_xz","3d_z^2","3d_x^2-y^2"]
                        orb_pairs.append(orb_temp)
                    elif sss in orb_f:  # such as 4f
                        orb_temp=[sss+i for i in orb_f2]  # 4f --> [...]
                        orb_pairs.append(orb_temp)
                    else:
                        orb_pairs.append(orb2[orb1.index(sss)])  # 2px --> 2p_x
            orbitals.append(orb_pairs)
        print(orbitals)
        return status,orbitals,inp2  # status,orbital name list in COHPCAR.losber,input orbital name list for each orbital pairs

    
alpha,beta=chr(945),chr(946)
num,titles,cohp = process_file("COHPCAR.lobster")

orbA = list(set([i[0] for i in titles[2:]]))  # all orbitals in atom 1
orbB = list(set([i[1] for i in titles[2:]]))  # all orbitals in atom 2

status1 = 1
while status1 == 1:
    status1,orb_cohp,orb_inp=get_orb()
    if status1 ==1:
        cohps=[]
        icohps=[]
        orbs=[]
        status2 = 1
        for i in range(len(orb_inp)):
            cohp_a,cohp_b,icohp_a,icohp_b=get_cohp(num,titles,cohp,orb_cohp[i])
            if cohp_a is None:
                status2 = 0
                break
            else:
                if len(cohp) == num * 4 + 1: # spin polarization
                    print("The ICOHP values of %s are \n   alpha \t beta\n %f\t%f" %("-".join(orb_inp[i]),icohp_a,icohp_b))
                    cohps.append(cohp_a)
                    cohps.append(cohp_b)
                    icohps.append(icohp_a)
                    icohps.append(icohp_b)
                    orbs.append("_".join(orb_inp[i])+"_"+alpha)
                    orbs.append("_".join(orb_inp[i])+"_"+beta)
                else:  # spin unpolarization
                    cohps.append(cohp_a)
                    icohps.append(icohp_a)
                    orbs.append("_".join(orb_inp[i]))
                    print("The ICOHP values of %s is %f.\n" %("-".join(orb_inp[i]),icohp_a))
        if status2 == 1:    
            df = pd.DataFrame(np.array(cohps).T,index = cohp[0].T,columns=orbs)
            df.plot()
            plt.xlim(-8,5)
            plt.xlabel(r"$E-E_f\ (eV)$",fontsize=12)
            plt.ylabel("-COHP",fontsize=12)
            plt.tight_layout()
            plt.show()
            check=input("Do you want to save the cohp data: y or [n]\n")
            if check=="y":
                cohp_set = np.vstack((cohp[0],np.array(cohps))).T
                header="{:^12s}".format("E-Ef(eV)") + " ".join(map(lambda x: "{:^12s}".format(x),orbs)) + "\n" \
                + "{:^12s}".format("ICOHP") + " ".join(map(lambda x: "{:^12.5f}".format(x),icohps))
                output="cohp_"+'_'.join('-'.join(inner) for inner in orb_inp)+".dat"
                np.savetxt(output,cohp_set,fmt="%12.5f",header=header)
                print("Selected COHP data has been written into %s file." %output)
            else:
                pass
        else:
            pass
    else:
        pass

