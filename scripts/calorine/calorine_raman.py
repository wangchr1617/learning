import numpy as np
from ase.io import read as ase_read

timestep = 1 #unit:fs
numba_installed = True

try:
    from numba import jit
except ImportError:
    numba_installed = False
    print("Warning: Numba is not installed, direct mode may be quite slow!")
    
result = np.loadtxt("polarizability_train.out")[:,:6]

frame = ase_read("train.xyz",index="0")
index=0
frames_numbers=0
with open("train.xyz") as f:
    for line in f:
        index+=1
        if "Lattice" in line:
            frames_numbers+=1

if index%(len(frame.get_atomic_numbers())+2)==0:
    atomic_numbers = [len(frame.get_atomic_numbers())]*frames_numbers
else:
    frames = ase_read("train.xyz",index=":")
    atomic_numbers = [len(i.get_atomic_numbers()) for i in frames]  
atomic_numbers = np.array(atomic_numbers)
tmp = (result*np.tile(atomic_numbers.reshape(-1,1),(1,6))).T

pol_data = np.array([tmp[0],tmp[3],tmp[5],tmp[3],tmp[1],tmp[4],tmp[5],tmp[4],tmp[2]]).T

# use only the first 10% of the data for the correlation function, as the rest is not statistically meaningful and our quadratic algorithm becomes too slow    
Nmax=len(pol_data)//10 

def conditional_decorator(condition):
    def decorator(func):
        if not condition:
            # Return the function unchanged, not decorated.
            return func
        return jit(func)
    return decorator
    
@conditional_decorator(numba_installed)
def calc_iso(data,Nmax):
    correlation = np.zeros(Nmax)
    for i in range(Nmax): #iter all correlation step
        for j in range(len(data)-Nmax): #number of origins
            correlation[i]+=(data[j])*(data[j+i])
    correlation/=correlation[0]
    return(correlation)

@conditional_decorator(numba_installed)
def calc_aniso(data,Nmax):
    correlation = np.zeros(Nmax) #the polarizability autocorrelation function (PACF)    
    for nc in range(Nmax): #loop over the correlation steps
        for m in range(len(data)-Nmax):
            origin = data[m]
            tau = data[m+nc]
            correlation[nc] += (2.0 * origin[0][0] * tau[0][0] - origin[0][0] * tau[1][1] 
                - origin[0][0] * tau[2][2] - origin[1][1] * tau[0][0] 
                - origin[2][2] * tau[0][0] + 6.0 * origin[0][1] * tau[0][1] 
                + 6.0 * origin[0][2] * tau[0][2] + 6.0 * origin[1][2] * tau[1][2] 
                + 2.0 * origin[1][1] * tau[1][1] + 2.0 * origin[2][2] * tau[2][2]
                - origin[1][1] * tau[2][2] - origin[2][2] * tau[1][1]) / (9.0) 

    correlation/=correlation[0] 
    return correlation

@conditional_decorator(numba_installed)
def compute_power_spectrum(acf,Nmax):
    raman = np.zeros(Nmax)
    for k in range(Nmax):
        raman[k] = (acf*np.cos(2*np.pi*k*np.arange(0,Nmax,1)/(2*Nmax-1))).sum()    
    return raman

############################  isotropic polarizability  ###########################
pol_data = pol_data.reshape((-1,3,3))
pol_data_iso = (pol_data[:,0,0]+pol_data[:,1,1]+pol_data[:,2,2])/3.0
correlation_iso = calc_iso(pol_data_iso-pol_data_iso.mean(),Nmax)

#write autocorrelation function
file_handle = open("iso_polar_correlation.time.txt",'w')
file_handle.write("# correlation in the time domain, first column time in fs\n")
for i in range(-Nmax+1,Nmax):
    file_handle.write("%f %f\n" %(i*timestep,correlation_iso[abs(i)]))
file_handle.close()

Kronecker_function = np.ones(Nmax) *2
Kronecker_function[0] = 1
Hann_window = (np.cos(np.pi*np.array(range(Nmax))/Nmax)+1)*0.5
correlation_iso = correlation_iso * Hann_window * Kronecker_function
raman_iso_intensity = compute_power_spectrum(correlation_iso,Nmax)

############################  anisotropic polarizability  ###########################
pol_data_aniso = pol_data - np.array([i *np.identity(3) for i in pol_data_iso])    
correlation_aniso = calc_aniso(pol_data_aniso-pol_data_aniso.mean(axis=0),Nmax)

#write autocorrelation function
file_handle = open("aniso_polar_correlation.time.txt",'w')
file_handle.write("# correlation in the time domain, first column time in fs\n")
for i in range(-Nmax+1,Nmax):
    file_handle.write("%f %f\n" %(i*timestep,correlation_aniso[abs(i)]))
file_handle.close()

correlation_aniso = correlation_aniso * Hann_window * Kronecker_function
raman_aniso_intensity = compute_power_spectrum(correlation_aniso,Nmax)

cm_array = np.arange(0,Nmax,1)/((2*Nmax-1)*timestep*10**(-15)*2.99792458*10**(10)) 
spectrum_iso = np.array([cm_array,1.0/3 * raman_iso_intensity*cm_array**2]).T
np.savetxt("frequency_reduced_raman_iso.txt",spectrum_iso)

spectrum_aniso = np.array([cm_array,raman_aniso_intensity*cm_array**2]).T
np.savetxt("frequency_reduced_raman_aniso.txt",spectrum_aniso)
