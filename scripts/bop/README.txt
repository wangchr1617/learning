To compile the BOP calculation code, simply type:

$make

Then execute the command:

$./bop

If you see the output of the usage message: 

Usage: bop cutoff average?(1:0)

It means that you have compiled the code successfully. To use this code, you have to first generate a configuration file in the XYZ format, then type in the execution command with two input arguments:

1) cutoff: this is the cutoff distance for the nearest neighbor. Namely, how far apart for two atoms you still consider as two neighboring atoms. Only the bonds between the neighboring atoms are included in the BOP calculations. I used 3.8 Angstrom for gold atoms.

2) average: whether you want to calculate the globally averaged BOPs (only one set of BOP values will be given) (1) or local BOP values for each atom (0). 

The input configuration file and the output data should be read and saved by the redirection operators.

For example, if I want to obtain the global BOP values of a gold nanocluster, I would type in:

./bop 3.8 1 < cube.1372.xyz > gbop.dat

If I want to obtain the local BOP values for each atom, I would type in:

./bop 3.8 0 < cube.1372.xyz > lbop.dat

I hope the above is helpful to you. Please do not forget to cite our JCP paper if you have publishable results by using our code:

Yanting Wang, S. Teitel, and Christoph Dellago, Melting of Icosahedral Gold Nanoclusters from Molecular Dynamics Simulations J. Chem. Phys. 122, 214722-214738 (2005).

Thanks a lot for using our BOP calculation code.

使用之前加载一下:module load gcc/9.3
