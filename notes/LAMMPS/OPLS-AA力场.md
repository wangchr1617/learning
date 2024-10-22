

### lammps中设置OPLS-AA立场

键结作用：harmonics键长和键角势能，Fourier扭矩势能

非键作用：LJ势和Coulomb势

分别通过以下的lammps参数设置：

```
atom_style full
bond_style harmonic
angle_style harmonic
dihedral_style opls
pair_style lj/cut/coul/long
pair_modify mix geometric
special_bonds lj/coul 0.0 0.0 0.5
```

