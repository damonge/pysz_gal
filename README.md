# pysz_z
This is the python code for the calculation of the thermal SZ power spectrum and its cross correlation with galaxies.
The tSZ Cl calculation is based on [the fortran code](https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumdks/) devoloped by Eiichiro Komatsu.

** NOTICE **
The code is under developement. Please ask me if you want to use the code in your research project.

# INSTALL
1. Move to source directory

```
cd pysz_gal/pysz_gal/source
```

2. Compile fortran codes

```
make clean && make
```

if you can not make it, modify Makefile in this directory


3. then go back to top directory of pysz_gal and install python modules

```
cd PYSZDIRECTORY
python setup.py install
```

4. sample code

```
import pysz_gal.pysz_gal as szgal
calc_cls = szgal.tsz_gal_cl() # init class
## Set params ##
ell = np.logspace(1,np.log10(1500),10)
dndz = np.loadtxt('./dndz_example.txt')
h0 = 0.67
ombh2 = 0.02225
omch2 = 0.1198
lnAs = 3.065
ns = 0.9645
mnu = 0.06
bias = 1.55
Mcut = 10**11.8
M1 = 10**11.73
kappa = 1.0
sigma_Ncen = 0.15
alp_Nsat = 0.77
rmax = 4.39
rgs = 1.17
pars = {'h0':h0, 'obh2':ombh2,'och2':omch2,\
        'As':np.exp(lnAs)*1e-10,'ns':ns,'mnu':mnu,\
        'mass_bias':bias,\
        'Mcut':Mcut,'M1':M1,'kappa':kappa,'sigma_Ncen':sigma_Ncen,'alp_Nsat':alp_Nsat,\
        'rmax':rmax,'rgs':rgs,\
        'flag_nu':True,'flag_tll':False}
## get Cls ##
gg, gy, tll = szgal_cl.get_tsz_cl(ell,pars,dndz,1e-5,1.0)
```
