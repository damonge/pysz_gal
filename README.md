# pysz_z
This is the python wrpaaer for the calculation of thermal SZ power spectrum and its cross correlation with galaxies.
The tSZ Cl code is originally devoloped by Eiichiro Komatsu.

# INSTALL
1. Move to source directory

```
cd pysz_z/pysz_z/source
```

2. Compile fortran codes

```
make clean && make
```

if you can not make it, modify Makefile in this directory


3. then go back to top directory of pysz and install python modules

```
cd PYSZDIRECTORY
python setup.py install
```
