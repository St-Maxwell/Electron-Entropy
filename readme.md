## Principle
* [_J. Chem. Theory Comput._ 2013, 9, 3165−3169](http://pubs.acs.org/doi/full/10.1021/ct400212t)
* [求解电子的热力学量](https://gensoukyo.me/-/eles)

## How to compile codes
```
gfortran -c m_global.f90 -o m_global.o
gfortran -c m_timer.f90 -o m_timer.o -fopenmp
gfortran -c m_function.f90 -o m_function.o -fopenmp
gfortran -c m_solve.f90 -o m_solve.o
gfortran -c ElectronEntropy.f90 -o ElectronEntropy.o
gfortran -o ElectronEntropy *.o -fopenmp
```