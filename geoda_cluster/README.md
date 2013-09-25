This directory holds testing code for exploring parallelization on the GeoDa Center cluster


## Software Installation

The following components are being explored

 - [Anaconda](http://continuum.io/downloads) Python Distribution
 - [Pypar](https://github.com/daleroberts/pypar)


Anaconda python can be installed in user directories.
```
chmod +x Anaconda-1.7.0-Linux-x86_64.sh
./Anaconda-1.7.0-Linux-x86_64.sh
```

After Anaconda is installed pypar can be installed as follows:

```
git clone https://github.com/daleroberts/pypar.git
cd pypar
use openmpi-1.6.4
python setup.py install
```






