##Building##

The cython code needs to be compiled on the system that it will run on.  The contents of this directory have been compiled on an Ubuntu 12.04, 64bit box and may be portable to other architectures.

To build: `python setup.py build_ext --inplace`

We see two build warning that can be safely ignored.

To run: `python lmorans.py`  

The code utilizes 
