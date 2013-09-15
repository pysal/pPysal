'''
This is for the exploration of dynamically computing
the required portions of the diameter matrix for the
dynamic population of the error matrix.

The idea is that the error matrix is only (nxk), so
scalability in the memory domain will be better than
storing the entire (nxn) distance matrix.
'''

import multiprocessing as mp
import ctypes
import numpy as np


def allocate(values, classes):
    numvalues = len(values)

    errctype = mp.RawArray(ctypes.c_double,numvalues * classes)
    errormatrix = np.frombuffer(errctype)
    errormatrix.shape = (classes, numvalues)

    arrRow = np.array([values])
    n = np.arange(numvalues) + 1
    errormatrix[0] =  ((np.cumsum(np.square(arrRow))) - \
                ((np.cumsum(arrRow)*np.cumsum(arrRow)) / (n)))

    pivotmatrix = np.ndarray((1, classes), dtype=np.float)

    initarr(errormatrix)
    return pivotmatrix

def err(row, y, step, lenrow):
    stop = y + step
    if stop + 1 > lenrow:
        stop = lenrow - 1
    while y <= stop:
        err = sharederror[row-1][row-1:y+row]
        print err
        y += 1

def initarr(errormatrix_):
    global sharederror
    sharederror = errormatrix_

def init_error_row(errorrow_):
    global sharedrow
    sharedrow = errorrow_

def main(values, classes, sort=True):
    if sort:
        values.sort()

    pivotmatrix = allocate(values, classes)
    print sharederror
    cores = mp.cpu_count()

    #Naive partitioning
    row = 1
    jobs = []
    for x in sharederror[1:]:
        errorrow = x[row:]
        init_error_row(errorrow)
        step = len(errorrow) // cores
        for y in range(0,len(errorrow), step+1):
            p = mp.Process(target=err,args=(row, y, step, len(errorrow)))
            jobs.append(p)

        for j in jobs:
            j.start()
        for j in jobs:
            j.join()
        del jobs[:], p, j
        row += 1
if __name__ == "__main__":
    values = [12,10.8, 11, 10.8, 10.8, 10.8, 10.6, 10.8, 10.3, 10.3, 10.3,10.4, 10.5, 10.2, 10.0, 9.9]
    main(values, 5)
