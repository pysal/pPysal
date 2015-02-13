import sys
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

shape = (40,2)
samplesize = 4

#Phase I: Initialize
if rank == 0:
    #Generate an array of points
    pts = np.random.random(shape)
else:
    pts = None

local_pts = np.empty((shape[0] / comm.size, 2))

#Phase II: Scatter the data, local quicksort, return regular sample
comm.Scatter([pts, MPI.DOUBLE], [local_pts, MPI.DOUBLE])
local_pts = np.sort(local_pts, axis=0)
samples = local_pts[::samplesize]
comm.Barrier()

#Phase III: Gather, merge, and broadcast samples
sampledpts = comm.gather(samples, root=0)  #Replace with Gather
if rank == 0:
    sampledpts = np.sort(np.concatenate(sampledpts), axis=0)  #Should be a multimerge
    pivots = sampledpts[::3][1:]
    print pivots
else:
    pivots = np.empty((6), dtype=np.float64)
comm.Barrier()
#Ndarrays need to be 'flat' to communicate
comm.Bcast([pivots.ravel(), MPI.DOUBLE])

pivots = pivots.reshape(-1,2)

#Phase IV: Partition local data into p classes using pivots
indices = local_pts[:,1].searchsorted(pivots[:,1]) #Hard coded for y-axis
partitions = np.split(local_pts, indices, axis=0)
partitionsizes = np.array([i.shape[0] for i in partitions])
comm.Barrier()

#Phase V: Gather the ith partition onto this core
sizes = np.empty(comm.size, dtype=np.int)
for i in range(comm.size):
    comm.Gather([partitionsizes[i], MPI.INTEGER],
                [sizes, MPI.INTEGER], root=i)

pointcount = np.sum(sizes) * 2  # 2 vertices
#Create the memory space to collect the ith lists
data = np.empty(pointcount, dtype=np.float64)
print rank, data.size
for i in range(comm.size):
    trans = partitions[i].ravel()
    print trans.shape
    comm.Gather([trans, MPI.DOUBLE],
                [data, MPI.DOUBLE], root=i)

#print rank, len(partitions)
