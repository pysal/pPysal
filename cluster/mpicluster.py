from mpi4py import MPI

#MPI Boilerplate
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()
status = MPI.Status()
host = MPI.Get_processor_name()

group = comm.Get_group()


for r in range(ncores):
    if r == rank:
        print "I am core {} of {} on {} in group {}".format(r, ncores, host, group)
