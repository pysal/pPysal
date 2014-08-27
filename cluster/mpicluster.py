from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()
status = MPI.Status()
name = MPI.Get_processor_name()

for r in range(ncores):
    if r == rank:
        print "I am core {} of {} on {}".format(r, ncores, name)
