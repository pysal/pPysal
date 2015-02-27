import itertools
import sys
import numpy as np
import pysal as ps
from mpi4py import MPI

from globalsort import globalsort



if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    rstate = np.random.RandomState(123456)
    horizontal_pivots = np.linspace(0, 1.0, comm.size + 1)

    #Phase I: Compute Hi bounds and sort the points - this get the points local to the cores as well
    if rank == 0:
        fname = sys.argv[1]
        shpfileobj = ps.open(fname)
        geoms = []
        x = []
        y = []
        for i, poly in enumerate(shpfileobj):
            for j in poly.vertices[:-1]:
                geoms.append(i)
                x.append(j[0])
                y.append(j[1])
        nvertices = len(x)
        pts = np.empty((nvertices, 3))
        pts[:,0] = x
        pts[:,1] = y
        pts[:,2] = geoms
    else:
        nvertices = None

    nvertices = comm.bcast(nvertices)
    shape = (nvertices, 3)
    comm.Barrier()

    if rank == 0:
        local_hi = globalsort(comm, rank, shape, pts=pts, axis='y')
    else:
        local_hi = globalsort(comm, rank, shape, pts=None, axis='y')

    '''
    for i in range(comm.size):
        if rank == i:
            print i, local_hi[local_hi[:,0].argsort()]
            sys.exit()
    '''
    comm.Barrier()

    #Compute coincident x geometries
    local_hi = local_hi[np.lexsort((local_hi[:,0], local_hi[:,1]))]
    coincidentx = []
    seed = local_hi[0][:2]
    clist = set([])
    for i in local_hi:
        if np.array_equal(i[:2], seed):
            clist.add(i[2])
        else:
            coincidentx.append(clist)
            clist = set([i[2]])
            seed = i[:2]

    comm.Barrier()

    if rank == 0:
        for j in coincidentx:
            print j
        #I need the pivots here to do a boundary check.  Then aggregate.
    sys.exit()

    if rank == 0:
        local_vj = globalsort(comm, rank, shape, pts=pts, axis='x')
    else:
        local_vj = globalsort(comm, rank, shape, pts=None, axis='x')

    #Compute coincident y geometries
    coincidenty = []
    seed = local_vj[0][1]
    clist = set([])
    for i in local_vj:
        if i[1] == seed:
            clist.add(i[2])
        else:
            coincidenty.append(clist)
            clist = set([i[2]])
            seed = i[1]

    comm.Barrier()

    for i in range(comm.size):
        if i == rank:
            print i, coincidentx[0], coincidenty[0]
            sys.exit()
    '''
    for i in range(comm.size):
        if rank == i:
            print i, local_vj
    '''

    #Compute the Vj nearest neighbor set
    vj = {}
    for i in local_vj[:,2]:
        vj[i] = (np.inf, -1)
    nearestneighbor(vj, local_vj)

    '''
    for i in range(comm.size):
        if rank == i:
            print i, vj
    '''

    #Compute Iij
    cornercoords = np.array(list(itertools.product(vertical_pivots[rank:rank+2],
                                          horizontal_pivots)))

    p = comm.size
    iij = []
    for i in range(p):
        iij.append([cornercoords [i],
                   cornercoords [i+1],
                   cornercoords [i+p + 1],
                   cornercoords [i+p+2]])

    '''
    for i in range(comm.size):
        if rank == i:
            print i, iij, '\n'
    '''

    comm.Barrier()

    #Compute Cij
    cij = []
    for i in range(p):
        ci = []
        seen = []
        bottom = horizontal_pivots[i]
        top = horizontal_pivots[i+1]
        candidates = local_vj[np.logical_or(local_vj[:,1] < bottom, local_vj[:,1] > top)]

        for points in itertools.product(candidates, iij[i]):
            a, b = points
            dist = euclidean(a[:2], b)
            if dist < vj[a[2]][0]:
                #A pseudo set operation to not get duplicates of unhashable ndarrays
                if a[2] not in seen:
                    ci.append(a)
                seen.append(a[2])
        cij.append(np.array(ci))

    cijsizes = np.zeros(comm.size, dtype=np.int)
    for i, j in enumerate(cij):
        if len(j) > 0:
            cijsizes[i] = j.shape[0] * j.shape[1]
    comm.Barrier()

    #Gather Ci by first collecting the number of pts, and then the pts
    sizes = np.empty(comm.size, dtype=np.int)
    for i in range(comm.size):
        comm.Gather([cijsizes[i], MPI.INTEGER],
                    [sizes, MPI.INTEGER], root=i)

    pointcount = np.sum(sizes)

    candidates = np.empty(pointcount, dtype=np.float64)
    for i in range(comm.size):
        comm.Gatherv([cij[i].ravel(), MPI.DOUBLE],
                     [candidates, (sizes, None), MPI.DOUBLE],
                     root=i)
    comm.Barrier()

    '''
    for i in range(comm.size):
        if rank == i:
            print i, ci.reshape(-1,3), local_hi, '\n'
    '''

    candidates = candidates.reshape(-1, 3)

    ci = {}
    keys = set(candidates[:,2]).union(set(local_hi[:,2]))
    for i in keys:
        ci[i] = (np.inf, -1)

    for points in itertools.product(candidates, local_hi):
        a, b = points
        if a[2] == b[2]:
            continue
        else:
            dist = euclidean(a[:2], b[:2])
            if dist < ci[a[2]][0]:
                ci[a[2]] = (dist, b[2])
            if dist < ci[b[2]][0]:
                ci[b[2]] = (dist, a[2])

    '''
    for i in range(comm.size):
        if rank == i:
            print ci
    '''

    #If hi, vj, and ci are ndarrays to start, just stack (faster)
    keys = set(vj.keys()).union(set(hi.keys()).union(set(ci.keys())))
    keys = sorted(list(keys))
    idx = np.arange(len(keys))
    lookup = {keys[i]:idx[i] for i in range(len(keys))}
    local_distances = np.zeros((len(keys), 3))
    local_distances[:,2] = np.inf

    for k, v in hi.iteritems():
        pos = lookup[k]
        if v[0] < local_distances[pos][2]:
            local_distances[pos] = [k, v[1], v[0]]

    for k, v in vj.iteritems():
        pos = lookup[k]
        if v[0] < local_distances[pos][2]:
            local_distances[pos] = [k, v[1], v[0]]

    for k, v in ci.iteritems():
        pos = lookup[k]
        if v[0] < local_distances[pos][2]:
            local_distances[pos] = [k, v[1], v[0]]

    #Collect distances
    step = shape[0] / comm.size
    id_pivots = range(step, shape[0], step)
    indices = local_distances[:,0].searchsorted(id_pivots)

    partitions = np.split(local_distances, indices, axis=0)
    partitionsizes = np.array([i.shape[0] * shape[1] for i in partitions])
    comm.Barrier()

    #Gather Ci by first collecting the number of pts, and then the pts
    sizes = np.empty(comm.size, dtype=np.int)
    for i in range(comm.size):
        comm.Gather([partitionsizes[i], MPI.INTEGER],
                    [sizes, MPI.INTEGER], root=i)

    comm.Barrier()

    pointcount = np.sum(sizes)
    nearest = np.empty(pointcount, dtype=np.float64)
    for i in range(comm.size):
        comm.Gatherv([partitions[i].ravel(), MPI.DOUBLE],
                     [nearest, (sizes, None), MPI.DOUBLE],
                     root=i)
    comm.Barrier()

    nearest = nearest.reshape(-1, 3)
    nearest = nearest[nearest[:,0].argsort()]

    near = np.empty((shape[0] / comm.size, 3))
    near[:,2] = np.inf

    for i in nearest:
        cidx = i[0] - (rank * (shape[0] / comm.size))
        if i[2] < near[cidx][2]:
            near[cidx] = i
    comm.Barrier()

    for i in range(comm.size):
        if i == rank:
            print i, near
            comm.Barrier()

    sys.exit()
    if rank == 0:
        import pylab as plt
        labels = np.arange(len(pts))
        x, y = np.meshgrid(horizontal_pivots, vertical_pivots)
        plt.plot(x, y)
        plt.plot(y, x)
        plt.plot(x,y, 'ro')
        plt.plot(pts[:,0], pts[:,1], 'o')

        for label, x, y in zip(labels, pts[:, 0], pts[:, 1]):
            plt.annotate(label,
                     xy = (x, y), xytext = (6, 3),
                     textcoords = 'offset points', ha = 'right', va = 'bottom')

        plt.xticks([])
        plt.yticks([])

        plt.xlabel(r'$V_{j}$', fontsize=18)
        plt.ylabel(r'$H_{i}$', fontsize=18)

        plt.show()
