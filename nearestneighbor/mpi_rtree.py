
#Standard Lib. Imports
from copy_reg import pickle
import itertools
import multiprocessing as mp
import sys
import time
from types import MethodType
from mpi4py import MPI

import numpy as np

import pysal as ps
import pysal.cg.rtree as rtree
from pysal.cg.standalone import get_shared_segments

QUEEN = 1
ROOK = 2
class _PolyQ(dict):
    def __init__(self):
        dict.__init__(self)
        self.size = 20  # use the first 20 objects to calculate the average Size.
        self.ids = []

    def __checkSize(self):
        """
        Use the objects in the Q to calculate the average size of the objects
        Adjust Q.size to hold Q_TARGET_MEM_SIZE/avgSize object
        This is as many average size object that fit into Q_TARGET_MEM_SIZE
        """
        if len(self.ids) > 50:
            return True
        return False

    def add(self, poly):
        if poly.id not in self:
            if len(self.ids) >= self.size:
                if self.__checkSize():
                    del self[self.ids.pop(0)]
            self[poly.id] = poly
            self.ids.append(poly.id)

class ContiguityWeights_rtree:
    def __init__(self, geoObj, joinType=ROOK):
        self.index = rtree.Rtree()
        self.geoObj = geoObj
        self.joinType = joinType
        self.w = {}
        self.Q = _PolyQ()
        self.cache_hits = 0
        self.cache_misses = 0
        self.ids = set([])  #Hacks this this is why you don't override append...
        self.polys = []

        self.create()
        print self.polys
        #print "Misses: ",self.cache_misses
        #print "Hits: ",self.cache_hits

    def create(self):
        for id, poly in enumerate(self.geoObj):
            self.ids.add(id)
            poly.id = id
            self.append2(poly)
            self.polys.append(poly)

        self.geoObj.close()

    def append2(self, poly):
        self.Q.add(poly)
        b = poly.bounding_box
        bbox = [b.left, b.lower, b.right, b.upper]
        for id in self.index.intersection(bbox): #Queries the rtree
            id = int(id)
            if self.check(id, poly) >= self.joinType:
                self.setW(id, poly.id)
        if poly.id not in self.w:  # add the null cases
            self.w[poly.id] = set()
        self.index.add(poly.id, bbox)

    def setW(self, id0, id1):
        "updates the W matrix seting two polygon's as neighbors"
        w = self.w
        if id0 not in w:
            w[id0] = set()
        if id1 not in w:
            w[id1] = set()
        w[id0].add(id1)
        w[id1].add(id0)

    def check(self, id0, poly1):
        "Check's if two polygon's are neighbors"
        if id0 in self.Q:
            self.cache_hits += 1
            poly0 = self.Q[id0]
        else:
            self.cache_misses += 1
            poly0 = self.geoObj.get(id0)
            poly0.id = id0
            self.Q.add(poly0)
        common = set(poly0.vertices).intersection(set(poly1.vertices))
        if len(common) > 1 and self.joinType == ROOK:
            #double check rook
            if get_shared_segments(poly0, poly1, True):
                return ROOK
            return False
            #for vert in common:
            #    idx = poly0.vertices.index(vert)
            #    IDX = poly1.vertices.index(vert)
            #    try:
            #        if poly0.vertices[idx+1] == poly1.vertices[IDX+1] or poly0.vertices[idx+1] == poly1.vertices[IDX-1]\
            #        or poly0.vertices[idx-1] == poly1.vertices[IDX+1] or poly0.vertices[idx-1] == poly1.vertices[IDX-1]:
            #            return ROOK
            #    except IndexError:
            #        pass
            #return False
        elif len(common) > 0:
            return QUEEN
        else:
            return False

#MPI Boilerplate
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


if rank == 0:
    t1 = time.time()
    shp = ps.open(sys.argv[1])
    t2 = time.time()
    print "File I/O took {} seconds".format(t2 - t1)
    tree_class = ContiguityWeights_rtree(shp, joinType=QUEEN)
    t3 = time.time()
    print "Generating tree took {} seconds.".format(t3 - t2)
    tree = MPI._p_pickle.dumps(tree_class)
else:
    tree = None

tree_pickle = comm.bcast(tree, root=0)
tree = MPI._p_pickle.loads(tree_pickle)

comm.Barrier()
if rank == 0:
    t4 = time.time()
    print "Communication of the rTRee took {} seconds".format(t4 - t3)
for r in range(comm.size):
    if rank == r:
        print tree.Q.ids

sys.exit()
#Compute the offsets in the kdtree.data structure
quotient, remainder = divmod(npoints, comm.size)
scattersize = list([(quotient + remainder)  ]) +\
              [quotient for i in range(comm.size - 1)]
scatteroffsets = [0] + (np.cumsum(scattersize)[:-1].tolist())
comm.Barrier()

'''
for r in range(comm.size):
    if r == rank:
        print scatteroffsets
'''

start = scatteroffsets[rank]
if rank == comm.size - 1:
    stop = None
else:
    stop = scatteroffsets[rank + 1]

local_pts = kdt.data[start:stop]

'''
for r in range(comm.size):
    if rank == r:
        print rank, local_pts.shape
'''
ni, di = kdt.query(local_pts, k=2)
nn_result = np.column_stack((ni, di[:,1]))

'''
for r in range(comm.size):
    if rank == r:
        print nn_result
'''
if rank == 0:
    t6 = time.time()
    print "KDQuery took {} seconds".format(t6 - t5)
    print "Total runtime was {} seconds".format(t6 - t1)
