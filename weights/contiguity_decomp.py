from IPython.parallel import Client
from itertools import combinations

client = Client()

print client.ids

ncpus = len(client.ids)

import pysal as ps
import numpy as np

sf = ps.open(ps.examples.get_path("nat.shp"))
#sf = ps.open(ps.examples.get_path("columbus.shp"))
#sf = ps.open(ps.examples.get_path("sids2.shp"))
bb = sf.bbox


w = (bb[0] - bb[2]) / (ncpus-1)

w = np.abs(w)

bounds = np.arange(bb[0], bb[2]+w, w)[1:]
bounds[-1] += 0.0001*w
bins = {}
ids = {} 
for b in range(len(bounds)):
    bins[b] = []
    ids[b] = []

shps = []

for i,shp in enumerate(sf):
    shps.append(shp)

    bbi = shp.bounding_box
    left = np.nonzero((bounds > bbi.left))[0][0]
    right = np.nonzero((bounds > bbi.right))[0][0]
    bids = range(left, right+1)
    for bid in bids:
        bins[bid].append(shp)
        ids[bid].append(i)



#sf.close()

def bf_queen(shps):
    n = len(shps)
    w = {}
    coords = {}

    for i in range(n-1):
        si = shps[i]
        vertsi = si.vertices
        for vi in vertsi:
            if vi not in coords:
                coords[vi] = set([i])
            else:
                coords[vi] = coords[vi].union(set([i]))
        for j in range(i+1,n):
            sj = shps[j]
            vertsj = sj.vertices
            for vj in vertsj:
                if vj not in coords:
                    coords[vj] = set([j])
                else:
                    coords[vj] = coords[vj].union(set([j]))
    for coord in coords:
        if len(coords[coord]) > 1:
            pairs = combinations(coords[coord],2)
            for pair in pairs:
                l,r = pair
                if l not in w:
                    w[l] = set([r])
                else:
                    w[l] = w[l].union(set([r]))
                if r not in w:
                    w[r] = set([l])
                else:
                    w[r] = w[r].union(set([l]))
    return w, coords

import time
t1 = time.time()
res, c = bf_queen(shps)
t2 = time.time()

print 'Sequential: ',t2-t1

# parallel
t1 = time.time()
view = client[0:-1]
with client[:].sync_imports():
    from itertools import combinations
results = view.map(bf_queen, bins.values())
t2 = time.time()

#map process ids back to original ids
neighbors = {}
for i,result in enumerate(results.result):
    neigh,c = result
    for key in neigh:
        idi = ids[i][key]
        idjs = [ ids[i][j] for j in neigh[key]]
        if idi not in neighbors:
            neighbors[idi] = set(idjs)
        else:
            neighbors[idi] = neighbors[idi].union(set(idjs))

t3 = time.time()
print 'Parallel: ', t3-t1


