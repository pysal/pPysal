import pysal
from pysal.cg.standalone import get_shared_segments
import numpy as np
from collections import defaultdict
from itertools import combinations



# delta to get buckets right
DELTA = 0.000001

QUEEN = 1
ROOK = 2

# constants for bucket sizes
BUCK_SM = 8
BUCK_LG = 80
SHP_SMALL = 1000


def bbcommon(bb, bbother):
    """
    Checks for overlaps of bounding boxes. First, east-west, then north-south.
    Element 0 is west, element 2 is east, element 1 is north?, element 3 is
    south?
    All four checks must be false for chflag to be true, meaning the two
    bounding boxes do not overlap.
    """
    chflag = 0
    if not ((bbother[2] < bb[0]) or (bbother[0] > bb[2])):
        if not ((bbother[3] < bb[1]) or (bbother[1] > bb[3])):
            chflag = 1
    return chflag


def bin_shapefile(shpFile, wtype='rook', n_cols=10, n_rows=10, buff=1.0001):

    # set up bins

    shpfr = pysal.open(shpFile)

    n_polygons = len(shpfr)


    c_bins = np.linspace(shpfr.bbox[0], shpfr.bbox[2]*buff, n_cols)
    r_bins = np.linspace(shpfr.bbox[1], shpfr.bbox[3]*buff, n_rows)

    cells = dict([((i,j), set()) for i in range(n_rows+1) \
            for j in range(n_cols+1)])
    poly2cells= dict([ (i,set()) for i in xrange(n_polygons) ])


    shps = []

    for i in xrange(n_polygons):
        shpObj = shpfr.get(i)
        shps.append(shpObj)
        sbbox = shpObj.bounding_box[:]
        l = np.digitize([sbbox[0]], c_bins)[0]
        r = np.digitize([sbbox[2]], c_bins)[0] + 1
        b = np.digitize([sbbox[1]], r_bins)[0]
        t = np.digitize([sbbox[3]], r_bins)[0] + 1
        for row in range(l,r):
            for col in range(b,t):
                cells[(row,col)].add(i)
                poly2cells[i].add((row,col))
    shpfr.close()
    results = {}
    results['n_polygons'] = n_polygons
    results['cells'] = cells
    results['poly2cells'] = poly2cells
    results['shapes'] = shps
    return results


def check_joins(poly2cells, cells, shapes, weight_type='ROOK',
        polygon_ids = []):
    w = {}
    weight_type = weight_type.upper()
    if polygon_ids == []:
        polygon_ids = range(len(poly2cells))
    for poly_id in polygon_ids:
        poly_cells = poly2cells[poly_id]
        candidates = set()
        for cell in poly_cells:
            candidates = candidates.union(cells[cell])
        candidates.discard(poly_id)
        # only check half-set
        candidates = [ poly for poly in candidates if poly > poly_id ] 
        p0 = shapes[poly_id]
        verts0 = set(p0.vertices)
        if poly_id not in w:
            w[poly_id] = set()
        for id1 in candidates:
            join = False
            p1 = shapes[id1]
            if bbcommon(p0.bounding_box,p1.bounding_box):
                common = verts0.intersection(p1.vertices)
                n_common = len(common)
                if n_common > 1 and weight_type == 'ROOK':
                    d0 = defaultdict(list)
                    d1 = defaultdict(list)
                    for i,v in enumerate(p0.vertices[:-1]):
                        d0[v].append(p0.vertices[i+1])
                        d0[p0.vertices[i+1]].append(v)
                    for i,v in enumerate(p1.vertices[:-1]):
                        d1[v].append(p1.vertices[i+1])
                        d1[p1.vertices[i+1]].append(v)
                    for seq in combinations(common,2):
                        l,r = seq
                        if l in d0[r] and l in d1[r]:
                            join = True
                            break
                if n_common > 0 and weight_type == 'QUEEN':
                    join = True
                if join:
                    w[poly_id].add(id1)
                    if id1 not in w:
                        w[id1] = set()
                    w[id1].add(poly_id)
    return w



if __name__ == "__main__":
    import time
    #fname = pysal.examples.get_path('10740.shp')
    fname = pysal.examples.get_path('nat.shp')
    t0 = time.time()
    c = pysal.weights.Contiguity.ContiguityWeights(pysal.open(fname), ROOK)
    t1 = time.time()
    print "using " + str(fname)
    print "time elapsed for ... using bins: " + str(t1 - t0)
    res= bin_shapefile(fname)



    t2 = time.time()
    res= bin_shapefile(fname)
    w = check_joins(res['poly2cells'], res['cells'], res['shapes'])
    t3 = time.time()
    print 'time refactored prior to parallelization: ', str(t3-t2)

    print c.w == w

    keys = c.w.keys()
    for key in keys:
        if c.w[key] != w[key]:
            print key, c.w[key], w[key]

