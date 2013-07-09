
"""
Brute force contiguity builders used for testing results of parallel
algorithms
"""
_author_ = "Serge Rey <sjsrey@gmail.com>"


import pysal as ps
import numpy as np
from itertools import combinations

#sf = ps.open(ps.examples.get_path("nat.shp"))
#sf = ps.open(ps.examples.get_path("columbus.shp"))
#sf = ps.open(ps.examples.get_path("sids2.shp"))
#bb = sf.bbox


def bf_contiguity(shps, wttype = "QUEEN"):
    """
    Brute force contiguity builder

    Arguments
    ---------

    shps: list of pysal.cg.Polygon shapes

    wttype: string
            contiguity type

    Returns
    -------
    neighbors: dict
               key is id, value is list of neighboring ids

    """
    neighbors = {}
    if wttype.upper() == "QUEEN":
        vertices = {}
        for i, shp in enumerate(shps):
            si = set([i])
            for vertex in shp.vertices:
                if vertex not in vertices:
                    vertices[vertex] = set()
                vertices[vertex] = vertices[vertex].union(si)
        for vertex in vertices:
            pairs = combinations(vertices[vertex], 2)
            for l,r in pairs:
                if l not in neighbors:
                    neighbors[l] = set()
                if r not in neighbors:
                    neighbors[r] = set()
                neighbors[l] = neighbors[l].union([r])
                neighbors[r] = neighbors[r].union([l])
        return neighbors
    elif wttype.upper() == 'ROOK':
        vertices = {}
        for i, shp in enumerate(shps):
            si = set([i])
            for vertex in shp.vertices:
                if vertex not in vertices:
                    vertices[vertex] = set()
                vertices[vertex] = vertices[vertex].union(si)
        checked = []
        edges = {}
        for vertex in vertices:
            pairs = combinations(vertices[vertex], 2)
            for l,r in pairs:
                if l not in neighbors:
                    neighbors[l] = set()
                if r not in neighbors:
                    neighbors[r] = set()
                lr = set([l,r])
                if checked.count(lr) == 0:
                    checked.append(lr)
                    # check for shared edge
                    vl = shps[l].vertices
                    vr = shps[r].vertices
                    vls = set(vl)
                    vrs = set(vr)
                    common = vls.intersection(vrs)
                    if len(common) >= 2:
                        # at least two vertices in common
                        # now check for sequential positioning
                        good_left = []
                        for el in combinations(common, 2):
                            # check left poly for an edge defined 
                            o,d = el
                            # how many times does point appear
                            co = vl.count(o) 
                            cd = vl.count(d)
                            sc = 0
                            for po in range(co):
                                lo = vl.index(o, sc)
                                sd = 0
                                for pd in range(cd):
                                    ld = vl.index(d, sd)
                                    diff = np.abs(lo - ld)
                                    sd = ld + 1
                                    if diff == 1:
                                        good_left.append(el)
                                sc = lo + 1
                        for el in good_left:
                            o,d = el
                            co = vr.count(o)
                            cd = vr.count(d)
                            sc = 0
                            for po in range(co):
                                ro = vr.index(o, sc)
                                sd = 0
                                for pd in range(cd):
                                    rd = vr.index(d, sd)
                                    diff = np.abs(ro - rd)
                                    sd = rd + 1
                                    if diff == 1:
                                        neighbors[l] = neighbors[l].union(set([r]))
                                        neighbors[r] = neighbors[r].union(set([l]))
                                sc = ro + 1
        return neighbors
    else:
        print "Weight type not supported: ", wttype


def qf_shapefile(sf):
    shps = []
    f = ps.open(sf)
    for shp in f:
        shps.append(shp)
    f.close()
    return  bf_contiguity(shps, wttype = 'QUEEN')

def rf_shapefile(sf):
    shps = []
    f = ps.open(sf)
    for shp in f:
        shps.append(shp)
    f.close()
    return  bf_contiguity(shps, wttype = 'ROOK')


if __name__ == '__main__':

    sf = ps.examples.get_path("columbus.shp")
    queen_col = qf_shapefile(sf)
    rook_col = rf_shapefile(sf)
    wrc = ps.W(rook_col)
    print wrc.histogram

    import time
    sf = ps.examples.get_path("nat.shp")
    t1 = time.time()
    queen = qf_shapefile(sf)
    t2 = time.time()
    print "National queen: ", t2-t1
    sf = ps.examples.get_path("nat.shp")
    t1 = time.time()
    rook = rf_shapefile(sf)
    t2 = time.time()
    print "National rook: ", t2-t1

