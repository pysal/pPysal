
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
                    vertices[vertex] = si
                else:
                    vertices[vertex] = vertices[vertex].union(si)
        for vertex in vertices:
            pairs = combinations(vertices[vertex], 2)
            for l,r in pairs:
                if l not in neighbors:
                    neighbors[l] = set([r])
                if r not in neighbors:
                    neighbors[r] = set([l])
                neighbors[l] = neighbors[l].union([r])
                neighbors[r] = neighbors[r].union([l])
        return neighbors

    elif wttype.upper() == "ROOK":
        # check for shared (non-directed) edges
        edges = {}
        for i, shp in enumerate(shps):
            nvi = len(shp.vertices)
            si = set([i])
            for oi in range(nvi-1):
                for di in range(oi+1, nvi):
                    edge = shp.vertices[oi], shp.vertices[di]
                    if edge not in edges:
                        edges[edge] = si
                    twin = shp.vertices[di], shp.vertices[oi]
                    edges[edge] = edges[edge].union(si)
                    if twin not in edges:
                        edges[twin] = si
                    edges[twin] = edges[twin].union(si)
        neighbors = {}
        for edge in edges:
            pairs = combinations(edges[edge], 2)
            for l,r in pairs:
                if l not in neighbors:
                    neighbors[l] = set([r])
                if r not in neighbors:
                    neighbors[r] = set([l])
                neighbors[l] = neighbors[l].union([r])
                neighbors[r] = neighbors[r].union([l])
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
    queen = qf_shapefile(sf)
    rook = rf_shapefile(sf)

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
