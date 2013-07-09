
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


def bf_contiguity(shpfile, wttype = "QUEEN"):
    f = ps.open(shpfile)
    neighbors = {}
    if wttype.upper() == "QUEEN":
        vertices = {}
        for i, shp in enumerate(f):
            si = set([i])
            for vertex in shp.vertices:
                if vertex not in vertices:
                    vertices[vertex] = si
                else:
                    vertices[vertex] = vertices[vertex].union(si)
        f.close()
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
        edges = {}
        for i, shp in enumerate(f):
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
        f.close()
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
        f.close()
        print "Weight type not supported: ", wttype



if __name__ == '__main__':

    sf = ps.examples.get_path("columbus.shp")
    queen = bf_contiguity(sf, wttype='queen')
    rook = bf_contiguity(sf, wttype="ROOK")
