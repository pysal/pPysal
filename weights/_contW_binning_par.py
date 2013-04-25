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

    t1 = time.time()

    shpFileObject = pysal.open(shpFile)

    if shpFileObject.type != pysal.cg.Polygon:
        return False

    shapebox = shpFileObject.bbox      # bounding box

    numPoly = len(shpFileObject)
    shapes = [[]] * numPoly

    # bucket size
    if (numPoly < SHP_SMALL):
        bucketmin = numPoly / BUCK_SM + 2
    else:
        bucketmin = numPoly / BUCK_LG + 2
    # bucket length
    lengthx = ((shapebox[2] + DELTA) - shapebox[0]) / bucketmin
    lengthy = ((shapebox[3] + DELTA) - shapebox[1]) / bucketmin

    # initialize buckets
    columns = [set() for i in range(bucketmin)]
    rows = [set() for i in range(bucketmin)]

    minbox = shapebox[:2] * \
        2                                  # minx,miny,minx,miny
    binWidth = [lengthx, lengthy] * \
        2                              # lenx,leny,lenx,leny
    bbcache = {}
    poly2Column = [set() for i in range(numPoly)]
    poly2Row = [set() for i in range(numPoly)]
    for i in range(numPoly):
        shpObj = shpFileObject.get(i)
        bbcache[i] = shpObj.bounding_box[:]
        shapes[i] = shpObj
        projBBox = [int((shpObj.bounding_box[:][j] -
                         minbox[j]) / binWidth[j]) for j in xrange(4)]
        for j in range(projBBox[0], projBBox[2] + 1):
            columns[j].add(i)
            poly2Column[i].add(j)
        for j in range(projBBox[1], projBBox[3] + 1):
            rows[j].add(i)
            poly2Row[i].add(j)
    # loop over polygons rather than bins
    w = {}
    for polyId in xrange(numPoly):
        idRows = poly2Row[polyId]
        idCols = poly2Column[polyId]
        rowPotentialNeighbors = set()
        colPotentialNeighbors = set()
        for row in idRows:
            rowPotentialNeighbors = rowPotentialNeighbors.union(rows[row])
        for col in idCols:
            colPotentialNeighbors = colPotentialNeighbors.union(
                columns[col])
        potentialNeighbors = rowPotentialNeighbors.intersection(
            colPotentialNeighbors)
        if polyId not in w:
            w[polyId] = set()
        for j in potentialNeighbors:
            if polyId < j:
                if bbcommon(bbcache[polyId], bbcache[j]):
                    w[polyId].add(j)
    t2 = time.time()
    print t2-t1

    results = {}
    results['n_polygons'] = numPoly
    results['potential_neighbors'] = w
    results['shapes'] = shapes
    t2 = time.time()
    print t2 - t1
    return results


def check_joins(potential_neighbors, shapes, weight_type='ROOK',
        polygon_ids = []):
    w = {}
    weight_type = weight_type.upper()

    if not polygon_ids:
        polygon_ids = xrange(len(shapes))

    if weight_type == 'QUEEN':
        # check for a shared vertex
        vertCache = {}
        for polyId in polygon_ids:
            iVerts = shapes[polyId].vertices
            nbrs = potential_neighbors[polyId]
            if polyId not in vertCache:
                vertCache[polyId] = set(iVerts)
            if polyId not in w:
                w[polyId] = set()
            for j in nbrs:
                join = False
                if j not in vertCache:
                    vertCache[j] = set(shapes[j].vertices)
                common = vertCache[polyId].intersection(vertCache[j])
                if len(common) > 0:
                    join = True
                if join:
                    w[polyId].add(j)
                    if j not in w:
                        w[j] = set()
                    w[j].add(polyId)
        return w
    elif weight_type == 'ROOK':
        # check for a shared edge
        edgeCache = {}
        for polyId in polygon_ids:
            if polyId not in edgeCache:
                iEdges ={}
                iVerts = shapes[polyId].vertices
                nv = len(iVerts)
                ne = nv - 1
                for i in range(ne):
                    l = iVerts[i]
                    r = iVerts[i+1]
                    iEdges[(l,r)] = []
                    iEdges[(r,l)] = []
                edgeCache[polyId] = iEdges
            nbrs = potential_neighbors[polyId]
            if polyId not in w:
                w[polyId] = set()
            for j in nbrs:
                join = False
                if j not in edgeCache:
                    jVerts = shapes[j].vertices
                    jEdges = {}
                    nv = len(jVerts)
                    ne = nv - 1
                    for e in range(ne):
                        l = jVerts[e]
                        r = jVerts[e+1]
                        jEdges[(l,r)] = []
                        jEdges[(r,l)] = []
                    edgeCache[j] = jEdges
                for edge in edgeCache[j]:
                    if edge in edgeCache[polyId]:
                        join = True
                        w[polyId].add(j)
                        if j not in w:
                            w[j] = set()
                        w[j].add(polyId)
                        break
        return w
    else:
        print 'unsupported weight type'
        return None



def bin_shapefile1(shpFile, wtype='rook', n_cols=10, n_rows=10, buff=1.0001):
        shpFileObject = pysal.open(shpFile)

        if shpFileObject.type != pysal.cg.Polygon:
            return False

        shapebox = shpFileObject.bbox      # bounding box

        numPoly = len(shpFileObject)

        # bucket size
        if (numPoly < SHP_SMALL):
            bucketmin = numPoly / BUCK_SM + 2
        else:
            bucketmin = numPoly / BUCK_LG + 2
        # bucket length
        lengthx = ((shapebox[2] + DELTA) - shapebox[0]) / bucketmin
        lengthy = ((shapebox[3] + DELTA) - shapebox[1]) / bucketmin

        # initialize buckets
        columns = [set() for i in range(bucketmin)]
        rows = [set() for i in range(bucketmin)]

        minbox = shapebox[:2] * \
            2                                  # minx,miny,minx,miny
        binWidth = [lengthx, lengthy] * \
            2                              # lenx,leny,lenx,leny
        bbcache = {}
        poly2Column = [set() for i in range(numPoly)]
        poly2Row = [set() for i in range(numPoly)]
        for i in range(numPoly):
            shpObj = shpFileObject.get(i)
            bbcache[i] = shpObj.bounding_box[:]
            projBBox = [int((shpObj.bounding_box[:][j] -
                             minbox[j]) / binWidth[j]) for j in xrange(4)]
            for j in range(projBBox[0], projBBox[2] + 1):
                columns[j].add(i)
                poly2Column[i].add(j)
            for j in range(projBBox[1], projBBox[3] + 1):
                rows[j].add(i)
                poly2Row[i].add(j)
        # loop over polygons rather than bins
        w = {}
        for polyId in xrange(numPoly):
            idRows = poly2Row[polyId]
            idCols = poly2Column[polyId]
            rowPotentialNeighbors = set()
            colPotentialNeighbors = set()
            for row in idRows:
                rowPotentialNeighbors = rowPotentialNeighbors.union(rows[row])
            for col in idCols:
                colPotentialNeighbors = colPotentialNeighbors.union(
                    columns[col])
            potentialNeighbors = rowPotentialNeighbors.intersection(
                colPotentialNeighbors)
            if polyId not in w:
                w[polyId] = set()
            for j in potentialNeighbors:
                if polyId < j:
                    if bbcommon(bbcache[polyId], bbcache[j]):
                        w[polyId].add(j)
        shpFileObject.close()
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
    w = check_joins(res['potential_neighbors'], res['shapes'])
    t3 = time.time()
    print 'time refactored prior to parallelization: ', str(t3-t2)

    print c.w == w

    #keys = c.w.keys()
    #for key in keys:
    #    if c.w[key] != w[key]:
    #        print key, c.w[key], w[key]

