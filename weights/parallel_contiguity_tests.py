import pysal 
from pysal.cg.standalone import get_shared_segments
import numpy as np
from collections import defaultdict
from itertools import combinations
import multiprocessing as mp
import time

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

    results = {}
    results['n_polygons'] = numPoly
    results['potential_neighbors'] = w
    results['shapes'] = shapes
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

def pcheck_joins2(q,resultq, weight_type='ROOK'):
    pid=mp.current_process()._identity[0]
    while True:
        work = q.get()
        if work == None:
            #print "Got the pill."
            q.task_done()
            break
        #Unpack the args from q
        potential_neighbors = work[0]
        shapes = work[1]
        polygon_ids = work[2]
        mdict = {}
        weight_type = weight_type.upper()
        #print "Process {} working on polygons {} - {}.".format(pid, polygon_ids[0], polygon_ids[-1])
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
                if polyId not in mdict:
                    mdict[polyId] = []           
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
                            d = mdict[polyId]
                            d.append(j)
                            mdict[polyId] = d
                            if j not in mdict:
                                mdict[j] = []                       
                            k = mdict[j]
                            k.append(polyId)
                            mdict[j] = k
                            break
        #Put the resultant dict back into the queue and alert that the work is done.           
        resultq.put(mdict)
        q.task_done()
    return

def bbcommon(bb, bbother):
            """
            Checks for overlaps of bounding boxes. First, east-west, then north-south.
            Element 0 is west, element 2 is east, element 1 is north, element 3 is
            south
            All four checks must be false for chflag to be true, meaning the two
            bounding boxes do not overlap.
            """
            chflag = 0
            if not ((bbother[2] < bb[0]) or (bbother[0] > bb[2])):
                if not ((bbother[3] < bb[1]) or (bbother[1] > bb[3])):
                    chflag = 1
            return chflag
        
if __name__ == "__main__":

    fnames = ['1024_lattice.shp', '10000_lattice.shp', '50176_lattice.shp', '100489_lattice.shp', '1000_poly.shp', '10000_poly.shp', '50000_poly.shp', '100000_poly.shp']
    for fname in fnames:
 
        '''Optimized Serial Code'''
        t2 = time.time()
        res = bin_shapefile('TestData/'+fname)
        w = check_joins(res['potential_neighbors'], res['shapes'])
        t3 = time.time()
        print str(fname)
        print 'time refactored prior to parallelization: ', str(t3-t2)   
        
        '''PARALLEL CODE'''
        t6 = time.time()
        res = bin_shapefile('TestData/'+fname)
        cores = mp.cpu_count()
        #Create a joinable queue from which to draw cells and a solution queue to get results
        q = mp.JoinableQueue()
        resultq = mp.Queue()
        
        #Start up a number of child workers equal to the number of cores
        #This is a great way to manage a web service.
        jobs = [mp.Process(target=pcheck_joins2, args=(q, resultq)) for x in range(cores)]
        for job in jobs:
            job.start()
        
        #Sending too small of a job still causes issues with performance.    
        n = len(res['shapes'])
        starts = range(0,n,n/cores*2)
        ends = starts[1:]
        ends.append(n)        
        offsets = [ range(z[0],z[1]) for z in zip(starts, ends) ]
    
        #Load the work into the queue
        #As the jobs are loaded, they start, so we avoid some of the packing overhead.
        for i in offsets:
            args = []
            args.append(res['potential_neighbors'])
            args.append(res['shapes'])
            args.append(i)
            #args.append(weight_type='Queen')
            q.put(args)
        
        #Load a poison pill into the queue to kill the children when work is done
        for i in range(cores):
            q.put(None)    
    
        results = []
        for i in range(len(offsets)):
            results.append(resultq.get())
        t7 = time.time()
        
        ddict = defaultdict(set)
        for d in (results):
            for key, value in d.items():
                for v in value:
                    ddict[key].add(v)
        t8 = time.time()
        
 
        print "MP using Process: {}".format(t8-t6)
        #print "Combining results took {} seconds".format(t8-t7)
        print "Are the results the same? {}".format(ddict == w)
        print
        