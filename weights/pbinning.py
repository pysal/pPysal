#Binning
import pysal
import multiprocessing as mp
import time
from collections import defaultdict

# delta to get buckets right
DELTA = 0.000001

QUEEN = 1
ROOK = 2

# constants for bucket sizes
BUCK_SM = 8
BUCK_LG = 80
SHP_SMALL = 1000

global shpFileObject


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


def get_bbox(offset,shpRows, shapes, minbox, binWidth):
    #Allocate in the child instead of in the parent
    poly2Column = {}
    poly2Row = {}
    bbcache = {}
    columns = {}
    rows = {}
    for i in range(len(shapes)):
        columns[i] = set()
        rows[i] = set()

    for i in offset:
        poly2Column[i] = set()
        poly2Row[i] = set()

        shpObj = shpRows[i]
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
    return [poly2Column, poly2Row,columns, rows, bbcache ]

#Callback setup here.
poly2Column = {}
poly2Row = {}
columns_cb = defaultdict(set)
rows_cb = defaultdict(set)
bbox_cb = {}

def bbox_callback(r):
    poly2Column.update(r[0])
    poly2Row.update(r[1])
    #columns
    for key, value in r[2].items():
        for v in value:
            columns_cb[key].add(v)      
    #rows
    for key, value in r[3].items():
        for v in value:
            rows_cb[key].add(v)   
    #bbcache
    bbox_cb.update(r[4])
    
def bin_shapefile(shpFile, wtype='rook', n_cols=10, n_rows=10, buff=1.0001):
    
    shpFileObject = pysal.open(shpFile)

    if shpFileObject.type != pysal.cg.Polygon: #This is really slow...
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
    #columns = {}
    #rows = {}
    #for i in range(bucketmin):
        #columns[i] = set()
        #rows[i] = set()
    #rows = [set() for i in range(bucketmin)]

    minbox = shapebox[:2] * \
        2                                  # minx,miny,minx,miny
    binWidth = [lengthx, lengthy] * \
        2                              # lenx,leny,lenx,leny
    #bbcache = {}
    #poly2Column = [set() for i in range(numPoly)]
    #poly2Row = [set() for i in range(numPoly)]
    
    #Parallelize the generation of the bbox here.
    cores = mp.cpu_count()
    pool = mp.Pool(cores) #This defaults to all cores for now, double for Cortez
    

    starts = range(0,numPoly,numPoly/cores)
    ends = starts[1:]
    ends.append(numPoly)    
    offsets = [ range(z[0],z[1]) for z in zip(starts, ends) ]

    #Load the shapefile object into a list
    shpRows = []
    for x in range(numPoly):
        shpRows.append(shpFileObject.get(x))

    for offset in offsets:
        pool.apply_async(get_bbox, args=(offset,shpRows, shapes, minbox, binWidth), callback=bbox_callback)
    pool.close() 
    pool.join() 

    # loop over polygons rather than bins
    #The union calls are what is really slow here.
    w = {}

    #Parallelize the intersection here.
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
