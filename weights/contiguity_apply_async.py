#Contiguity using map_async

import pysal
from binning import bin_shapefile, bbcommon
from collections import defaultdict
import multiprocessing as mp
import time
import sys

def pcheck_joina(polygon_ids, weight_type='ROOK'):

    weight_type = weight_type.upper()
    w = {}
    
    if weight_type == 'ROOK':
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
                w[polyId] = []           
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
                        d = w[polyId]
                        d.append(j)
                        w[polyId] = d
                        if j not in w:
                            w[j] = []                       
                        k = w[j]
                        k.append(polyId)
                        w[j] = k
                        break
        return w
    else:
        print 'unsupported weight type'
        return None     


def async_apply_w_callback(res,cores): 
    
    ddict = defaultdict(set)
    def callback_dict(w):
        for key, value in w.items():
            for v in value:
                ddict[key].add(v)     
    
    t1 = time.time()
    #cores = mp.cpu_count()
    pool = mp.Pool(cores)
    
    #Get offsets
    n = len(res['shapes'])
    starts = range(0,n,n/cores)
    ends = starts[1:]
    ends.append(n)        
    offsets = [ range(z[0],z[1]) for z in zip(starts, ends) ] 

    for offset in offsets:
        pool.apply_async(pcheck_joina, args=(offset,), callback=callback_dict)
    pool.close()
    pool.join()
    
    t2 = time.time()
    
    return (t2 - t1)

if __name__ == "__main__":

    cores = int(sys.argv[1])
    print "This version uses apply_async with a callback function and {} cores.".format(cores)
    
    fnames = ['1024_lattice.shp', '10000_lattice.shp', '50176_lattice.shp', '100489_lattice.shp', '1000_poly.shp', '10000_poly.shp', '50000_poly.shp', '100000_poly.shp']
    
    for fname in fnames:
        res = bin_shapefile('TestData/'+fname)
        global shapes
        global potential_neighbors 
        shapes = res['shapes']
        potential_neighbors = res['potential_neighbors']   
        
        t = async_apply_w_callback(res,cores)
        
        print "{} required {} seconds.".format(fname, t)