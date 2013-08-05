# function that bins and will be mapped in parallel

def equalPysal(wps, wother):
    """
    Test if two w objects are the same
    """
    nDiff = 0
    for key in wps.neighbors:
        if set(wps.neighbors[key]) != wother.neighbors[key]:
            nDiff += 1
    if not nDiff:
        return True
    else:
        return  False
        

def binShapes(shapes, bBoxes, minBox, binWidth, bucketMin, ids = []):
    if not ids:
        ids = range(len(shapes))
    nShapes = len(ids)
    poly2Columns = dict([ (ids[i],set()) for i in range(nShapes) ])
    poly2Rows = dict([ (ids[i],set()) for i in range(nShapes) ])
    columns = dict([ (i,set()) for i in range(bucketMin + 2) ])
    rows = dict([ (i,set()) for i in range(bucketMin + 2)])
    
    for i in range(nShapes):
        idI = ids[i]
        shapeI = shapes[i]
        bBoxI = bBoxes[i][:]
        projBox = [int((bBoxI[j] - minBox[j]) / binWidth[j]) for j in range(4)]
        #print i, projBox, bucketMin
        for j in range(projBox[0], projBox[2] + 1  ):
            columns[j].add(idI)
            poly2Columns[idI].add(j)
        for j in range(projBox[1], projBox[3] + 1):
            rows[j].add(idI)
            poly2Rows[idI].add(j)
    results = {}
    results['poly2Columns'] = poly2Columns
    results['poly2Rows'] = poly2Rows
    results['columns' ] = columns
    results['rows' ] = rows
    return results


def check_contiguity(shapes, bBoxes, poly2Rows, poly2Columns, rows, columns, ids = []):
    # note that ids is a list of particular shapes to check but all other arguments are not sliced
    
    
    def bbcommon(bb, bbother):
        """
        Checks for overlaps of bounding boxes. First, east-west, then north-south.
        Element 0 is west, element 2 is east, element 3 is north, element 1 is
        south.
        All four checks must be false for chflag to be true, meaning the two
        bounding boxes do overlap.
        """
        chflag = 0
        if not ((bbother[2] < bb[0]) or (bbother[0] > bb[2])):
            if not ((bbother[3] < bb[1]) or (bbother[1] > bb[3])):
                chflag = 1
        return chflag

    def queen(shapeA, shapeB):
        """
        Check if shapeA and shapeB are queen neighbors
    
        Arguments
        =========
    
        shapeA: pysal polygon object
    
        shapeB: pysal polygon object
    
        Returns
        =======
    
        1 if true, 0 if false
    
    
        Examples
        ========
    
        >>> sf = pysal.open(pysal.examples.get_path("columbus.shp"))
        >>> p0 = sf.get(0)
        >>> p1 = sf.get(1)
        >>> p3 = sf.get(3)
        >>> import binning
        >>> binning.queen(p0,p1)
        1
        >>> binning.queen(p0,p3)
        0
        >>> binning.queen(p1,p3)
        1
    
    
        """
        if bbcommon(shapeA.bounding_box, shapeB.bounding_box):
            a = set(shapeA.vertices)
            b = set(shapeB.vertices)
            common = a.intersection(b)
            if len(common) > 0:
                return 1
            else:
                return 0
        else:
            return 0



    if not ids:
        ids = poly2Rows.keys()
    
    neighbors = {}
    
    for polyId in ids:
        idRows = poly2Rows[polyId]
        idColumns = poly2Columns[polyId]
        rNeighbors = set()
        cNeighbors = set()
        for row in idRows:
            rNeighbors = rNeighbors.union(rows[row])
        for column in idColumns:
            cNeighbors = cNeighbors.union(columns[column])
        neighborsPolyId = rNeighbors.intersection(cNeighbors)
        if polyId not in neighbors:
            neighbors[polyId] = set()
        for j in neighborsPolyId:
            if polyId < j:
                if queen(shapes[polyId], shapes[j]):
                    if j not in neighbors:
                        neighbors[j] = set()
                    neighbors[j].add(polyId)
                    neighbors[polyId].add(j)
    return neighbors

if __name__ == '__main__':

    import pysal as ps
    import numpy as np
    import multiprocessing as mp
    import time

    def p2(sfname):

        """
        First stage parallel, second stage parallel
        """
        sf = ps.examples.get_path(sfname)
        shpFileObject = ps.open(sf)
        shapeBox = shpFileObject.bbox
        nShapes = len(shpFileObject)
        shapes = []
        bBoxes = []
        for shape in shpFileObject:
            shapes.append(shape)
            bBoxes.append(shape.bounding_box[:])

        # figure out grid that will be used by all processes
        DELTA = 0.000001
        bucketMin = nShapes / 80 + 2
        lengthX = ((shapeBox[2] + DELTA) - shapeBox[0]) / bucketMin
        lengthY = ((shapeBox[3] + DELTA) - shapeBox[1]) / bucketMin
        minBox = shapeBox[:2] * 2 # [minx, miny, minx, miny]
        binWidth = [ lengthX, lengthY] * 2 # [lenX, lenY, lenX, lenY]

        # parallel

        cores = mp.cpu_count()
        pool = mp.Pool(cores)

        step = nShapes / (cores - 1)
        start = range(0, nShapes, step)[0:-1]
        end = start[1:]
        end.append(nShapes)
        slices = zip(start, end)
        stage1 = {}
        for c in range(cores-1):
            pids = range(slices[c][0], slices[c][1])
            stage1[c] = pool.apply_async(binShapes, args=(shapes, bBoxes, minBox,
                binWidth, bucketMin, pids))
        pool.close()
        pool.join()

        results1 = {}
        for c in stage1:
            results1[c] = stage1[c].get()

        # combine results from stage 1

        columns = dict([ (i,set()) for i in range(bucketMin + 2) ])
        rows = dict([ (i,set()) for i in range(bucketMin + 2)])
        poly2Rows = {}
        poly2Columns = {}
        for c in  results1:
            result = results1[c]
            for row in rows:
                rows[row] = rows[row].union(result['rows'][row])
            for column in columns:
                columns[column] = columns[column].union(result['columns'][column])
            for polyId in result['poly2Columns']:
                poly2Columns[polyId] = result['poly2Columns'][polyId]
                poly2Rows[polyId] = result['poly2Rows'][polyId]
        
        # stage two
        pool = mp.Pool(cores)
        cuts = np.cumsum(np.arange(nShapes-1, 0, -1))
        cuts = cuts / (cuts[-1] / (cores - 1))
        start = [ np.nonzero(cuts==c)[0][0] for c in range((cores-1)) ]
        end = start[1:]
        end.append(nShapes)
        slices = zip(start,end)

        
        stage2 = {}
        for c in range(cores - 1):
            pids = range(slices[c][0], slices[c][1])
            stage2[c] = pool.apply_async(check_contiguity, args=(shapes, bBoxes,
                poly2Rows, poly2Columns, rows, columns, pids))
        pool.close()
        pool.join()
        results = {}
        for c in stage2:
            results[c] = stage2[c].get()


        neighbors = dict([(i,set()) for i in range(nShapes) ])

        for c in results:
            for key in results[c]:
                neighbors[key] = neighbors[key].union(results[c][key])

        return ps.W(neighbors,silent_island_warning=True)


    def s1p2(sfname):
        """
        First stage sequential, second stage parallel
        """

        sf = ps.examples.get_path(sfname)
        shpFileObject = ps.open(sf)
        shapeBox = shpFileObject.bbox
        nShapes = len(shpFileObject)
        shapes = []
        bBoxes = []
        for shape in shpFileObject:
            shapes.append(shape)
            bBoxes.append(shape.bounding_box[:])

        # figure out grid that will be used by all processes
        DELTA = 0.000001
        bucketMin = nShapes / 80 + 2
        lengthX = ((shapeBox[2] + DELTA) - shapeBox[0]) / bucketMin
        lengthY = ((shapeBox[3] + DELTA) - shapeBox[1]) / bucketMin
        minBox = shapeBox[:2] * 2 # [minx, miny, minx, miny]
        binWidth = [ lengthX, lengthY] * 2 # [lenX, lenY, lenX, lenY]


        stage1Seq = binShapes(shapes, bBoxes, minBox, binWidth, bucketMin)

        cores = mp.cpu_count()
        pool = mp.Pool(cores)
        cuts = np.cumsum(np.arange(nShapes-1, 0, -1))
        cuts = cuts / (cuts[-1] / (cores - 1))
        start = [ np.nonzero(cuts==c)[0][0] for c in range((cores-1)) ]
        end = start[1:]
        end.append(nShapes)
        slices = zip(start,end)

        rows = stage1Seq['rows']
        columns = stage1Seq['columns']
        poly2Rows = stage1Seq['poly2Rows']
        poly2Columns = stage1Seq['poly2Columns']

        
        stage2 = {}
        for c in range(cores - 1):
            pids = range(slices[c][0], slices[c][1])
            stage2[c] = pool.apply_async(check_contiguity, args=(shapes, bBoxes,
                poly2Rows, poly2Columns, rows, columns, pids))
        pool.close()
        pool.join()
        results = {}
        for c in stage2:
            results[c] = stage2[c].get()


        neighbors = dict([(i,set()) for i in range(nShapes) ])

        for c in results:
            for key in results[c]:
                neighbors[key] = neighbors[key].union(results[c][key])

        return ps.W(neighbors, silent_island_warning=True)


    t1 = time.time()
    sf = ps.examples.get_path('columbus.shp')
    t2 = time.time()
    print 'pysal time: ', t2-t1
    wps = ps.queen_from_shapefile(sf)
    t1 = time.time()
    wp2 = p2(sf)
    t2 = time.time()
    print 'columbus p2: ', t2 - t1
    print equalPysal(wps, wp2)
    t1 = time.time()
    ws1p2 = s1p2(sf)
    t2 = time.time()
    print 'columbus s1p2', t2 - t1
    print equalPysal(wps, ws1p2)


    sf = ps.examples.get_path('nat.shp')
    t1 = time.time()
    wps = ps.queen_from_shapefile(sf)
    t2 = time.time()

    print 'pysal time: ', t2-t1
    t1 = time.time()
    wp2 = p2('nat.shp')
    t2 = time.time()
    print 'nat p2: ', t2 - t1
    print equalPysal(wps, wp2)

    t1 = time.time()
    ws1p2 = s1p2('nat.shp')
    t2 = time.time()
    print 'nat s1p2', t2 - t1


    print equalPysal(wps, ws1p2)



   
    """
    # sequential step 1, parallel step 2

    stage1Seq = binShapes(shapes, bBoxes, minBox, binWidth, bucketMin)

    pool = mp.Pool(cores)
    cuts = np.cumsum(np.arange(nShapes-1, 0, -1))
    cuts = cuts / (cuts[-1] / (cores - 1))
    start = [ np.nonzero(cuts==c)[0][0] for c in range((cores-1)) ]
    end = start[1:]
    end.append(nShapes)
    slices = zip(start,end)
    print slices

    rows = stage1Seq['rows']
    columns = stage1Seq['columns']
    poly2Rows = stage1Seq['poly2Rows']
    poly2Columns = stage1Seq['poly2Columns']
    stage2 = {}
    for c in range(cores - 1):
        pids = range(slices[c][0], slices[c][1])
        print pids
        stage2[c] = pool.apply_async(check_contiguity, args=(shapes, bBoxes,
            poly2Rows, poly2Columns, rows, columns, pids))
    pool.close()
    pool.join()
    results2 = {}
    for c in stage2:
        results2[c] = stage2[c].get()

    neighbors = dict([(i,set()) for i in range(nShapes) ])

    for c in results2:
        for key in results2[c]:
            neighbors[key] = neighbors[key].union(results2[c][key])

    wp2 = ps.W(neighbors)

    if wp2.neighbors == wseq.neighbors:
        print 'Parallel == Sequential'
    else:
        print 'Parallel != Sequential'
    
    sfName = "nat.shp"

    t1 = time.time()
    sf = ps.examples.get_path(sfName)
    shpFileObject = ps.open(sf)
    shapeBox = shpFileObject.bbox
    nShapes = len(shpFileObject)
    shapes = []
    bBoxes = []
    for shape in shpFileObject:
        shapes.append(shape)
        bBoxes.append(shape.bounding_box[:])

    # figure out grid that will be used by all processes
    DELTA = 0.000001
    bucketMin = nShapes / 80 + 2
    lengthX = ((shapeBox[2] + DELTA) - shapeBox[0]) / bucketMin
    lengthY = ((shapeBox[3] + DELTA) - shapeBox[1]) / bucketMin
    minBox = shapeBox[:2] * 2 # [minx, miny, minx, miny]
    binWidth = [ lengthX, lengthY] * 2 # [lenX, lenY, lenX, lenY]

    # sequential
    stage1Seq = binShapes(shapes, bBoxes, minBox, binWidth, bucketMin)
    stage2Seq = check_contiguity(shapes, bBoxes, stage1Seq['poly2Rows'],
            stage1Seq['poly2Columns'], stage1Seq['rows'], stage1Seq['columns'])


    wseq = ps.W(stage2Seq)
    t2 = time.time()
    print 'Sequential time: ', t2 - t1

    t1 = time.time()
    sf = ps.examples.get_path(sfName)
    shpFileObject = ps.open(sf)
    shapeBox = shpFileObject.bbox
    nShapes = len(shpFileObject)
    shapes = []
    bBoxes = []
    for shape in shpFileObject:
        shapes.append(shape)
        bBoxes.append(shape.bounding_box[:])

    # figure out grid that will be used by all processes
    DELTA = 0.000001
    bucketMin = nShapes / 80 + 2
    lengthX = ((shapeBox[2] + DELTA) - shapeBox[0]) / bucketMin
    lengthY = ((shapeBox[3] + DELTA) - shapeBox[1]) / bucketMin
    minBox = shapeBox[:2] * 2 # [minx, miny, minx, miny]
    binWidth = [ lengthX, lengthY] * 2 # [lenX, lenY, lenX, lenY]


    stage1Seq = binShapes(shapes, bBoxes, minBox, binWidth, bucketMin)

    cores = mp.cpu_count() - 1
    pool = mp.Pool(cores)
    cuts = np.cumsum(np.arange(nShapes-1, 0, -1))
    cuts = cuts / (cuts[-1] / (cores - 1))
    start = [ np.nonzero(cuts==c)[0][0] for c in range((cores-1)) ]
    end = start[1:]
    end.append(nShapes)
    slices = zip(start,end)

    rows = stage1Seq['rows']
    columns = stage1Seq['columns']
    poly2Rows = stage1Seq['poly2Rows']
    poly2Columns = stage1Seq['poly2Columns']
    stage2 = {}
    for c in range(cores - 1):
        pids = range(slices[c][0], slices[c][1])
        stage2[c] = pool.apply_async(check_contiguity, args=(shapes, bBoxes,
            poly2Rows, poly2Columns, rows, columns, pids))
    pool.close()
    pool.join()
    results2 = {}
    for c in stage2:
        results2[c] = stage2[c].get()

    neighbors = dict([(i,set()) for i in range(nShapes) ])

    for c in results2:
        for key in results2[c]:
            neighbors[key] = neighbors[key].union(results2[c][key])

    wp2 = ps.W(neighbors)
    t2 = time.time()
    print 'Parallel (seq stage 1, par stage 2)', t2 - t1

    if wp2.neighbors == wseq.neighbors:
        print 'Same weights'
    else:
        print 'different weights'

    t1 = time.time()
    wpysal = ps.queen_from_shapefile(sf)
    t2 = time.time()
    print 'Sequential pysal: ', t2 - t1
    nBad = 0
    for key in wpysal.neighbors:
        if set(wpysal.neighbors[key]) != wp2.neighbors[key]:
                nBad += 1
    if nBad == 0:
        print 'Same weights as pysal: ', True
    else:
        print 'Same weights as pysal: ', False

    """


