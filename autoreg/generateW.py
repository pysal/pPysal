import sys
import os
import cPickle
import pysal as ps

#Hack to get the KDTree used in GWK to pickle.
# KDTree uses nested classes that the pickler can not find.
from scipy.spatial import kdtree
kdtree.node = kdtree.KDTree.node
kdtree.leafnode = kdtree.KDTree.leafnode
kdtree.innernode = kdtree.KDTree.innernode

def main(datafile, adjacency):
    import time
    t1 = time.time()
    if adjacency.lower() == 'queen':
        w = ps.weights.user.queen_from_shapefile(datafile)
    elif adjacency.lower() == 'rook':
        w = ps.weights.user.rook_from_shapefile(datafile)
    w.transform = 'r'
    w.adjacency = adjacency
    t2 = time.time()
    print "Time W: ", t2 - t1

    gwk = ps.weights.user.kernelW_from_shapefile(datafile, 5, diagonal=True)
    gwk.adjacency = adjacency
    t3 = time.time()
    print "Time GWK: ", t3 - t2
    #Write W to disk
    wout = datafile.split('.')[0] + '_W.pkl'
    gwkout = datafile.split('.')[0] + '_GWK.pkl'

    try:
        os.remove(wout)
    except OSError:
        pass

    try:
        os.remove(gwkout)
    except OSError:
        pass

    t1 = time.time()
    with open(wout, 'wb') as f:
        cPickle.dump(w, f, cPickle.HIGHEST_PROTOCOL)
    with open(gwkout, 'wb') as f:
        cPickle.dump(gwk, f, cPickle.HIGHEST_PROTOCOL)
    t2 = time.time()
    print "Dumping: ",t2 - t1

    return wout, gwkout
    """
    #Uncomment to test
    with open(wout, 'rb') as f:
        w2 = cPickle.load(f)
    with open(gwkout, 'rb') as f:
        gwk2 = cPickle.load(f)

    assert w.cardinalities == w2.cardinalities
    """

if __name__ == '__main__':
    datafile = sys.argv[1]
    adjacency = sys.argv[2]

    main(datafile, adjacency)


