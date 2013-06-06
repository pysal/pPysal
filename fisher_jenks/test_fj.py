from fj_refactored import fisher_jenks
import numpy
import time

cores = [1,2,4,16,32]
classes = [5,6,7]
data_sizes = [1000,2000,4000,8000,16000,24000,32000,40000,48000,56000,64000,72000,80000]

for c in cores:
    for d in data_sizes:
        for k in classes:
            data = numpy.random.ranf(size=d)
            try:
                t1 = time.time()
                #wrapped in try since we will blow out RAM at some point
                classification = fisher_jenks(data, k)
                t2 = time.time()
                print "Processed {} data points in {} classes using {} cores.  Total time: {}".format(d, k, c, t2-t1)
            except:
                print "FAILURE: {} Data Points.".format(d)