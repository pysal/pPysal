import collections
import multiprocessing as mp
import random
from random import randint, uniform

import pysal as ps
import numpy as np
from numpy.random import RandomState

class LocalSearch(mp.Process):

    def __init__(self ,attribute, w, lock = None, pid=None, floor=3,
            maxfailures=50, maxiterations=15):
        mp.Process.__init__(self)
        self.index = pid
        self.lock = lock
        self.w = w
        self.z = attribute
        self.floor = floor
       
        #Shared memory setup
        self.solnspace = np.frombuffer(shared_solnspace.get_obj(), dtype=np.int32)
        self.solnspace.shape = (-1, self.w.n + 1)
        self.solnspacesize = self.solnspace[:,0].size
        self.nregions = self.solnspace[0,0]
        
        #Tabu parameters
        self.maxfailures = maxfailures + int(maxfailures * uniform(-1.1, 1.2))
        self.maxiterations = 15
        self.tabulength = self.computetabulength()
        self.tabulist = collections.deque(maxlen=self.tabulength)

        #Compute the current objective function value
        self.getglobalobjective()

        #Seed the python random number generator
        random.seed(self.index)

    def __repr__(self):
        return """
        The current state of index {} is:
            working on soln {} / {}
            current region membership: {}
            current tabu list length: {}
            current maximum number of failures: {}
        """.format(self.index, self.index, self.solnspacesize,
                self.solncolumn, self.tabulength, self.maxfailures)

    def run(self):
        while self.maxiterations > 0:
            self.solncolumn = self.solnspace[self.index][1:]
            sln = self.solncolumn
            failures = 0
            while failures <= self.maxfailures:
                regionids = range(1,self.nregions + 1)
                
                #Seed the random number generator
                randomstate = RandomState(self.index)
                randomstate.shuffle(regionids)
                print regionids
                #Iterate through the regions, checking potential swaps
                
                nset = set()
                for r in regionids:
                   neighbors = np.where(sln == r)[0]
                   

                #Placeholder to escape the while loop
                failures += 1

            print self.__repr__()
            self.maxiterations -= 1

    def computetabulength(self):
        '''Talliard 1990'''
        smin = (self.nregions - 1) * 0.9
        smax = (self.nregions - 1) * 1.1
        tabu_length = 6 + (randint(0, int(smax - smin)))
        return tabu_length

    def intensifysoln(self):
        """
        Method to populate some percentage of the best soln to the global space
        """
        pass

    def diversifysoln(self):
        """
        Method to randomize a non-improving solution
        """
        pass

    def getglobalobjective(self):
        """
        Computes the global objective function value
        """


def initshared_localsoln(_solnspace):
    global shared_solnspace
    shared_solnspace = _solnspace
