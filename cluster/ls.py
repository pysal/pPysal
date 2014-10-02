import collections
import multiprocessing as mp
from operator import gt, lt
import random
from random import randint, uniform

import pysal as ps
from pysal.region.components import check_contiguity

import numpy as np
from numpy.random import RandomState

class LocalSearch(mp.Process):
    """
    Attributes
    ----------

    failures                int     the current number of failures for this iteration
    intensificationsize     int     the size of the soln space to propagate the best soln
    tabulist                deque   of tuples in the form (unit, move to region)
    
    """
    def __init__(self ,attribute, w, nregions, lock = None, pid=None, floor=3,
            maxfailures=50, maxiterations=15, intensification=0.5):
        mp.Process.__init__(self)
        self.index = pid
        self.lock = lock
        self.w = w
        self.z = attribute
        self.floor = floor
        self.maxiterations = maxiterations
        self.wss = 0

        #Shared memory setup
        self.solnspace = np.frombuffer(shared_solnspace.get_obj(), dtype=np.float32)
        self.solnspace.shape = (-1, self.w.n + 1)
        self.solnspacesize = self.solnspace[:,0].size
        self.nregions = nregions
        #Work on a copy of the shared memory space.
        self.solncolumn = np.empty(self.w.n)
        self.unitchooser = np.arange(len(self.solncolumn))

        #Setup for intensification and diversification
        self.intensificationsize = int(self.solnspacesize * intensification)

        #Tabu parameters
        self.failures = 0
        self.maxfailures = maxfailures + int(maxfailures * uniform(-1.1, 1.2))
        self.maxiterations = 15
        self.tabulength = self.computetabulength()
        self.tabulist = collections.deque(maxlen=self.tabulength)

        #Seed the python random number generator
        self.randstate = RandomState(self.index)
        random.seed(self.index)
        
    def __repr__(self):
        return """
        The current state of index {} is:
            working on soln {} / {}
            current region membership: 
{}
            current tabu list length: {}
            current maximum number of failures: {}
            current obj. func value: {}
        """.format(self.index, self.index, self.solnspacesize,
                self.solncolumn.reshape(8,8), self.tabulength, self.maxfailures, self.wss)

    def localsearch(self):
        while self.failures < self.maxfailures: 
            sln = self.solncolumn
            ids = np.arange(len(sln))
            #Select a random atomic unit
            selectedunit = self.randstate.choice(self.unitchooser)
            
            #Check the floor constraint
            region_membership = sln[selectedunit]
            
            mask = sln[sln == region_membership]
            unit_region_count = len(mask)
            #units_in_region = np.where(sln == region_membership)[0].tolist()
            units_in_region = ids[sln == region_membership].tolist()
            if unit_region_count == self.floor:
                #A move from this region would break the floor constraint
                self.failures += 1
                continue
            
            #Get the neighbors and remove the current region
            neighborunits = self.w.neighbors[selectedunit]  #List of neighbors from the W Obj
            neighborregions = set(sln[neighborunits]) #List of neighbor regions
            #print "Regions: {}, Unit: {}".format(neighborregions, sln[selectedunit])
            neighborregions.remove(sln[selectedunit]) #Remove the region the unit is currently in
            
            if not neighborregions:
                #The unit is internal to the region
                self.failures += 1
                continue

            possibleswaps = {}    
            for n in neighborregions:
                sln[selectedunit] = n
                wss = self.objective_func(regions=sln)
                possibleswaps[wss] = n
                sln[selectedunit] = region_membership
                
            #If no swaps exist, get another atomic unit
            if not possibleswaps:
                self.failures += 1
                continue
            #Iterate through the swaps
            for k in sorted(possibleswaps):
                #Check contiguity
                contigious = check_contiguity(self.w, units_in_region, selectedunit)
                if not contigious:
                    self.failures += 1
                    continue
                #Check the tabu list
                if (selectedunit, possibleswaps[k]) in self.tabulist:
                    #TODO: Add an aspiration function here that allows a tabu move to procedd
                    # if it is 'good enough'
                    self.failures += 1
                    continue

                #The move is valid
                self.tabulist.appendleft((selectedunit, possibleswaps[k]))
                sln[selectedunit] = possibleswaps[k]  #where the value is the new region id
                self.wss = k  #where k is the current wss
                self.failures = 0
        
        with self.lock:
            currentwss = self.solnspace[:,0][self.index]
            maxwss = np.max(self.solnspace[:,0])
                    
            if self.wss < currentwss:
                self.solnspace[self.index][0] = self.wss
                self.solnspace[self.index][1:] = sln

                #Check to see if this is also a global best
                if self.wss < maxwss:
                    pass
                    #If so, propagate to some % of the solution space
        return

    def run(self):
        #Populate the initial objective function value
        with self.lock:
            self.wss = self.solncolumn[0]
        while self.maxiterations > 0: 
            with self.lock:
                #Populate the local working space
                self.solncolumn[:] = self.solnspace[self.index][1:]
                #Compute the current objective function value
                self.wss = self.solnspace[self.index,0]
            #Diversification occurs here, before local search
            self.failures = 0  #Reset the failure counter before each iteration
            self.localsearch()
            
            #This is a constant contiguity check that can be removed once validated.
            cont = True
            for i in range(1, int(self.nregions) + 1):
                region = np.where(self.solncolumn == i)[0].tolist()

                if test_region(self.w, region) == False:
                    cont = False
            if cont == False:
                print "ERROR: ", self.__repr__()
            
            #Increment the index counter to step around the solution space
            self.index += 1
            if self.index > self.solnspacesize:
                self.index = 0
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

    def objective_func(self, regions=None):
        """
        Computes the global objective function value
        """
        wss = 0
        for r in range(1, int(self.nregions) + 1):
            ids = np.where(regions == r)[0]
            m = self.z[ids]
            var = m.var()
            wss += np.sum(var * len(ids))
        return wss 


def initshared_localsoln(_solnspace):
    global shared_solnspace
    shared_solnspace = _solnspace

def test_region(w,neighbors):
    d={}
    g=Graph()
    for i in neighbors:
        d[i]=[j for j in w.neighbors[i] if (j in neighbors)]
    for i in d:
        for j in d[i]:
            g.add_edge(i,j,1.0)
    cc=g.connected_components(op=gt)
    if len(cc)==1:
        return True
    else:
        return False

class Graph(object):
    def __init__(self):
        self.nodes=set()
        self.edges={}
        self.cluster_lookup={}
        self.no_link={}

    def add_edge(self,n1,n2,w):
        self.nodes.add(n1)
        self.nodes.add(n2)
        self.edges.setdefault(n1,{}).update({n2:w})
        self.edges.setdefault(n2,{}).update({n1:w})

    def connected_components(self,threshold=0.9, op=lt):
        nodes = set(self.nodes)
        components,visited =[], set()
        while len(nodes) > 0:
            connected, visited = self.dfs(nodes.pop(), visited, threshold, op)
            connected = set(connected)
            for node in connected:
                if node in nodes:
                    nodes.remove(node)
            subgraph=Graph()
            subgraph.nodes = connected
            subgraph.no_link = self.no_link
            for s in subgraph.nodes:
                for k,v in self.edges.get(s,{}).iteritems():
                    if k in subgraph.nodes:
                        subgraph.edges.setdefault(s,{}).update({k:v})
                if s in self.cluster_lookup:
                    subgraph.cluster_lookup[s] = self.cluster_lookup[s]
            components.append(subgraph)
        return components
    
    def dfs(self, v, visited, threshold, op=lt, first=None):
        aux=[v]
        visited.add(v)
        if first is None:
            first = v
        for i in (n for n, w in self.edges.get(v,{}).iteritems() \
                  if op(w, threshold) and n not in visited):
            x,y=self.dfs(i,visited,threshold,op,first)
            aux.extend(x)
            visited=visited.union(y)
        return aux, visited
