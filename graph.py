#!/bin/python3

from exceptions import EdgeDoesNotExists, NoCompatibleConnectionFound
from itertools import combinations

import pprint 

class Weight:

    def __init__(self, w):
       self.vector = w
       
    def __add__(self, other):
        sum = self.vector + other.vector
        return Weight(sum)

    def __sub__(self, other):
        sum = (self.vector[0] - other.vector[0], self.vector[1] - other.vector[1])
        return Weight(sum)

    def __rmul__(self, other):
        return Weight([other*entry for entry in self.vector])

    def __repr__(self) -> str:
        return f"({self.vector[0]}, {self.vector[1]})"

def GKMCondition(w1: Weight, w2: Weight, w3: Weight) -> bool:
    """
    Consider w3 = nabla_w1 w2 = w2 + c*w1 
    (overloaded notation, w are weights and edges...)
    Return True if such a c exists, otherwise False

    AT THE MOMENT, IT WORKS ONLY FOR WEIGHTS IN Z^2
    """
    
    w = w3 - w2
    det = w.vector[0]*w1.vector[1] - w.vector[1]*w1.vector[0] 

    if not det == 0: return False

    try:
        if w.vector[0] % w1.vector[0] == 0:
            return True
        else:
            return False
    except ZeroDivisionError:
        if w.vector[1] % w1.vector[1] == 0:
            return True
        else:
            return False

    
class Edge:

    def __init__(self, v1, v2, w = None) -> None:
        self.v1 = v1
        self.v2 = v2
        self.weight = w
        self.name = (v1,v2)
        self.hash = f"{v1}{v2}"

    def setWeight(self, w : Weight):
        self.weight = w
        return self

    def __repr__(self) -> str:
        r = ""

        return f"({self.v1},{self.v2}){r}"

    def __eq__(self, other) -> bool:
        if (self.v1 == other.v1 and self.v2 == other.v2):
            return True
        else:
            return False

    def __hash__(self):
        return hash(self.hash)

def findEdge(edges: list, v1, v2)->Edge:

    for e  in edges:
        if (e.v1 == v1 and e.v2 == v2) or (e.v1 == v2 and e.v2 == v1):
            return e
    
    raise EdgeDoesNotExists("Edge does not exist")

class Connection:

    def __init__(self) -> None:
        self.con = {}

    def setConnection(self, e1: Edge, e2: Edge, e3: Edge):
        self.con[e1] = { e2 : e3 }
        return self

        
class Graph:

    """
    In the configuration file the first column represents a vertex und the other
    columns indicate the vertex to which the former one is connected.
    """

    def __init__(self, configuration, connection = None) -> None:
        self.vertices = []
        self.edges = []
        self.connection = {}

        with open(configuration , "r") as file:
            for line in file.readlines():
                items = line.strip().split(" ")
                for item in items:
                    self.vertices.append(item) if not item in self.vertices else self.vertices
                    if not item == items[0]:
                        self.edges.append(Edge(items[0],item))

        if connection:
            with open(connection, "r") as file:
                for line in file.readlines():
                    i1 = line.split(";")[0].split(",")

                    for edge in self.edges:
                        if edge.v1 == i1[0] and edge.v2 == i1[1]:
                            e1 = edge

                    if not e1 in self.connection:
                        self.connection[e1] = {}

                    i2 = line.split(";")[1].split(",")
                    i3 = line.split(";")[2].split(",")

                    for edge in self.edges:
                        if edge.v1 == i2[0] and edge.v2 == i2[1]:
                            e2 = edge
                        if edge.v1 == i3[0] and edge.v2 == i3[1].strip():
                            e3 = edge

                    self.connection[e1][e2] = e3
                    self.connection[e1][e3] = e2
                    

    def findEdge(self, v1, v2) -> Edge:
        e = Edge(str(v1),str(v2))
        for edge in self.edges:
            if edge == e:
                return edge
        else:
            raise EdgeDoesNotExists(f"Edge: ({e.v1},{e.v2}) does not exists!")
         
    def returnEdges(self, name):
        with open(f"edges_{name}.gkm", "w") as file:
            for edge in self.edges:
                file.write(f"{edge.v1},{edge.v2};\n")

    def loadWeights(self, name):
        with open(f"edges_{name}.gkm", "r") as file:
            for line in file.readlines():
                items = line.strip().split(";")
                e1 = items[0].split(",")[0]
                e2 = items[0].split(",")[1]
                w1 = items[1].split(",")[0]
                w2 = items[1].split(",")[1]

                e = self.findEdge(e1, e2)
                e.setWeight(Weight((int(w1),int(w2))))

        return self
                
    def connectionPath(self, e1: Edge, e2: Edge):
        """
        We start with nabla_e1 (e2)
        """
        path = []
        start = e1
        end = e2
        next = self.connection[e1][e2]
        prev = e1
        path.append(e1)

        while True:

            if next == end and self.connection[next][prev] == start:
                path.append(end)
                break

            path.append(next)
            p = next
            # nabla_next (prev)
            next = self.connection[next][prev]
            prev = p

        return path

    def emanatingEdges(self, vertex):
        edges = []

        for edge in self.edges:
            if edge.v1 == vertex or edge.v2 == vertex:
                edges.append(edge)

        return edges

    def getEpsilon(self, e1: Edge, e2: Edge):
        e3 = self.connection[e1][e2]

        w1 = e1.weight
        w2 = e2.weight
        w3 = e3.weight

        if GKMCondition(w1,w2,w3):
            return 1
        if GKMCondition(w1,(-1)*w2,w3):
            return -1 

    def computeEta(self, vertex, pathEdge: Edge):
        edges = self.emanatingEdges(vertex)
        eta = -1

        for edge in edges:
            if not edge == pathEdge:
                eps = self.getEpsilon(pathEdge, edge)
                eta = eta * eps

        return eta

    def computeOrientationPath(self, path: list):
        r = 1
        for e in path:
            eta = self.computeEta(e.v1, e)
            r = r * eta
        return r

    def computeAllConnectionPaths(self):
        paths = []
        set_paths = []

        for vertex in self.vertices:
            edges = self.emanatingEdges(vertex)
            pairs = list(combinations(edges, 2))
            for pair in pairs:
                new = self.connectionPath(pair[0],pair[1]) 
                if not set(new) in set_paths:
                    paths.append(new)
                    set_paths.append(set(new))
        
        return paths

    def createConnection(self):
        for edge in self.edges:
            self.connection[edge] = {}
            initial_edges = self.emanatingEdges(edge.v1)
            terminal_edges = self.emanatingEdges(edge.v2)

            initial_edges.remove(edge)
            terminal_edges.remove(edge)

            for initial in initial_edges:
                for terminal in terminal_edges:
                    if GKMCondition(edge.weight, initial.weight,
                                    terminal.weight):
                        #print(edge, initial, terminal)
                        self.connection[edge][initial] = terminal
                        self.connection[edge][terminal] = initial
                        terminal_edges.remove(terminal)
                        break
                    elif GKMCondition(edge.weight, (-1)*initial.weight,
                                    terminal.weight):
                        #print(edge, initial, terminal)
                        self.connection[edge][initial] = terminal
                        self.connection[edge][terminal] = initial
                        terminal_edges.remove(terminal)
                        break
                    else:
                        #print(f"nabla_{edge}{initial}={terminal} violates GKM condition")
                        pass

            if self.connection[edge] == {}:
                raise NoCompatibleConnectionFound(f"No compatible connection found for edge{edge}")

def vgraph():
    vg = Graph("vater_graph.gkm")

    vg.loadWeights("vater_graph")
    #p = [ (1,3), (3,4), (4,6), (6,7), (7,8), (8,11), (11,12), (1,12) ]
    p = [(4,5), (5,6), (6,7), (7,9), (8,9), (7,8), (6,7), (4,6)]
    outer = [ vg.findEdge(e[0],e[1]) for e in p]
    vg.createConnection()
    print(p, vg.computeOrientationPath(outer))
    #print(vg.computeOrientationPath(outer))
    #for path in vater_graph.computeAllConnectionPaths():
    #    print(path, vater_graph.computeOrientationPath(path))

def cp3():
    cp3 = Graph("cp3.gkm")
    cp3.loadWeights("cp3")
    cp3.createConnection()
    paths = cp3.computeAllConnectionPaths()
    for path in paths:
        print(path, cp3.computeOrientationPath(path))
    #connection_paths = cp3.computeAllConnectionPaths()
    #for path in connection_paths:
    #    print(cp3.computeOrientationPath(path))

def tolman():
    tolman = Graph("tolman.gkm")
    tolman.loadWeights("tolman")
    tolman.createConnection()
    pprint.pprint(tolman.connection)
    for path in tolman.computeAllConnectionPaths():
        print(path, tolman.computeOrientationPath(path))



def main():
    #cp3()
    #tolman()
    vgraph()
    

if __name__ == "__main__":
  main()

        
