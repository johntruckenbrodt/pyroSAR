#
# Modul zur Verwaltung von Graphen (Adjazenzlistendarstellung)
#
# 10. Januar 2013
# Robert Zeranski
#
# edited March 2016
# John Truckenbrodt


class Graph(object):
    def __init__(self, n):
        self.nodeList = [{} for _ in range(0, n + 1)]

    @property
    def numNodes(self):
        return len(self.nodeList)-1

    @property
    def v(self):
        return range(0, self.numNodes)

    @property
    def e(self):
        return [(i, j) for i in self.v for j in self.v if self.isEdge(i, j)]

    def isNode(self, u):
        return 0 <= u < self.numNodes

    def addNode(self, x):
        for i in range(0, x):
            self.nodeList.append({})

    def isEdge(self, u, v):
        return False if not self.isNode(u) and not self.isNode(v) else v in self.nodeList[u].keys()

    def addEdge(self, u, v):
        if not self.isNode(u) or not self.isNode(v):
            raise IndexError("node {0} or {1} does not exist".format(u, v))
        else:
            if not self.isEdge(u, v):
                self.nodeList[u][v] = float("inf")

    def setWeight(self, u, v, weight):
        if self.isEdge(u, v):
            self.nodeList[u][v] = float(weight)
        else:
            raise IndexError("edge ({0}, {1}) does not exist".format(u, v))

    def w(self, u, v):
        return self.nodeList[u][v] if self.isEdge(u, v) else float("inf")

    def n(self, u):
        return self.nodeList[u].keys() if self.isNode(u) else []

    def output(self):
        print self.v
        print self.e
        print [self.w(i, j) for i in self.v for j in self.v if self.isEdge(i, j) and self.w(i, j) < float("inf")]
