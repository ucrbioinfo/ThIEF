import numpy as np
from copy import deepcopy


# note the assumption that location of the nucleosome is the first element in
# the corresponding vector

class layerOfNodes:
    # list of nodes in the layer
    nodes = []

    def __init__(self, listNodes=None):
        '''(listNodes=None)
        Initialize node layer with nodes from the listNodes
        '''
        self.nodes = []
        if listNodes is None:
            # self.nodes = []
            pass
        else:
            for node in listNodes:
                self.nodes.append(node)

    def __del__(self):
        print 'deleting layer'
        try:
            self.nodes = []
        except:
            print 'deleting layer failed'
            print len(self.nodes)

    def GetNodes(self, id=None):
        if id is None:
            return self.nodes
        else:
            try:
                return self.nodes[id]
            except:
                print 'smth went wrong while getNodes(id)'

    def Location(self, listIDs=None):
        '''(listIDs)
        Returns a list of locations for the nodes in the listIDs.
        If listIDs is not provided all of the locations are returned
        '''
        if listIDs is None:
            output = []
            try:
                for x in self.nodes:
                    # assuming location is the first element in the node vector
                    output.append(x[0])
            except:
                print "smth is wrong with the list of nodes in the layer"
            return output
        else:
            output = []
            try:
                for x in listIDs:
                    # assuming location is the first element in the node vector
                    output.append(self.nodes[x][0])
            except:
                print "smth went wrong with layerOfNodes location"
            return output

    def Vicinity(self, loc, size):
        '''(loc, size)
        returns list of IDs for nodes in the layer that are trapped in vicinity
        (loc, size)
        '''
        res = []
        idx = -1
        for node in self.nodes:
            idx += 1
            if abs(node[0] - loc) < size:
                res.append(idx)
        return res

    def Correct(self):
        '''
        Returns whether layerOfNodes is correct
        '''
        try:
            ndim = len(self.nodes)
        except:
            return False
        if ndim > 0:

            for i in xrange(ndim):
                # assuming index of node location is first element in the
                # vector
                if self.nodes[i][0] < 0:
                    return False
                else:
                    return True
        else:
            return False


class hyperGraph:
    # number of layers in the graph
    nLayers = 0
    # list of nodes
    nodes = []
    # list of layers
    layers = []
    # list of hyperedges (simple list of node IDs from list nodes)
    edges = []
    # list of corresponding hyperedge costs
    costs = []
    # list of corresponding hyperedge locations
    edgeLocs = []
    # vicinity Radius
    okil = 100
    # penalty for missing node
    missPenalty = 50

    def Clean(self):
        '''
        cleans the graph. This is patch code that is not suppose to be here
        Don't do code like this
        '''
        self.nLayers = 0
        self.nodes = []
        self.layers = []
        self.edges = []
        self.costs = []
        self.edgeLocs = []

    def __init__(self, layer=None):
        '''(listNodes=None)
        Initialize hyperGraph with nodes from layer

        '''
        self.nLayers = 0
        self.nodes = []
        self.layers = []
        self.edges = []
        self.costs = []
        self.edgeLocs = []
        if layer is None:
            pass
        else:
            self.layers.append(layer)
            ind = -1
            for node in layer.GetNodes():
                ind += 1
                self.edges.append([ind])
                # assuming location is the first element in the vector
                if hasattr(node,"__len__"):
                    tmp = node[0]
                else:
                    tmp = node
                self.edgeLocs.append(tmp)
                self.costs.append(0)
            self.nLayers = self.nLayers + 1
            # self.nodes = np.array(nodes) #break functionality @addLayer()

    def __del__(self):
        print 'deleting graph'
        try:
            self.nodes = []
        except:
            print 'deleting graph.nodes failed'
        try:
            for layer in self.layers:
                del layer
            self.layers = []
        except:
            print 'deleting graph.layers failed'
        try:
            self.edges = []
        except:
            print 'deleting graph.edges failed'
        try:
            self.costs = []
        except:
            print 'deleting graph.costs failed'
        try:
            self.edgeLocs = []
        except:
            print 'deleting graph.edgeLocs failed'
        self.Clean()

    def EdgeNode(self, edgeID, edgePos):
        '''(int ID, int POS)
        returns node attributes for the ID hyperedge's node at position POS
        '''
        try:
            nodeID = self.edges[edgeID][edgePos]
        except IndexError:
            print "node index out of bounds for current edge"
            nodeID = -1
        if nodeID >= 0:
            return self.layers.getNodes(nodeID)
        else:
            return -1

    def AddLayer(self, layerOfNodes):
        '''(layerOfNodes)
        adds a new layer of nodes to a hypergraph
        '''
        if layerOfNodes.Correct():
            # add new layer of nodes to the list
            self.layers.append(layerOfNodes)
            # prolong each hyperedge with dummy node
            edgesToDull = deepcopy(self.edges)
            for x in edgesToDull:
                x.append(-1)

            # create new edges for every node in the new layer
            dullToLayer = []
            count = 0
            for n in layerOfNodes.GetNodes():
                tmpList = []
                for i in xrange(self.nLayers):
                    tmpList += list([-1])
                dullToLayer.append(tmpList + [count])
                # assuming location is the first elem in node vector
                self.edgeLocs.append(n[0])
                self.costs.append(0)
                count += 1

            # create new edges for extension of existing edges
            newEdges = []
            ind = -1
            for t in self.edges:
                ind += 1
                candidateNodes = layerOfNodes.Vicinity(
                    self.edgeLocs[ind], self.okil)
                if len(candidateNodes) > 0:
                    for nodeID in candidateNodes:
                        newEdge = t + [nodeID]
                        newEdges.append(newEdge)
                        newEdgeLocation = self.EdgeAveLoc(newEdge)
                        self.edgeLocs.append(newEdgeLocation)
                        self.costs.append(0)
            self.edges = edgesToDull + dullToLayer + newEdges
            self.nLayers += 1
        else:
            print 'smth is wrong with the layer to add'

    def EdgeAveLoc(self, edge):
        '''(edge)
        returns an average hyper edge location
        '''
        cumLoc = 0
        cumCount = 0
        curLayer = -1
        for edgeID in edge:
            curLayer += 1
            if edgeID >= 0:
                try:
                    loc = self.layers[curLayer].Location([edgeID])
                    increment = 1
                except:
                    print 'smth wrong with extracting edge locations', curLayer, edgeID
                    loc = 0
                    increment = 0
                cumLoc += loc[0]
                cumCount += increment
        if cumCount > 0:
            return cumLoc / cumCount
        else:
            return -1000

    def EdgeLoc(self, edge):
        '''(edge)
        returns a vector of hyper edge locations (NB!!! order doesn't matter)
        '''
        cumLoc = []
        cumCount = 0
        curLayer = -1
        for edgeID in edge:
            curLayer += 1
            if edgeID >= 0:
                try:
                    loc = self.layers[curLayer].Location([edgeID])
                    increment = 1
                except:
                    print 'smth wrong with extracting edge locations', curLayer, edgeID
                    loc = 0
                    increment = 0
                cumLoc.append(loc[0])
                cumCount += increment
        if cumCount > 0:
            return cumLoc
        else:
            return None

    def EdgeCostComputation(self):
        '''
        Updates the cost for every hyper edge in the graph
        '''
        try:
            count = -1
            for edge in self.edges:
                count += 1
                loc = self.edgeLocs[count]
                locations = self.EdgeLoc(edge)
                cumCost = 0
                for x in locations:
                    cumCost += abs(x - loc)
                self.costs[count] = cumCost + (
                    self.nLayers - len(locations)) * self.missPenalty
        except:
            print 'smth went wrong with updating edge costs'

    def NodeToEdgeList(self, layerID, nodeID):
        '''(layerID, nodeID)
        given the node report all edgeIDs that contain that node
        '''
        edgeList = []
        countID = -1
        for edge in self.edges:
            countID += 1
            if edge[layerID] == nodeID:
                edgeList.append(countID)
        return edgeList

    def GetTrackStat(self, listOfEdgeIDs, numOfAttributes=1):
        '''(listOfEdgeIDs)
        given the list of hyperedge IDs return a full table of nucleosome stats
        '''
        table = []
        emptyPlace = -1 * np.ones(numOfAttributes)
        for edgeID in listOfEdgeIDs:
            currentEdge = self.edges[edgeID]
            curLayer = -1
            trackAttr = np.array([])
            for nodeID in currentEdge:
                curLayer += 1
                if nodeID >= 0:
                    nodeAttr = self.layers[curLayer].nodes[nodeID]
                else:
                    nodeAttr = emptyPlace 
                trackAttr = np.concatenate((trackAttr, nodeAttr))
            table.append(trackAttr)
        return np.array(table)


def CPLEXprint(graph, fileName):
    numVars = len(graph.edges)
    numLayers = graph.nLayers
    tmpSum = 0
    for layer in graph.layers:
        tmpSum += len(layer.Location())
    numNodes = tmpSum
    with open(fileName, 'w') as fout:
        # code for printing objective function goes here
        print >> fout, 'Minimize'
        objectString = ''
        for i in xrange(numVars):
            objectString += ' +' + str(graph.costs[i]) + ' x' + str(i)
        print >> fout, 'object:', objectString
        # code for printing constraints goes here
        print >> fout, 'Subject To'
        countLayer = -1
        for layer in graph.layers:
            countLayer += 1
            countNode = -1
            for node in layer.GetNodes():
                countNode += 1
                boundaryString = 'l' + str(
                    countLayer) + 'n' + str(countNode) + ':'
                coveringEdges = graph.NodeToEdgeList(countLayer, countNode)
                if len(coveringEdges) > 0:
                    for ind in coveringEdges:
                        boundaryString += ' + x' + str(ind)
                    boundaryString += ' = 1'
                    print >> fout, boundaryString
                else:
                    print >> 'l' + str(countLayer) + 'n' + str(
                        countNode) + ' has no covering edges'
        # code for printint variable bountadires goes here
        print >> fout, 'Bounds'
        for i in xrange(numVars):
            print >> fout, '0 <= x' + str(i) + ' <= 1'


def main():
    # test = hyperGraph()
    # test.init([11, 22, 33])
    return 0


if __name__ == '__main__':
    # main()
    pass
