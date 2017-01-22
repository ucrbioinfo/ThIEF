from random import randrange, random, gauss
import numpy as np
import sys
import TrNM as hg
import os
import LPsolMerge as lp


def ShakeLayer(layer0=None, shakeAmplitude=25, probToDisappear=0.1):
    """
    Function takes input vector/matrix layer0 and randomly "shakes" the 0 coordinate with shakeAmplitude
    returns a matrix (location, isPresentFlag), where first coordinate is location, second - is it present
    """
    if layer0 is None:
        print "no base layer"
        return -1
    else:
        nLen = len(layer0)
        layer1 = []
        for i in range(nLen):
            x = layer0[i, 0] + gauss(0,shakeAmplitude)
            y = (random() > probToDisappear) * 1
            layer1.append([int(x), int(y)])
        return np.array(layer1)


def SaveToFile(layer, fileName):
    from os.path import isfile
    # if isfile(fileName):
    #     print "file " + fileName + " already exists"
    #     print "output will go to " + fileName + "_1"
    #     fileName = fileName + '_1'
    with open(fileName, 'w') as fout:
        for line in layer:
            if line[1] == 1:
                print >> fout, line[0]


def SaveTracksToFile(listOfLayers, fileName):
    from os.path import isfile
    # if isfile(fileName):
    #     print "file " + fileName + " already exists"
    #     print "output will go to " + fileName + "_1"
    #     fileName = fileName + '_1'
    numLayer = len(listOfLayers)
    numLines = len(listOfLayers[0])
    with open(fileName, 'w') as fout:
        for i in range(numLines):
            strToWrite = ''
            for j in range(numLayer):
                strToWrite += str(listOfLayers[j][i][0])+'\t'+ str(listOfLayers[j][i][1]) + '\t'
            print >> fout, strToWrite


def main(fileName='default.out', numPoints=1000):
    """function to generate a set of points on a plane that are at least "linker" apart
    returns a matrix (location, isPresentFlag), where first coordinate is location, second - boolen True if the point is present

    """
    probToDisappear = 0.1  # probablity of a node to disappear from the plane
    linker = 146  # size of the linker region between points
    if len(sys.argv) > 1:
        fileName = sys.argv[1]
    if len(sys.argv) > 2:
        linker = int(sys.argv[2])
    else:
        linker = 146 
    if len(sys.argv) > 3:
        shakeAmplitude = int(sys.argv[3])
    else:
        shakeAmplitude = 25
    if len(sys.argv) > 4:
        missPenalty = int(sys.argv[4])
    else:
        missPenalty = 30
    node0 = []
    x = 500  # starting point
    for i in range(numPoints):
        x += randrange(linker, linker * 2)
        y = (random() > probToDisappear) * 1
        node0.append([x, y])
    node0 = np.array(node0)
    node1 = ShakeLayer(node0, shakeAmplitude, probToDisappear)
    node2 = ShakeLayer(node1, shakeAmplitude, probToDisappear)
    SaveTracksToFile([node0, node1, node2], 'OUTPUT/'+fileName + '_3tracks')
    
## Uncomment this part in case you want tacks consisting of 4 layers 
    # node3 = ShakeLayer(node2, shakeAmplitude, probToDisappear)
    # SaveTracksToFile([node0, node1, node2, node3], 'OUTPUT/'+fileName + '_4tracks')


## Uncomment this part in case you want 5 layers of tracks instead of 4
#     if ('node3' not in locals()):
#         node3 = ShakeLayer(node2, shakeAmplitude, probToDisappear)
#     node4 = ShakeLayer(node3, shakeAmplitude, probToDisappear)
#     SaveTracksToFile([node0, node1, node2, node3, node4], 'OUTPUT/'+fileName + '_5tracks')


# ##this part for building solution
# ##NB!!! It requires gpsol installed (no check for this) 
   
#     layer0 = hg.layerOfNodes(node0[node0[:, 1] == 1, :])
#     layer1 = hg.layerOfNodes(node1[node1[:, 1] == 1, :])
#     layer2 = hg.layerOfNodes(node2[node2[:, 1] == 1, :])
#     flag3 = False
#     if ('node3' in locals()):
#        layer3 = hg.layerOfNodes(node3[node3[:, 1] == 1, :])
#        flag3 = True
#     flag4 = False
#     if ('node4' in locals()):
#         layer4 = hg.layerOfNodes(node4[node4[:, 1] == 1, :])
#         flag4 = True
#     print 'layers processed'
#     graph = hg.hyperGraph(layer0)
#     graph.missPenalty = missPenalty
#     graph.okil = missPenalty*2
#     print "layer0 done"
#     graph.AddLayer(layer1)
#     print "layer1 done"
#     graph.AddLayer(layer2)
#     print "layer2 done"
#     if flag3:
#         graph.AddLayer(layer3)
#         print "layer3 done"
#     if flag4:
#         graph.AddLayer(layer4)
#         print "layer4 done"
#     graph.EdgeCostComputation()
#     print 'done building graph'
#     hg.CPLEXprint(graph, fileName+'_tmp.lp')
#     # run linear solver gpsol and parse it's output
#     print 'start linear solver'
#     os.system('./runLS.sh '+fileName+'_tmp.lp')
#     print 'linear solution done'
#     lpSol = lp.ReadColumn(fileName+'_tmp.csv')
#     print 'solution read'
#     os.system('rm '+fileName+'_tmp.lp')
#     os.system('rm '+fileName+'_tmp.csv')
#     os.system('rm '+fileName+'_tmp.sol')
#     table = graph.GetTrackStat(lpSol)
#     # print table
#     np.savetxt('OUTPUT/'+fileName+'_thief.csv',
#                table, delimiter='\t', fmt='%.2f')

if __name__ == '__main__':
    main()
