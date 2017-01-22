from random import randrange, random, gauss
import numpy as np
import sys
import TrNM as hg
import os
import LPsolMerge as lp
import csv

def ReadNucTable(inputStr, headerFlag=None):
    '''
    Takes table in csv format, skips first line (headerFlag) if needed and creates a table of float numbers
    '''
    if headerFlag is None:
        headerFlag = False
    points = []
    with open(inputStr, 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        if (headerFlag):
            print 'skip header'
            csvreader.next()
        for row in csvreader:
            points.append(map(float, row))
    return np.array(points)

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


def main(fileName='default.out', missPenalty=30):
    """function to generate a set of points on a plane that are at least "linker" apart
    returns a matrix (location, isPresentFlag), where first coordinate is location, second - boolen True if the point is present

    """
    numTracks = 0
    flag1 = False
    flag2 = False
    flag3 = False
    flag4 = False
    flag5 = False
    flag6 = False
    flag7 = False
    flag8 = False
    flag9 = False
    flag10 = False
    flag11 = False
    flag12 = False
    flag13 = False
    flag14 = False
    flag15 = False
    flag16 = False
    flag17 = False
    flag18 = False
    flag19 = False
    flag20 = False
    flag21 = False
    flag22 = False
    flag23 = False
    flag24 = False
    flag25 = False
    flag26 = False
    flag27 = False
    flag28 = False
    flag29 = False
    flag30 = False
    flag31 = False
    if len(sys.argv) > 1:
        fileName = sys.argv[1]
    else:
        print ">python process.py <output file> <track 1> <track 2> [<track 3>] [<track 4>] [<track 5>]"
        return 0
    if len(sys.argv) > 2:
        fileNameIn1 = sys.argv[2]
        flag1 = True
    else:
        print "Not enough input files to create tracks..."
        print ">python process.py <output file> <track 1> <track 2> [<track 3>] [<track 4>] [<track 5>]"
        return 0 
    if len(sys.argv) > 3:
        fileNameIn2 = sys.argv[3]
        numTracks = numTracks + 1
        flag2 = True
    else:
        print "Only one input file present, just use it as your track :) ..."
        return 0 
    if len(sys.argv) > 4:
        fileNameIn3 = sys.argv[4]
        numTracks = numTracks + 1
        flag3 = True

    if len(sys.argv) > 5:
        fileNameIn4 = sys.argv[5]
        numTracks = numTracks + 1
        flag4 = True

    if len(sys.argv) > 6:
        fileNameIn5 = sys.argv[6]
        numTracks = numTracks + 1
        flag5 = True 

    if len(sys.argv) > 7:
        fileNameIn6 = sys.argv[7]
        numTracks = numTracks + 1
        flag6 = True

    if len(sys.argv) > 8:
        fileNameIn7 = sys.argv[8]
        numTracks = numTracks + 1
        flag7 = True  

    if len(sys.argv) > 9:
        fileNameIn8 = sys.argv[9]
        numTracks = numTracks + 1
        flag8 = True

    if len(sys.argv) > 10:
        fileNameIn9 = sys.argv[10]
        numTracks = numTracks + 1
        flag9 = True  

    if len(sys.argv) > 11:
        fileNameIn10 = sys.argv[11]
        numTracks = numTracks + 1
        flag10 = True  

    if len(sys.argv) > 12:
        fileNameIn11 = sys.argv[12]
        numTracks = numTracks + 1
        flag11 = True

    if len(sys.argv) > 13:
        fileNameIn12 = sys.argv[13]
        numTracks = numTracks + 1
        flag12 = True  
    
    if len(sys.argv) > 14:
        fileNameIn13 = sys.argv[14]
        numTracks = numTracks + 1
        flag13 = True

    if len(sys.argv) > 15:
        fileNameIn14 = sys.argv[15]
        numTracks = numTracks + 1
        flag14 = True

    if len(sys.argv) > 16:
        fileNameIn15 = sys.argv[16]
        numTracks = numTracks + 1
        flag15 = True

    if len(sys.argv) > 17:
        fileNameIn16 = sys.argv[17]
        numTracks = numTracks + 1
        flag16 = True

    if len(sys.argv) > 18:
        fileNameIn17 = sys.argv[18]
        numTracks = numTracks + 1
        flag17 = True

    if len(sys.argv) > 19:
        fileNameIn18 = sys.argv[19]
        numTracks = numTracks + 1
        flag18 = True 

    if len(sys.argv) > 20:
        fileNameIn19 = sys.argv[20]
        numTracks = numTracks + 1
        flag19 = True 

    if len(sys.argv) > 21:
        fileNameIn20 = sys.argv[21]
        numTracks = numTracks + 1
        flag20 = True 

    if len(sys.argv) > 22:
        fileNameIn21 = sys.argv[22]
        numTracks = numTracks + 1
        flag21 = True 

    if len(sys.argv) > 23:
        fileNameIn22 = sys.argv[23]
        numTracks = numTracks + 1
        flag22 = True 

    if len(sys.argv) > 24:
        fileNameIn23 = sys.argv[24]
        numTracks = numTracks + 1
        flag23 = True 

    if len(sys.argv) > 25:
        fileNameIn24 = sys.argv[25]
        numTracks = numTracks + 1
        flag24 = True

    if len(sys.argv) > 26:
        fileNameIn25 = sys.argv[26]
        numTracks = numTracks + 1
        flag25 = True

    if len(sys.argv) > 27:
        fileNameIn26 = sys.argv[27]
        numTracks = numTracks + 1
        flag26 = True

    if len(sys.argv) > 28:
        fileNameIn27 = sys.argv[28]
        numTracks = numTracks + 1
        flag27 = True

    if len(sys.argv) > 29:
        fileNameIn28 = sys.argv[29]
        numTracks = numTracks + 1
        flag28 = True

    if len(sys.argv) > 30:
        fileNameIn29 = sys.argv[30]
        numTracks = numTracks + 1
        flag29 = True

    if len(sys.argv) > 31:
        fileNameIn30 = sys.argv[31]
        numTracks = numTracks + 1
        flag30 = True

    if len(sys.argv) > 32:
        fileNameIn31 = sys.argv[32]
        numTracks = numTracks + 1
        flag31 = True

    node0 = ReadNucTable(fileNameIn1)
    numAttributes = len(node0[0])
    node1 = ReadNucTable(fileNameIn2)
    if (flag3):
        node2 = ReadNucTable(fileNameIn3)
    if (flag4):
        node3 = ReadNucTable(fileNameIn4)
    if (flag5):
        node4 = ReadNucTable(fileNameIn5)
    if (flag6):
        node5 = ReadNucTable(fileNameIn6)
    if (flag7):
        node6 = ReadNucTable(fileNameIn7)
    if (flag8):
        node7 = ReadNucTable(fileNameIn8)
    if (flag9):
        node8 = ReadNucTable(fileNameIn9)
    if (flag10):
        node9 = ReadNucTable(fileNameIn10)
    if (flag11):
        node10 = ReadNucTable(fileNameIn11)
    if (flag12):
        node11 = ReadNucTable(fileNameIn12)
    if (flag13):
        node12 = ReadNucTable(fileNameIn13)
    if (flag14):
        node13 = ReadNucTable(fileNameIn14)
    if (flag15):
        node14 = ReadNucTable(fileNameIn15)
    if (flag16):
        node15 = ReadNucTable(fileNameIn16)
    if (flag17):
        node16 = ReadNucTable(fileNameIn17) 
    if (flag18):
        node17 = ReadNucTable(fileNameIn18) 
    if (flag19):
        node18 = ReadNucTable(fileNameIn19)
    if (flag20):
        node19 = ReadNucTable(fileNameIn20) 
    if (flag21):
        node20 = ReadNucTable(fileNameIn21) 
    if (flag22):
        node21 = ReadNucTable(fileNameIn22)
    if (flag23):
        node22 = ReadNucTable(fileNameIn23)
    if (flag24):
        node23 = ReadNucTable(fileNameIn24)
    if (flag25):
        node24 = ReadNucTable(fileNameIn25)
    if (flag26):
        node25 = ReadNucTable(fileNameIn26)
    if (flag27):
        node26 = ReadNucTable(fileNameIn27)
    if (flag28):
        node27 = ReadNucTable(fileNameIn28)
    if (flag29):
        node28 = ReadNucTable(fileNameIn29)
    if (flag30):
        node29 = ReadNucTable(fileNameIn30)
    if (flag31):
        node30 = ReadNucTable(fileNameIn31)
    #modifications to add more layers go here
    
##this part for building solution
##NB!!! It requires gpsol installed (no check for this) 
   
    layer0 = hg.layerOfNodes(node0)
    layer1 = hg.layerOfNodes(node1)
    if (flag3):
        layer2 = hg.layerOfNodes(node2)
    if (flag4):
        layer3 = hg.layerOfNodes(node3)
    if (flag5):
        layer4 = hg.layerOfNodes(node4)
    if (flag6):
        layer5 = hg.layerOfNodes(node5)
    if (flag7):
        layer6 = hg.layerOfNodes(node6)
    if (flag8):
        layer7 = hg.layerOfNodes(node7)
    if (flag9):
        layer8 = hg.layerOfNodes(node8)
    if (flag10):
        layer9 = hg.layerOfNodes(node9)
    if (flag11):
        layer10 = hg.layerOfNodes(node10)
    if (flag12):
        layer11 = hg.layerOfNodes(node11)
    if (flag13):
        layer12 = hg.layerOfNodes(node12)
    if (flag14):
        layer13 = hg.layerOfNodes(node13)
    if (flag15):
        layer14 = hg.layerOfNodes(node14)
    if (flag16):
        layer15 = hg.layerOfNodes(node15)
    if (flag17):
        layer16 = hg.layerOfNodes(node16)
    if (flag18):
        layer17 = hg.layerOfNodes(node17)
    if (flag19):
        layer18 = hg.layerOfNodes(node18)
    if (flag20):
        layer19 = hg.layerOfNodes(node19)
    if (flag21):
        layer20 = hg.layerOfNodes(node20)
    if (flag22):
        layer21 = hg.layerOfNodes(node21)
    if (flag23):
        layer22 = hg.layerOfNodes(node22)
    if (flag24):
        layer23 = hg.layerOfNodes(node23)
    if (flag25):
        layer24 = hg.layerOfNodes(node24)
    if (flag26):
        layer25 = hg.layerOfNodes(node25)
    if (flag27):
        layer26 = hg.layerOfNodes(node26)
    if (flag28):
        layer27 = hg.layerOfNodes(node27)
    if (flag29):
        layer28 = hg.layerOfNodes(node28)
    if (flag30):
        layer29 = hg.layerOfNodes(node29)
    if (flag31):
        layer30 = hg.layerOfNodes(node30)
    #modifications to add more layers go here
    print 'layers processed'
    graph = hg.hyperGraph(layer0)
    graph.missPenalty = missPenalty
    graph.okil = missPenalty*2
    print "layer0 done"
    graph.AddLayer(layer1)
    print "layer1 done"
    if (flag3):
        graph.AddLayer(layer2)
        print "layer2 done"
    if (flag4):
        graph.AddLayer(layer3)
        print "layer3 done"
    if (flag5):
        graph.AddLayer(layer4)
        print "layer4 done"
    if (flag6):
        graph.AddLayer(layer5)
        print "layer5 done"
    if (flag7):
        graph.AddLayer(layer6)
        print "layer6 done"
    if (flag8):
        graph.AddLayer(layer7)
        print "layer7 done"
    if (flag9):
        graph.AddLayer(layer8)
        print "layer8 done"
    if (flag10):
        graph.AddLayer(layer9)
        print "layer9 done"
    if (flag11):
        graph.AddLayer(layer10)
        print "layer10 done"
    if (flag12):
        graph.AddLayer(layer11)
        print "layer11 done"
    if (flag13):
        graph.AddLayer(layer12)
        print "layer12 done"
    if (flag14):
        graph.AddLayer(layer13)
        print "layer13 done"
    if (flag15):
        graph.AddLayer(layer14)
        print "layer14 done"
    if (flag16):
        graph.AddLayer(layer15)
        print "layer15 done"
    if (flag17):
        graph.AddLayer(layer16)
        print "layer16 done"
    if (flag18):
        graph.AddLayer(layer17)
        print "layer17 done"
    if (flag19):
        graph.AddLayer(layer18)
        print "layer18 done"
    if (flag20):
        graph.AddLayer(layer19)
        print "layer19 done"
    if (flag21):
        graph.AddLayer(layer20)
        print "layer20 done"
    if (flag22):
        graph.AddLayer(layer21)
        print "layer21 done"
    if (flag23):
        graph.AddLayer(layer22)
        print "layer22 done"
    if (flag24):
        graph.AddLayer(layer23)
        print "layer23 done"
    if (flag25):
        graph.AddLayer(layer24)
        print "layer24 done"
    if (flag26):
        graph.AddLayer(layer25)
        print "layer25 done"
    if (flag27):
        graph.AddLayer(layer26)
        print "layer26 done"
    if (flag28):
        graph.AddLayer(layer27)
        print "layer27 done"
    if (flag29):
        graph.AddLayer(layer28)
        print "layer28 done"
    if (flag30):
        graph.AddLayer(layer29)
        print "layer29 done"
    if (flag31):
        graph.AddLayer(layer30)
        print "layer30 done"
    #modifications to add more layers go here

    graph.EdgeCostComputation()
    print 'done building graph'
    hg.CPLEXprint(graph, fileName+'_tmp.lp')
    # run linear solver gpsol and parse it's output
    print 'start linear solver'
    os.system('./runLS.sh '+fileName+'_tmp.lp')
    print 'linear solution done'
    lpSol = lp.ReadColumn(fileName+'_tmp.csv')
    print 'solution read'
    os.system('rm '+fileName+'_tmp.lp')
    os.system('rm '+fileName+'_tmp.csv')
    os.system('rm '+fileName+'_tmp.sol')
    table = graph.GetTrackStat(lpSol, numAttributes)
    print table
    np.savetxt(fileName,
               table, delimiter='\t', fmt='%.2f')
    print "Finally Done"
    return 1

if __name__ == '__main__':
    main()
