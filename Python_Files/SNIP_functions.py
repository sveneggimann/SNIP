# ======================================================================================
# Copyright 2014  Swiss Federal Institute of Aquatic Science and Technology
#
# This file is part of SNIP (Sustainable Network Infrastructure Planning)
# SNIP is used for determining the optimal degree of centralization for waste
# water infrastructures. You find detailed information about SNIP in Eggimann et al. (2014).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#    
# The Djikstra and a* algorithm are adapted from Hetland (2010).
# The algorithm is developed for Python 2.7 and ArcGIS 10.2
#
# Literature
# ----------
# Eggimann Sven, Truffer Bernhard, Maurer Max (2015): To connect or not to connect? 
# Modelling the optimal degree of centralisation for wastewater infrastructures.   
# Water Research, 84, 218-231.
#
# Hetland M.L. (2010): Python Algorithms. Mastering Basic Algorithms in the Python Language. apress.
#
# Contact:   sven.eggimann@eawag.ch
# Version    1.0
# Date:      1.07.2015
# Autor:     Eggimann Sven
# ======================================================================================

# Imports
import arcpy
import os, math, sys, operator
from datetime import datetime
from SNIP_astar import *
from SNIP_costs import *


def distanceCalc2d(p0, p1):
    """
    This functions calculates the euclidian distance in 2d space.

    Input Arguments: 
    p0, p1                --    Points with coordinates. Form: (ID, X, Y, Z)
    
    Output Arguments:
    distance              --    Distance between the points not taking in account the height difference.
    """
    distance = math.hypot(p0[0] - p1[0], p0[1] - p1[1])         
    return distance

def distanceCalc3d(p0, p1):
    """
    This functions calculates the euclidian distance in 3d space.

    Input Arguments: 
    p0, p1                --    Two points
    
    Output Arguments:
    distance              --    Distance between the points
    slope                 --    Slope between the two points
    heightDiff            --    Height difference betwen the two points
    """
    xp0, xp1, yp0, yp1, zFrom, zTo = p0[0], p1[0], p0[1], p1[1], p0[2], p1[2]
    distancePlanar = math.hypot(xp0 - xp1, yp0 - yp1) 

    if distancePlanar == 0:
        raise Exception("ERROR: DISTANCE TO ITSELF is calculated")
        return 0, 0, 0  # Distance zero is returned

    heightDiff = zTo - zFrom
    distance = math.sqrt(pow(distancePlanar, 2) + pow(heightDiff, 2))
    slope = (float(heightDiff) / float(distancePlanar)) * 100           # Slope in %
    return distance, slope, heightDiff

def fastCopyNodes(nodes):
    '''
    This function speeds up the copying of the nodes (instead of a deep copy)
    '''
    copyList = []
    for i in nodes:
        l = []
        for e in i[9]:
            l.append(e)
        copyList.append([i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], l, i[10]])
    return copyList

def fastCopy(nodes):
    '''
    This function speeds up the copying.
    '''
    copyList = []
    for i in nodes:
        l = []
        for f in i:
            l.append(f)
        copyList.append(l)
    return copyList

def getCornerPoints(aggregatedNodes):
    """
    This function reads out the top left and bottom right coordinate of points in a list. Is used for average nearest neighbour distance
    
    Input Arguments:
    aggregatedNodes    --      List with aggregated points
    
    Output Arguments:
    startnode          --      Return top left ID coordinate of tile with highest density
    """
    xRight, yTop, xLeft, yBottom = 0, 0, 9999999999, 9999999999
    for i in aggregatedNodes:
        if i[0] < xLeft:
            xLeft = i[0] 
        if i[0] > xRight:
            xRight = i[0]
        if i[1] > yTop:
            yTop = i[1] 
        if i[1] < yBottom:
            yBottom = i[1] 
    tupleTopLef, tupleBottomRight = (xLeft, yTop), (xRight, yBottom)
    return tupleTopLef, tupleBottomRight

def getTopFleftBottomRightTuple(aggregatedNodes):
    """
    This function reads out the top left and bottom right coordinate of points in a list.
    
    Input Arguments:
    aggregatedNodes    --      List with aggregated points
    
    Output Arguments:
    startnode          --      Return top left ID coordinate of tile with highest density
    """
    xRight, yTop, xLeft, yBottom = 0, 0, 9999999999, 9999999999
    for i in aggregatedNodes:
        if i[1] < xLeft:
            xLeft = i[1] 
        if i[1] > xRight:
            xRight = i[1]
        if i[2] > yTop:
            yTop = i[2] 
        if i[2] < yBottom:
            yBottom = i[2] 
    tupleTopLef, tupleBottomRight = (xLeft, yTop), (xRight, yBottom)
    return tupleTopLef, tupleBottomRight

def assignHighAggregatedNodes(aggregatetPoints, rasterPoints, rasterSize, minTD):
    """
    This function assigns the correct height to aggregate nodes by searching the height of a DEM.
    
    Input Arguments: 
    aggregatetPoints      -    nodes
    rasterPoints          -    list with nodes to reset flow
    rasterSize            -    WWTP from which the breath search was made ??? Really Needed?
    minTD                 -    Minimum Trench Depth
    
    Output Arguments:
    withHeith             -    Aggregated nodes with correct height
    """
    withHeith = []
    for i in aggregatetPoints:
        X, Y = i[1], i[2]
        for f in rasterPoints:  # Assign height 
            if (f[1] <= X + rasterSize and f[1] >= X - rasterSize) and (f[2] <= Y + rasterSize and f[2] >= Y - rasterSize):
                heightSTART = f[3]
                break
        z = [i[0], i[1], i[2], heightSTART, i[4], i[5], i[6], i[7], i[8], i[9], heightSTART - minTD]
        withHeith.append(z)
    return withHeith


def getFlowInitial(nodesFromNetwork, nodesToNetwork, nodes, flowWWFrom, flowWWTO, notInaNetworkBeforeConnection):
    """
    This function calculates the initial flows of the networks where the wwpts need to be connected.  
    
    Input Arguments: 
    nodesFromNetwork                --    Nodes of fromNetwork
    nodesToNetwork                  --    Nodes of toNetwork
    nodes                           --    Nodes
    flowWWFrom                      --    Flow of wwtpFrom
    flowWWTO                        --    Flow of wwtpto
    notInaNetworkBeforeConnection   --    Nodes not in a network before connection

    Output Arguments:
    flowInitial_from                --    Initial flow of wwtp in network of from node
    flowInitial_to                  --    Initial flow of wwtp in network of to node
    """
    if len(nodesFromNetwork) > 1:
        _, _, intialFlowA, intialFlowB, _, _ = getPns(nodesFromNetwork[1], nodes)
        flowInitial_from = flowWWFrom - (intialFlowA + intialFlowB)  
    else:
        flowInitial_from = flowWWFrom 
    
    if len(nodesToNetwork) > 1:
        if nodesToNetwork[1] in notInaNetworkBeforeConnection:
            _, _, intialFlowA, _, _, _ = getPns(nodesToNetwork[0], nodes)
            flowInitial_to = flowWWTO - intialFlowA
        else:
            _, _, intialFlowA, intialFlowB, _, _ = getPns(nodesToNetwork[1], nodes)
            flowInitial_to = flowWWTO - (intialFlowA + intialFlowB)
    else:
        flowInitial_to = flowWWTO 
    return flowInitial_from, flowInitial_to

def breathSearch(ID, sewers):
    """
    This function searches all nodes flowing to a WWTP (a breath search of a graph). 
    
    Input Arguments: 
    ID                  --    ID
    sewers              --    Sewer Network
    
    Output Arguments:
    allNodesToDelet     --    All nodes flowing to a wwtp
    """
    initialnetWorkToRemove, allNodesToDelet, newScrapList = ID, [], [ID]

    while 1:
        scrapList = newScrapList
        newScrapList = []
        # Search all nodes flowing to this wastewater
        for foundNode in scrapList:
            for ID in sewers:
                if sewers[ID] != ():
                    if sewers[ID][0] == foundNode:
                        newScrapList.append(ID)
                        allNodesToDelet.append(ID)
        if len(newScrapList) == 0:
            break    
    allNodesToDelet.append(initialnetWorkToRemove)  # Add initial WWTP
    return allNodesToDelet

def addBuildingsFarFromRoadTo(aggregatetPoints, forSNIP):
    '''Add all buildings which are not part of an aggregated node to the nodes
    Input Arguments:
    aggregatetPoints       --      List with aggregated nodes
    forSNIP                --      list for SNIP Algorithm.     
    
    Output Arguments:
    forSNIP                --      list for SNIP Algorithm now added with the nodes not aggregated on the road
    '''
    for i in aggregatetPoints:
        alreadyInList = 0
        for e in forSNIP:
            if e[1] == i[1] and e[2] == i[2]:
                alreadyInList = 1
                break
        if alreadyInList == 0:
            z = [i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10]]
            forSNIP.append(z)
    return forSNIP

def densityBasedSelection(aggregatedNodes, tileSize):
    """
    Find top left node and divide space into squares. Count how many aggregated nodes are within each square.
    Then it selects randomly one node within the square with the most nodes.
    
    Input Arguments:
    aggregatedNodes   --    Aggregated points
    tileSize          --    Size of a tile [m]
    
    Output Arguments:
    squares           --    Raster with density
    startnode         --    Return top left ID coordinate of tile with highest density
    startX, startY    --    Coordinates of startnode
    """
    tupleTopLef, tupleBottomRight = getTopFleftBottomRightTuple(aggregatedNodes)
    xRight, yTop, xLeft, yBottom = tupleBottomRight[0], tupleTopLef[1], tupleTopLef[0], tupleBottomRight[1]

    width = xRight - xLeft  # in meters
    height = yTop - yBottom  # in meters
    
    nrOfTilesHorizontal = range(int(float(width) / float(tileSize)) + 1)
    nrOfTilesVertical = range(int(float(height) / float(tileSize)) + 1)

    squares, line, IdNr = {} , -1, -1  # Dictionary to store squres { sqrID: (x, y)},  Iterators

    for tile in nrOfTilesVertical:
        line += 1
        row = -1
        for entry in nrOfTilesHorizontal:
            IdNr += 1
            row += 1
            squares[IdNr] = [xLeft + row * tileSize, yTop - line * tileSize, 0]
    
    # check how many are in each tile
    for entry in squares:
        nrOfPnts = 0
        for i in aggregatedNodes:
            if i[1] >= squares[entry][0] and i[1] < (squares[entry][0] + tileSize) and i[2] > (squares[entry][1] - tileSize) and i[2] <= squares[entry][1]:
                nrOfPnts += 1
        squares[entry][2] = nrOfPnts
    
    # Get coordinates with highest density
    maxNrofTiles = 0
    for entry in squares:
        if squares[entry][2] > maxNrofTiles:
            topLefTileX, topLefTileY, maxNrofTiles = squares[entry][0], squares[entry][1], squares[entry][2]
    
    # Iterate list and select first point which is in highest density tile
    for i in aggregatedNodes:
        if i[1] >= topLefTileX and i[1] <= topLefTileX + tileSize and i[2] >= topLefTileY - tileSize and i[2] <= topLefTileY:
            startnode, startX, startY = i[0], i[1], i[2]
            break
    return squares, startnode, startX, startY       

def clearFlow(nodes, allNodesNotInNet):
    """
    This function sets the flow to zero in all nodes which are not in a network.
    
    Input Arguments: 
    nodes              --   List with all nodes, all nodes which are not in a network
    allNodesNotInNet   --   All nodes which are not in a network
    
    Output Arguments:
    nodes              --   List with all nodes with corrected flow 
    """
    for nodeID in allNodesNotInNet:
        pos = 0
        for i in nodes:
            if i[0] == nodeID:  
                break
            pos += 1
        nodes[pos][4] = 0    # Clear flow
    return nodes

def getNotInNetwork(pathBetweenWWTPs, network):
    """
    This functions searches for all nodes of a path which are not already part of an existing network.
    
    Input Arguments: 
    pathBetweenWWTPs    --    path
    network             --    existing network
    
    Output Arguments:
    nodesNotInNetwork   --    all nodes of path not in a network
    """
    nodesNotInNetwork = []
    for i in pathBetweenWWTPs:
        try:
            _ = network[i]
        except KeyError:
            nodesNotInNetwork.append(i)
    return nodesNotInNetwork
        
def checkIfIsPump(pumps, nodeID): 
    """
    Check if there is a pump at a node.
    
    Input Arguments: 
    nodeID           --    node ID
    pumps            --    list with pumps
    
    Output Arguments:
    isPump           --    criteria (0: There is no pump, 1: There is a pump)
    """
    isPump = False
    for pump in pumps:
        if pump[0] == nodeID:
            isPump = True
            break
    return isPump
                    
def InvertandswapID(listToInvert):
    """
    Invert a list and swap Id
    
    Input Arguments: 
    listToInvert           --    list to swap
    
    Output Arguments:
    invertList             --    swapped List
    """
    invertList = []
    for i in listToInvert:
        invertList.append([i[1][0], [i[0], i[1][1]]])   
    invertList = invertList[::-1]                       # Invert     
    return invertList
                 
def readOnlyNetwork(nodes, network):
    """
    Get only the nodes which are part of a network.
    
    Input Arguments: 
    nodes           --    nodes 
    network         --    Network
    
    Output Arguments:
    flowNodes       --    All Nodes of Network with flow
    cleanNetwork    --    Networks with only nodes where there is a flow.
    """
    flowNodes, cleanNetwork = [], {}

    for i in network:  
        for z in nodes:                             
            if z[0] == i:                            
                x, y, height, flow, forceKriteria, nodeFlow, trench = z[1], z[2], z[3], z[4], z[5], z[8], z[3] - z[10]      
                a = [i, x, y, height, flow, forceKriteria, nodeFlow, trench]
                break
        if flow > 0 and nodeFlow > 0:
            flowNodes.append(a)
            cleanNetwork[i] = network[i]
        if flow == 0 and nodeFlow > 0:
            flowNodes.append(a)
            cleanNetwork[i] = network[i]      
        if flow > 0 and nodeFlow == 0:
            flowNodes.append(a)
            cleanNetwork[i] = network[i]
    return cleanNetwork, flowNodes

def WWTPsToDraw(WWTPs, nodes):
    """
    Append coordinates to list with wwtps.
    
    Input Arguments: 
    WWTPs                  --    list with wwtps 
    nodes                  --    list with nodes
    
    Output Arguments:
    wtpstodraw             --    List with wwtps to draw.
    """
    wtpstodraw = []
    for i in WWTPs:
        for z in nodes:                                       # Get all infos from nodes
            if z[0] == i[0]:                                  # Read out correct list element
                wtpstodraw.append([i[0], i[1], z[1], z[2]])   # Append to wwtp list
                break    
    return wtpstodraw

def removePumpsWhereNoNetwork(pumps, allNodesToDelet):
    """
    Remove pumps of list with nodes which are removed.
    
    Input Arguments: 
    pumps              --    list with pumps
    allNodesToDelet    --    Nodes which are removed
    
    Output Arguments:
    pumps              --    pump list with removed IDs
    """
    for i in allNodesToDelet:
        position = 0
        for e in pumps:
            if e[0] == i:
                del pumps[position]
                break
            position += 1
    return pumps

def removePumpCheck(pumps, ID):
    """
    Remove pumps with certain ID.
    
    Input Arguments: 
    pumps                 --    list with pumps
    ID                    --    ID which is removed
    
    Output Arguments:
    pumps                 --    Pump list with removed ID
    """
    delPos = 0              
    for i in pumps:
        if i[0] == ID:
            del pumps[delPos]
            break
        delPos += 1
    return pumps
            
def removeInSewers(sewers, liste):
    """
    Clean sewers of nodes from a list.
    
    Input Arguments: 
    sewers             --    dictionary
    liste              --    list with elements to remove
    pathBetweenWWTPs   --    Path between WWTPS
    
    Output Arguments:
    sewers             --    dictionary with removed elements
    """
    for i in liste:
        if i in sewers:
            del sewers[i]
    return sewers
                              
def changeFlowInNode(nodeList, nodeID, newFlow):
    """
    Change flow in node.
    
    Input Arguments: 
    nodeList            --   list with all nodes
    nodeID              --   Id of node to change flow
    newFlow             --   new flow
    
    Output Arguments:
    nodeList            --   list with changed flow
    """
    for i in nodeList:
        if i[0] == nodeID:
            i[4] = newFlow
            return nodeList

def readPath(path, ID):
    """
    This functions adds an ID to a path and returns an inverted path.
    
    Input Arguments: 
    path               --   path
    ID                 --   Id
    
    Output Arguments:
    newPath            --   appended new path
    """
    newPath = []
    for i in path:
        newNode = i[0]
        newPath.append(newNode)
    newPath.append(ID)
    newPath = newPath[::-1]  # inverse
    return newPath

def getclosestWTP(ID, WWTPs, sewers):
    """
    Search for closest wwpt in Network
    
    Input Arguments: 
    ID                       --   ID
    WWTPs                    --   list with WWTPs
    sewers                   --   Network
    
    Output Arguments:
    waytoclosestWTP          --   path to nearest wwtp
    closestARAtraditionell   --   closest wwtp
    """
    waytoclosestWTP, element, hatARA = [], ID, 0

    while 1:
        waytoclosestWTP.append(element)
        if sewers[element][0] != ():
            nextElement = sewers[element][0]
        else:
            if len(waytoclosestWTP) == 1:
                nextElement = element
            break
        for i in WWTPs:
            if i[0] == nextElement:
                waytoclosestWTP.append(nextElement)  
                hatARA = 1
                break
        if hatARA == 1:
            break
        element = nextElement
    closestARAtraditionell = waytoclosestWTP[-1]
    return waytoclosestWTP, closestARAtraditionell  # Return: whole path, closest WWTP

def getclosestWTPreconnection(ID, WWTPs, sewers):
    """
    Search for closest wwtp in network.
    
    Input Arguments: 
    ID                              --    ID
    WWTPs                           --    list with WWTPs
    sewers                          --    Sewer Network
    
    Output Arguments: 
    waytoclosestWWTP      --    path to closest wwtp
    closestARAtraditionell          --    closest wwtp
    """
    waytoclosestWWTP, element = [], ID

    while 1:
        waytoclosestWWTP.append(element)
        
        if sewers[element][0] != ():  
            nextElement = sewers[element][0]
        else:
            waytoclosestWWTP.append(element)
            break
        element = nextElement

    closestARAtraditionell = waytoclosestWWTP[-1]
    return waytoclosestWWTP, closestARAtraditionell  

def readOnlyNodesID(archPathWWTP):
    """
    Read out only the nodeID of a path.
    
    Input Arguments: 
    archPathWWTP            --   Path between wwtps
    
    Output Arguments:
    pathBetweenWWTPs        --   Path between wwtps only with IDs
    """
    pathBetweenWWTPs = []
    for i in archPathWWTP:
        pathBetweenWWTPs.append(i[0])
    pathBetweenWWTPs.append(archPathWWTP[-1][1][0])  
    return pathBetweenWWTPs       
                    
def addToEdgeList(edgeList, distStartEnd, slopeDijkstra, idp0, idp1, p0, p1):
    """
    This function adds an edge in case it is not already existing to a list with the all edges.
    
    Input Arguments: 
    edgeList           --   Path between wwtps
    distStartEnd       --   Distance between start and end
    slopeDijkstra      --   Slope calculated from Djikstra Distance
    idp0, idp1         --   IDs of start and ending node
    p0, p1             --   coordinates of starting and ending node
    
    Output Arguments:
    edgeList           --   List with added edge
    
    """
    alreadyInEdgeList = 0
    for i in edgeList:
        if i[0][0] == idp0 and i[1][0] == idp1 or i[1][0] == idp0 and i[0][0] == idp1:
            alreadyInEdgeList = 1
            break
    if alreadyInEdgeList == 0:
        edgeList.append([[idp0, p0[0], p0[1], p0[2]], [idp1, p1[0], p1[1], p1[2]], distStartEnd, slopeDijkstra, 0, 0])    
    return edgeList
                                          
def costToWTP(pathToWTP, edgesID, nodes, pumps, minTD, toNode, sewers, discountYearsSewers, interestRate, stricklerC, operationCosts, f_SewerCost):
    """
    This function calculates the network sewer costs from a node to the closest wwtp. Calculations start from source.
    
    Input Arguments: 
    pathToWTP             -- path to closest wwtp 
    edgesID               -- edges 
    nodes                 -- nodes
    pumps                 -- list with pumps
    minTD                 -- min trench depth
    toNode                -- source ID
    sewers                -- sewer network
    discountYearsSewers   -- lifespan of sewers
    interestRate          -- real interest rate
    stricklerC            -- strickler coefficient
    operationCosts        -- Operation costs of sewer per meter [CHF]
    
    Output Arguments: 
    totCosts              -- Cost of network (new and already existing)
    edgesID               -- Edge list with adopted pipe diameters
    """
    cnt, totCosts = -1, 0
    flowCurrentNode = getSummedFlow(nodes, toNode)              # Get flow of the new connected node
    
    # Iterate path to WWTP
    for i in pathToWTP:
        cnt += 1      # Initialisation
        nextNode = i  # Initialisation
        
        # Get distance & flow of each segments        
        if cnt > 0:
            pos = 0
            for edge in edgesID:
                if edge[0][0] == oldNode and edge[1][0] == nextNode:                # Stored inverse, thus slope needs to get inverted
                    distanz, slope = edge[2], edge[3] * -1                          # slope needs to be inverted                               
                    break   
                if edge[1][0] == oldNode and edge[0][0] == nextNode:
                    distanz, slope = edge[2], edge[3]                               # distance, # slope stays the same      
                    break
                pos += 1   
            
            flowSegment, trenchDepthFROM = readFlow(nodes, oldNode)                 # get flow
            _, trenchDepthTO = readFlow(nodes, nextNode)                            # get new trenchDepth

            # check if FromNode is a pump. If yes, the Average Trench Depth is minimum Trench Depth (for both ways)
            isaPump = checkIfIsPump(pumps, nextNode)

            if isaPump == True:
                averageTD = minTD
            else:
                averageTD = (abs(trenchDepthFROM) + abs(trenchDepthTO)) / 2
            
            # Calculate pipe segment costs
            totalFlow = flowSegment + flowCurrentNode     
            
            # flow upstream used for estimating diameter
            pipeDiameter = getPipeDiameter(totalFlow, slope, stricklerC)                                                                            # get pipe diameter 
            segmentCosts = calculatePipeCosts(pipeDiameter, distanz, averageTD, discountYearsSewers, interestRate, operationCosts, f_SewerCost)     # segment pipe costs #could be added totalFlow, slope, 
            edgesID[pos][4] = pipeDiameter         # Add pipe diameter to edgeList
     
            # If nextNode not in sewers then this pipe has to be newly constructed with full new costs. If nextNode is in sewers, then this pipe is already constructed and only only partial cost arise
            totCosts = totCosts + segmentCosts
        oldNode = i
    return totCosts, edgesID

def costsBetweenWWTPs(pathToWTP, edgesID, nodes, inverse, pumps, minTD, discountYearsSewers, interestRate, stricklerC, operationCosts, f_SewerCost):
    """
    This function calculates the network sewer costs from a wwtp to another wwtp.  Start from back
    
    Input Arguments: 
    pathToWTP             --    path to closest wwtp 
    edgesID               --    edges 
    nodes                 --    nodes
    inverse               --    Direction criteria
    pumps                 --    list with pumps
    minTD                 --    minium trench depth
    discountYearsSewers   --    lifespan of sewers
    interestRate          --    real interest rate
    stricklerC            --    strickler coefficient
    operationCosts        --    Operation costs of sewer per meter [CHF]
    
    Output Arguments: 
    segmentCostsSUM       --    Cost of pipes
    edgesID               --    Edges ID
    """
    segmentCostsSUM, count = 0, -1
  
    # Path to WWTP
    for i in pathToWTP:
        count += 1      # Initialisation
        nextNode = i    # Get next node

        # Get distance & flow of each segments        
        if count > 0:
            pos = 0
            for edge in edgesID:  # distances from edge list. 
                if edge[0][0] == oldNode and edge[1][0] == nextNode:
                    distanz, slope = edge[2], edge[3] * -1  # distance, slope needs to be multiplyed by minus one as direction is inversed 
                    break
                if edge[1][0] == oldNode and edge[0][0] == nextNode:
                    distanz, slope = edge[2], edge[3]  # distance, slope stays the same# 
                    break
                pos += 1
                
            for punkt in nodes:
                if punkt[0] == oldNode:
                    flowSegment, toflow, trenchDepthFROM = punkt[4], punkt[8], punkt[3] - punkt[10]  # flow from Node up Stream
                    break
                
            # Calculate pipe segment costs
            totalFlow = flowSegment + toflow  
            _, trenchDepthTO = readFlow(nodes, nextNode)  # get new trenchDepth

            # check if FromNode is a pump. If yes, the Average Trench Depth is minimum Trench Depth (for both ways)
            isaPump = checkIfIsPump(pumps, nextNode)
                
            if isaPump == True:
                averageTD = minTD
            else:
                averageTD = (abs(trenchDepthFROM) + abs(trenchDepthTO)) / 2

            pipeDiameter = getPipeDiameter(totalFlow, slope, stricklerC)    
            edgesID[pos][4] = pipeDiameter   # Change pipe diameter      

            # Normal Way
            if inverse == 0:
                segmentCosts = calculatePipeCosts(pipeDiameter, distanz, averageTD, discountYearsSewers, interestRate, operationCosts, f_SewerCost)  
            else:                           # Inverse flow
                slope = slope * -1          # slope is inversed
                segmentCosts = calculatePipeCosts(pipeDiameter, distanz, averageTD, discountYearsSewers, interestRate, operationCosts, f_SewerCost)
            segmentCostsSUM = segmentCostsSUM + segmentCosts
        oldNode = i
        
    return segmentCostsSUM, edgesID

def appendDistances(pathtoOrigin, sewers):
    """
    This functions appends the distances and slope to a path 
    
    Input Arguments: 
    pathtoOrigin         --   Path between WWTPs
    sewers               --   Sewer Network
    
    Output Arguments:
    pathOrigin           --   Path with distance and slope
    """
    pathOrigin = []                     # list with nods and distances to WTP
    for to_node in pathtoOrigin:
        if sewers[to_node] is ():
            z = (0)                     # needed as the startnode key has no distance
        else:
            z = sewers[to_node][1]      # distance between nodes     
        pathOrigin.append([to_node, z])
    return pathOrigin
    
def appendListWTPsToNetwork(sewers, liste):
    """
    Add all WWTPs in in a list to the network.
    
    Input Arguments: 
    sewers            --   Sewer Network
    WWTPliste         --   list with wwtps
    
    Output Arguments:
    sewers            --   Sewers with added wwtps
    """
    for entry in liste:  
        sewers[entry] = (), 0
    return sewers

def getPopNodes(pathNearWTP, nodes, sewersNoC):
    """
    Gets all populated, unconnected nodes on the new path to be added to the network.
    
    Input Arguments: 
    pathNearWTP         --    Path to source
    nodes               --    nodes
    sewersNoC           --    network without added path
    
    Output Arguments:
    scrapPath          --    Path with distance and slope
    sumFlow             --    Summed flow along the archPath
    """
    scrapPath, sumFlow = [pathNearWTP[0][0]], 0
        
    # Check if on the way to the closest WWTP there are populated nodes with no flow
    for i in pathNearWTP:
        for entry in nodes:
            if entry[0] == i[0]:
                idid, flow = entry[0], entry[8]
                break                      
        AlreadyConnected = 0
                          
        # Check if node has been connected in this iteration
        if idid in sewersNoC:
            AlreadyConnected = 1
                                  
        # The the node is not connected, populated and has no flow, put into ScrapListL in order that the flow can iteratively be added.
        if flow > 0 and AlreadyConnected == 0:
            scrapPath.append(entry[0])
            sumFlow = sumFlow + flow
            
    # List with nodes from which to update flow         
    scrapPath = scrapPath[::-1]  # Invert scrapPath
    return scrapPath, sumFlow

def readFlow(nodes, ID):
    """
    Read flow of a node with ID.
    
    Input Arguments: 
    nodes              --    nodes
    ID                 --    ID
    
    Output Arguments:
    flowSegment        --    Flow
    trenchDepth        --    Trench Depth
    """
    for punkt in nodes:
        if punkt[0] == ID:
            flowSegment = punkt[4]              # flow from Node up Stream
            trenchDepth = punkt[3] - punkt[10]  # Trench depth
            return flowSegment, trenchDepth

def getPns(ID, nodes):
    """
    This function reads out the coordinates of a node.
    
    Input Arguments: 
    ID               --    ID
    nodes            --    list with the nodes
    
    Output Arguments:
    idNode           --    Flow
    coordinates      --    (X, Y, Z)
    flowToNode       --    Flow from node itself
    flowInNode       --    Flow flowing in node
    foreConCrit      --    Criteria whether connection is needed or not
    """
    for i in nodes:
        if i[0] == ID:                                            
            coordinates, flowInNode, foreConCrit, flowToNode, trenchHight =(i[1], i[2], i[3]), i[4], i[5], i[8], i[10]
            return ID, coordinates, flowToNode, flowInNode, foreConCrit, trenchHight

def appendToNetwork(sewers, liste):
    """
    Append list with nodes to sewer network.
    
    Input Arguments: 
    sewers            --    sewers
    liste             --    nodes
    
    Output Arguments:
    sewers           --    Network with added nodes.
    """
    for e in liste:
        if e[1][0] not in sewers:                                          
            sewers[e[1][0]] = e[0], e[1][1] 
    return sewers

def insertPathDirection(sewers, liste):
    """
    Check if path already exists in sewers and if no, add path to network.
    
    Input Arguments: 
    sewers            -    Network
    liste             -    nodes
    
    Output Arguments:
    sewers            -    Network with added path.
    """
    for entry in liste:
        sewers[entry[0]] = entry[1][0], entry[1][1]
    sewers[liste[-1][1][0]] = (), 0 
    sewers[liste[0][0]] = liste[0][1][0], liste[-0][1][1]
    return sewers

def appendToSewers(sewers, path):
    """
    Input Arguments: 
    sewers               -    Sewer Network
    path                 -    path
    
    Output Arguments:
    sewers               -    Network with added path.
    """
    for entry in path:  # insert the whole path to sewers
        if entry not in sewers: 
            sewers.append(entry)
    return sewers

def delFromSewers(sewers, path):
    """
    Remove entries in Network
    
    Input Arguments: 
    sewers            -    Network
    path              -    Path
    
    Output Arguments:
    sewers            -    Network with removed path.
    """
    for entry in path:  
        if entry in sewers: 
            del sewers[entry]
    return sewers
                 
def testExpansion(PN):
    '''
    This function tests wheter the EM still needs to be executed.
    
    Input:
    PN           -    Length of prim edges
    
    Output:
    expansion    -    Criteria
    '''
    if len(PN) > 0:    
        expansion = 1       # Return to expansion-module
    else:
        expansion = 0       # All nodes are connected 
    return expansion
        
def invertArchPath(archPathWWTP):
    """
    Invert a path in a list.
    
    Input Arguments: 
    archPathWWTP         -    Path between new connected node and network
    
    Output Arguments:
    invertedArchPath     -    Inverted path
    """
    invertedArchPath = []
    for i in archPathWWTP:
        invertedArchPath.append([i[1][0], [i[0], i[1][1]]])
    invertedArchPath = invertedArchPath[::-1]
    return invertedArchPath

def addToListWTPs(archPath, nodes, sewersNoC, WWTPs):
    """
    Search elements in a path and if a node is a wwtp, ad it to the list of wwtps.
    
    Input Arguments: 
    archPath                --    path
    nodes                   --    nodes
    sewersNoC               --    Network before adding path
    WWTPs                   --    list with wwtps
    
    Output Arguments:
    WWTPs                   --    Updated list with wwtps
    intermediateWWTPs       --    list with wwtp ID to add to sewer network 
    """
    intermediateWWTPs = []

    for entry in archPath:
        to = entry[1][0]
        
        # If there is flow itself the node becomes a wwtps
        for i in nodes:
            if i[0] == to:
                flowulate = i[8] + i[4] 
                break
                                  
        if to not in sewersNoC:
            if flowulate > 0:
                WWTPs.append([to, flowulate])
                intermediateWWTPs.append(to)
    return WWTPs, intermediateWWTPs

def addEdgesUpdateStreetNetwork(wayProperty, archPathListMST, nodes, boundingCandidates, edgeList, streetNetwork):
    """
    Add new edges and add the edges to the street network. If the new edges are made out of raster coordinates, search coordinates and add them as well.
    
    Input Arguments: 
    wayProperty             -    Criteria for sewer cost calulation. If a* was used for path calculation, wayproperty == means field, otherwise street.
    archPathListMST         -    new path
    nodes                   -    nodes
    boundingCandidates      -    All DEM-points in a bounding box
    edgeList                -    list with edges
    streetNetwork           -    graph with streets
    
    Output Arguments:
    edgeList                -    Inverted path
    streetNetwork           -    Street network
    """
    # Append edges
    for entry in archPathListMST:
        idp0, idp1, distanz3d, steigung, foundCorrIDp0, foundCorrIDp0 = entry[0], entry[1][0], entry[1][1], entry[1][2], 0, 0
        
        for i in nodes:
            if i[0] == idp0:
                corrXIDp0, corrYIDp0, corrZIDp0, foundCorrIDp0 = i[1], i[2], i[3], 1
                break

        for i in nodes:
            if i[0] == idp1:
                corrXIDp1, corrYIDp1, corrZIDp1, foundCorrIDp1 = i[1], i[2], i[3], 1
                break

        # If coordintaes are rasterCoordinates
        if foundCorrIDp0 == 0:
            for i in boundingCandidates:
                if i[0] == idp0:
                    corrXIDp0, corrYIDp0, corrZIDp0 = i[1], i[2], i[3]
                    break

        if foundCorrIDp1 == 0:
            for i in boundingCandidates:
                if i[0] == idp1:
                    corrXIDp1, corrYIDp1, corrZIDp1 = i[1], i[2], i[3]
                    break

        # Check if already in EdegList
        alreadyInEdgeList = 0
        for i in edgeList:
            if i[0][0] == idp0 and i[1][0] == idp1 or i[1][0] == idp0 and i[0][0] == idp1:
                alreadyInEdgeList = 1
            
                # Adapt edgeList
                pipeDiameter = 0
                i[2], i[3], i[4], i[5] = distanz3d, steigung, pipeDiameter, wayProperty
                break
        
        if alreadyInEdgeList == 0:
            pipeDiameter = 0
            edgeList.append([[idp0, corrXIDp0, corrYIDp0, corrZIDp0], [idp1, corrXIDp1, corrYIDp1, corrZIDp1], distanz3d, steigung, pipeDiameter, wayProperty])  # append new edge  ([p1], [p2], distanc, slope),...)
        
        # Change street graph and enter distances. Check if already in StreetGraph. If not, create {}
        try:
            _ = streetNetwork[idp0]
        except KeyError:
            streetNetwork[idp0] = {}
        try:
            _ = streetNetwork[idp1]    
        except KeyError:
            streetNetwork[idp1] = {}  
            
        # Street is added both ways
        streetNetwork[idp0][idp1] = distanz3d   
        streetNetwork[idp1][idp0] = distanz3d  
    return edgeList, streetNetwork

def getFlowWWTP(WWTPs, ID):
    """
    Read out flow from list with WWTP.
    
    Input Arguments: 
    WWTPs         -    List with wwtps
    ID            -    ID
    
    Output Arguments:
    flowInWWTP    -    flow in WWTP
    """
    for i in WWTPs:
        if i[0] == ID:
            flowInWWTP = i[1]
            return flowInWWTP
    raise Exception("ERROR: There is no wwtp with ID: " + str(ID))
    
def findInNetworkPath(archPathWWTP, sewers_Current, sewers):
    """
    Find the path already in a network. The function checks as well that there are no loops.
    In case the archPathWWTP creates loops, the last non looped node is used and the path
    calculated back to the nearest WWTP on the existing non-looped network.
    
    Input Arguments: 
    archPathWWTP          --    Path
    sewers_Current        --    Network before adding the nodes
    sewers                --    Network
    
    Output Arguments:
    lastPath              --    Path already in network
    """ 
    position, entryOld =  0, archPathWWTP[0][0]    # First entry in List
    nodesInNetwork = [entryOld]                    # List containing all nodes which are part of a network with largestWWTPID as closestWWTP
    
    # Iterate path until last element in list which is in sewers_Current is found
    for i in archPathWWTP[1:]:
        position += 1
        entry = i[0]
        if entry not in sewers_Current:
            break
        else:
            testDirection = sewers[entry][0]        # check if really flows on the path to origin
            if testDirection != entryOld:           # If the direction is not the last Node
                break
            else:
                nodesInNetwork.append(entry)
        entryOld = i[0]
    return nodesInNetwork

def getSummedFlow(nodes, ID):
    """
    This function sums the flow flowing into a node and the flow coming from the node itself.
    
    Input Arguments: 
    nodes             -    nodes
    ID                -    Id
    
    Output Arguments:
    summedFlow        -    Sum of flows
    """
    for i in nodes:  
        if i[0] == ID:
            summedFlow = i[4] + i[8]
            return summedFlow

def updateFlowA1(nodes, pathNearWTP, allPopNodesOntheWay):
    """
    This function updates the flow along a path to a wwtp.
    
    Input Arguments: 
    nodes                  --    nodes
    pathNearWTP            --    path
    allPopNodesOntheWay    --    all populated nodes on the path
    
    Output Arguments:
    allNodesToDelet        --    all nodes flowing to a wwtp
    """
    pathNearWTP, allPopNodesOntheWay = pathNearWTP[::-1], allPopNodesOntheWay[::-1]
    
    # Get all points in archPathList which are populated and already added to P
    for populatedNode in allPopNodesOntheWay[1:]: 
        noramlweiter, count = False, 0
        
        for entry in pathNearWTP:
            if entry[0] == populatedNode or noramlweiter == True:  # Go to populated node in path to closest WWTP
                noramlweiter = True
                
                # New additional flow is added along the tree
                for i in nodes:
                    if i[0] == entry[0]:
                        if count > 0:                
                            i[4] = FlowNewConnectedPoint + i[4] # add flow
                        if count == 0:                          # Sum Flow Menge
                            FlowNewConnectedPoint = i[8]        # Flow
                            count += 1
                        break
    return nodes

def getInflowingNodes(nodes, partOfNetwork, initialFlowWWTP, notInaNetwork):
    """
    The function iterates on the path between nodes and calculates the inflowing nodes from branches of this network to this path.
    
    Input Arguments: 
    sewersNoC            --    path before adding connection
    nodes             --    nodes
    partOfNetwork     --    path betwwen two wwtps
    initialFlowWWTP      --    source of initial wwtp
    notInaNetwork        --    Nodes not in a network
    
    Output Arguments:
    nodes_BI             --    Nodes with only inflowing flow.
    """
    # the way on which the nodes get only the inlet flow always needs to start at the WWTP
    if len(partOfNetwork) > 1:                      # only one edge
        if len(partOfNetwork) == 2:
            for i in nodes:
                if i[0] == partOfNetwork[1]:
                    flowUpSTrem = i[8] + i[4]
                    break
            for i in nodes:
                if i[0] == partOfNetwork[0]:
                    i[4] = i[4] - flowUpSTrem
                    break  
        else:                                       # more than one edge
            counter = 0
            for i in partOfNetwork: 
                if counter == 2:
                    for s in nodes:
                        if s[0] == i:
                            thirdID, thirdNodeFlow = s[0], round(s[4] + s[8], 7) # Round to 7 digits
                            break
                    
                    # change flow in node that only the inflow stays in i[4]
                    for f in nodes:
                        if f[0] == secondID:
                            currFlow = round(f[4], 7)               # Round to 7 digits
                            if f[4] != 0:                           # If current flow is zero, don't change. this means that the connection between two wwtps were unpouplated nodes
                                f[4] = currFlow - thirdNodeFlow     # If current flow is zero, don't change. this means that the connection between two wwtps were unpouplated nodes
                            break
                    
                    firstNodeFlow = secondNodeFlow
                    secondNodeFlow = thirdNodeFlow
                    secondID = thirdID
    
                if counter == 1:
                    counter = 2
                    for s in nodes:
                        if s[0] == i:
                            secondNodeFlow, secondID, flowForFirstElement = s[4], s[0], s[4] + s[8]
                            break
                if counter == 0:   
                    counter = 1

            # Assign initialFlow to wwtp        
            for i in nodes:
                if i[0] == partOfNetwork[0]:
                    if i[4] != 0:  
                        if partOfNetwork[1] in notInaNetwork:    # Check if against flow in path
                            i[4] = initialFlowWWTP  
                        else:
                            i[4] = i[4] - flowForFirstElement
                        break
        return nodes
    else:
        return nodes

def findNetworkToRemove(archPathWWTP, nodesToNetwork, nodesFromNetwork, sewers, WWTPs):
    """
    This function checks if there are intermediate nodes which are not in the beginning
    or starting network. If yes, search the wwtp of this node and add to list in order that
    all nodes flowing to these networks can be later removed.

    Input Arguments: 
    archPathWWTP               -    Path betwen two nodes
    nodesToNetwork             -    Nodes in network of ending node
    nodesFromNetwork           -    Nodes in network of starting node
    sewers                     -    Network before iteration
    WWTPs                      -    list of wwtps

    Output Arguments:
    wwtpsOnTheWayBetweenWTPs   -    wwtps on the way between two nodes.
    needsNetworkRemoving       -    Criteria whether a network needs to be removed or not (0: No, 1: Yes)
    """
    noNetworkNode, wwtpsOnTheWayBetweenWTPs, needsNetworkRemoving = [], [], 0
    fromNodeWWTP = archPathWWTP[0][0]               # Starting node which might be a wwtp
    toNodeWWTP = archPathWWTP[-1][1][0]             # Starting node which might be a wwtp
                                
    # Get a list with all nodes not in the from or to Network
    for i in archPathWWTP[1:]:
        if i[0] not in nodesFromNetwork and i[0] not in nodesToNetwork:
            noNetworkNode.append(i[0])

    # Check if in sewers and thus in a network. If in a network, go to closestWWTP and give out WWTP (in oder to remove whole network)
    for i in noNetworkNode:
        try:
            inNetworkNode = sewers[i][0]                                                    # Check to which wWTP the node flows
            if inNetworkNode == ():                                                         # If node itself is a wwtp
                if i not in wwtpsOnTheWayBetweenWTPs:
                    wwtpsOnTheWayBetweenWTPs.append(i)
            else:
                _, closestARA = getclosestWTPreconnection(inNetworkNode, WWTPs, sewers)
                if fromNodeWWTP != closestARA and toNodeWWTP != closestARA:                 # Check if closest ARA is a wwtp
                    if closestARA not in wwtpsOnTheWayBetweenWTPs:                          # If closest ARA is not the starting or ending node(which can be a wwpt)
                        wwtpsOnTheWayBetweenWTPs.append(closestARA)
        except:
            _ = 0
            
    if len(wwtpsOnTheWayBetweenWTPs) > 0:                                                   # Check if there needs to be a removal of points
        needsNetworkRemoving = 1
    else:
        needsNetworkRemoving = 0
    return wwtpsOnTheWayBetweenWTPs, needsNetworkRemoving

def updateFlowInWWTP(WWTPs, nodes, ID):
    """
    This function updates the flow in a wwtp.

    Input Arguments: 
    WWTPs       --    List of wwtps
    nodes       --    nodes
    ID          --    ID

    Output Arguments:
    WWTPs       --    Updated list of wwtps
    """
    for i in nodes:
        if i[0] == ID:
            newFlow = i[4] + i[8]
            break
    
    for i in WWTPs:
        if i[0] == ID:
            i[1] = newFlow
            break      
    return WWTPs

def delEntry(WWTPs, ID):        
    """
    This function deletes a wwtps with an ID.

    Input Arguments: 
    WWTPs            --    List of wwtps
    ID               --    ID of wwtp to delete

    Output Arguments:
    WWTPs            --    Updated list of wwtps
    """
    delPos = 0
    for i in WWTPs:
        if i[0] == ID:
            del WWTPs[delPos]
            break
        delPos += 1
    return WWTPs
         
def changeFlowAlongPath(nodes_BI, pathBetweenWWTPs):
    """
    Change flow along the path. Just sums the flow and previous flow. Starts at source.

    Input Arguments: 
    nodes_BI             --    nodes
    pathBetweenWWTPs     --    path between two wwtps

    Output Arguments:
    nodes_BI             --    Updated nodes
    """
    cnt = 0
    for i in pathBetweenWWTPs:
        if cnt == 1:
            for eintrag in nodes_BI:
                if eintrag[0] == i:
                    eintrag[4] = flowAbove + eintrag[4]  # add flowFromAbove plus inletFlow
                    flowAbove = eintrag[4] + eintrag[8]
                    break
        else:
            for s in nodes_BI:
                if s[0] == i:
                    flowAbove = s[4] + s[8]
                    cnt = 1
                    break         
    return nodes_BI

def addDEMPntstoNodes(nodes, archPathListMST, boundingCandidates, fromNode, toNode, minTD):
    """
    Add ponts on path to nodes. In case the node is a DEM-Points, get DEM Coordinates. .
    
    Input Arguments: 
    nodes                 --    nodes
    archPathListMST       --    path
    boundingCandidates    --    All dem points within a bounding box
    fromNode              --    From node
    toNode                --    To node
    minTD                 --    Minium Trench Depth

    Output Arguments:
    nodes                 --    Updated nodes
    """
    if len(archPathListMST) > 1:
        for edge in archPathListMST:
            append, newPoint1 = 1, edge[0]

            for i in nodes:  # Check wheter already in nodes
                if i[0] == newPoint1:
                    append = 0
                    break

            if append == 1:
                newPoint2, heightP1, heightP2 = edge[1][0], edge[2], edge[3]
                for cell in boundingCandidates:
                    if newPoint1 != fromNode and newPoint1 != toNode and cell[0] == newPoint1:
                        nodes.append([newPoint1, cell[1], cell[2], heightP1, 0, 0, 0, 0, 0, [], heightP1 - minTD])

            # Add last Element    
            append, newPoint2 = 1, edge[3]
            
            if newPoint2 == newPoint1:
                append = 0
            if append == 1:
                newPoint2, heightP1, heightP2 = edge[1][0], edge[2], edge[3]
                for cell in boundingCandidates:
                    if newPoint2 != fromNode and newPoint1 != toNode and cell[0] == newPoint2:
                        nodes.append([newPoint2, cell[1], cell[2], heightP2, 0, 0, 0, 0, 0, [], heightP1 - minTD])
    return nodes

def getClosestNetworkWWTP(nodes, WWTPs, sewers, pZero, checkBackConnectionID):
    """
    This Function calculates the closeset network node (merging heuristic)
    Then it iterates on the found network until a wwtp is reached

    Input Arguments: 
    nodes                     --    nodes
    WWTPs                     --    List with wwtps
    sewers                    --    Sewer Network nodes
    pZero                     --    from Node
    checkBackConnectionID     --    ID of WWTP

    Output Arguments:
    closestWWTPinNet  --    Closest wwpt in a network
    cordWWTPinNet     --    Coordinate of wwtp in network
    """
    closestNetworkNode,  closestDistance = False, 9999999999
    currentNetworkNodes = breathSearch(checkBackConnectionID, sewers)      # find all nodes in current network
    
    for i in sewers:
        if sewers[i][0] != ():                                  # Check that not wwtp
            if i not in currentNetworkNodes:                    # In order that not own network node is found
                networkNode, corNet, _, _, _, _ = getPns(i, nodes)
                distanz3d, _, _ = distanceCalc3d(pZero, corNet) 
                if distanz3d < closestDistance:                 # Get closest unweighted wwtp
                    closestNetworkNode = networkNode
                    closestDistance = distanz3d                 # Replace distance
    
    # Closeset Network node
    if closestNetworkNode == False:
        return None, None, None                                 # There is no closeset network because all other wwtps are don't have networks
    
    # Iterate along path to find wwtp
    pathInClosestNetwork, _ = getclosestWTP(closestNetworkNode, WWTPs, sewers)
    closestWWTPinNet, cordWWTPinNet, _, _, _, _ = getPns(closestNetworkNode, nodes)    
    return closestWWTPinNet, cordWWTPinNet, pathInClosestNetwork


def flowIfSameCheck(WWTPs, aggregatetPoints):
    '''
    This function checks wheter after the SNIP calculation there was lost any flow. This is a controlling function.
    In case the flow is no the same throw an error (except when a connection is enforced and the flow cannot be the same).
    
    Input:
    WWTPs                 -    Waste water treatment plants
    aggregatetPoints      -    All aggregated nodes 
    '''
    sumFlowWWTP, sumFlowAggregated = 0, 0
    
    for i in WWTPs:
        sumFlowWWTP +=i[1]
    
    for i in aggregatetPoints:
        sumFlowAggregated += i[8]

    sumFlowWWTP = round(sumFlowWWTP, 3)
    sumFlowAggregated = round(sumFlowAggregated, 3)
    
    if sumFlowAggregated != sumFlowWWTP:
        arcpy.AddMessage("sumFlowWWTP: " + str(sumFlowWWTP))
        arcpy.AddMessage("sumFlowAggregated: " + str(sumFlowAggregated))
        raise Exception("Error: The total flow before SNIP and after SNIP is not identical.")

def calcZValue(aggregatetPoints, forceZ, iterativeCostCalc, hypoWWTPcorrectFlow, WWTPs):
    '''
    This function calculates the Z value according to Ambros (2000).
    
    Input:
    aggregatetPoints         -    Aggregated Nodes
    forceZ                   -    Criteria wheter a hihger Z value is enforced
    iterativeCostCalc        -    Criteria wheter costs are calculated in each iteration
    hypoWWTPcorrectFlow      -    Flow in case a higher Z value is inforced
    WWTPs                    -    waste water treatment plants
    '''
    sources = len(aggregatetPoints)
    if forceZ == True and iterativeCostCalc == 1:
        sinks = len(hypoWWTPcorrectFlow)
    else:
        sinks = len(WWTPs)
    degCen = round((float(sources) - float(sinks)) / float(sources), 4)   
    return degCen

def calcZweighted(forceZ, iterativeCostCalc, hypoWWTPcorrectFlow, sewers, aggregatetPoints, WWTPs):
    '''
    This function calculates the Z value according to Eggimann et al. (2015).
    
    Input:
    forceZ                     -    Criteria wheter a connection is enforced
    iterativeCostCalc          -    List storing costs
    hypoWWTPcorrectFlow        -    Flow in case a connection is enforced
    sewers                     -    Sewer Network
    aggregatetPoints           -    Aggregated nodes
    WWTPs                      -    Waste water treatment plants
    
    '''
    if forceZ == True and iterativeCostCalc == 1:
        listWWTPwithAggregatedNodes = getAggregatedNodesinListWWTP(hypoWWTPcorrectFlow, sewers, aggregatetPoints) # Get aggregated nodes
    else:
        listWWTPwithAggregatedNodes = getAggregatedNodesinListWWTP(WWTPs, sewers, aggregatetPoints) # Get aggregated nodes
    
    sumFlow, weightedTerm = 0, 0
    for i in listWWTPwithAggregatedNodes:
        sumFlow += i[1]
        weightedTerm += float((i[1]/i[2]))
        
    degCenWeighted = (sumFlow - weightedTerm) / sumFlow
    return degCenWeighted

def turnNodesIntoWWTP(WWTPs, nodes, aggregatetPoints, sewers):
    '''
    This function searchs all not yet connected nodes and adds them to the WWTP. This is used in case SNIP is aborted before reaching the MM.
    
    Input:
    WWTPS             -    WWTPs
    nodes             -    Nodes
    aggregatetPoints  -    Aggreagted nodes
    sewers            -    Sewers
    
    Output:
    final_wwtps       -    All nodes added to WWTPs
    '''
    final_wwtps = []
    for i in WWTPs:
        for e in nodes:
            if i[0] == e[0]:
                final_wwtps.append([i[0], i[1], e[1], e[2], ])   # ID, flow, X, Y
                break

    # ADd not yet connected nodes
    for i in aggregatetPoints:
        if i[0] not in final_wwtps:
            isInSewer = 0
            try:
                _ = sewers[i[0]][0]
                isInSewer = 1
            except:
                isInSewer = 0
                    
            for z in nodes:                                                         # Get all infos from nodes
                if z[0] == i[0] and z[8] > 0 and z[4] == 0 and isInSewer == 0:      # Read out correct list element
                    final_wwtps.append([i[0], z[8], z[1], z[2]])                    # ID, flow, x, y
                    break    
    return final_wwtps

def connectivityPotential(nodes, WWTPs, idWWTP, pZero, f_merge): 
    """
    Find node with highest connectivity-potential in wwtps. Plus the function finds the closest (euclidian distance) wwtp

    Input Arguments: 
    nodes                      --    nodes
    WWTPs                      --    List with wwtps
    idWWTP                 --    To Node
    pZero                      --    From Node
    f_merge                    --    Puts the size in relation ot the distance. If large, the size gets more important

    Output Arguments:
    potentialNode              --    Potential node with Index
    potentialNodeCoordinates   --    Coordinates of potential node
    closestWWTP                --    Id of closest WWTP
    closestCoordinates         --    Coordinates of closest WWTP
    """
    potIndex, closestDistance, potentialNode, closestWWTP = 9999999999, 9999999999, None, None
    
    for f in WWTPs:
        if f[0] != idWWTP:    
            idWWTPONE, pOne, _, _, forceConnection, _ = getPns(f[0], nodes)

            if forceConnection == 0:
                distanceinclSlope, _, _ = distanceCalc3d(pZero, pOne)           # calculate 3d distance
                size = f[1]               # size is in m3 
                
                # Get closest unweighted wwtp
                if distanceinclSlope < closestDistance:
                    closestWWTP, closestCoordinates, closestDistance = idWWTPONE, pOne, distanceinclSlope
                connectivityIndex = distanceinclSlope * size**(-1 * f_merge)
                    
                # returns the node with the hightes potential
                if connectivityIndex < potIndex:
                    potentialNode, potentialNodeCoordinates, potIndex = idWWTPONE, pOne, connectivityIndex
    if potentialNode == None:                                           # No Connection was found
        if closestWWTP == None:                                         # No Connection was found
            return None, None, None, None
        else:
            return None, None, closestWWTP, closestCoordinates
    else:
        return potentialNode, potentialNodeCoordinates, closestWWTP, closestCoordinates
        
def addOriginTolistWWTPs(nodes, WWTPs, origin):
    """
    Add starting node to listwwtps 

    Input Arguments: 
    nodes           --    nodes
    WWTPs           --    List with wwtps
    origin          --    starting node
   
    Output Arguments:
    WWTPs           --    List of wwtps
    """
    flowCheck = 0
    for i in nodes:
        if i[0] == origin:
            WWTPs.append([origin, i[4] + i[8]])
            flowCheck = i[4] + i[8]
            break
    if flowCheck == 0:
        raise Exception("ERROR: INITIAL NODE HAS NO FLOW. Select Different Starting Node")
    return WWTPs

def updatePathLong(nodes, nodes_invert, path, sewersNoC, toNode):
    """
    This function updates the correct flow in the nodes on a path.

    Input Arguments: 
    nodes                      --    nodes
    nodes_invert               --    nodes inversed direction flow
    path                       --    Inverted path to nearest wwtp
    sewersNoC                  --    Network before added new path
    toNode                     --    Source
    
    Output Arguments:
    nodes_invert              --    Updated nodes
    """
    counter, finishReadingPath = 0, 0
    
    # Change Flow of first node
    for f in nodes_invert:
        if f[0] == path[0]:
            flowStartNode = f[4]                            # Flow before [4] is set to zero in starting node
            
            # Against Flow
            for i in nodes:
                if i[0] == path[1]:                         # next element in path
                    if i[4] != 0:
                        flowAgainstDirection = i[4] + i[8]
                        break
                    else: 
                        if i[0] in sewersNoC:
                            flowAgainstDirection = i[8]
                        else:
                            flowAgainstDirection = 0        # node was node connected, meaning that there is no againstFlow
                        break
                                    
            # Flows in Node
            for t in nodes:
                if t[0] == path[0]:
                    flowInNode = t[4]
                    break

            startFlow = flowInNode - flowAgainstDirection
            f[4] = startFlow 
            break
                        
    # Iterate path to swap wwtp
    for entry in path: 
        counter += 1  # get next element
        if counter == 1:
            newest, secondNewest = 0, 0
        if counter == 2:
            secondNewest = 0  
        thirdNewest = secondNewest
        secondNewest = newest

        for i in nodes_invert: 
            if i[0] == entry:
                newest = i  
                break
                      
        # Needed to get three points
        if counter > 2 and finishReadingPath == 0:
            newest, second, third = newest, secondNewest, thirdNewest

            # Change flow along path to wwtp
            for i in nodes_invert:
                if i[0] == second[0]:                                                                                       # Change flow in secondNewest
                    if newest[0] == toNode:                                                                                 # Last 3 nodes
                        if secondNewest[8] == 0:                                                                            # Arrived at the end of path or only three entries
                            if secondNewest[4] == 0:                                                                        # Connection to archPoint. Take over flow from node below
                                updatedFlowSecondLast = third[4] + third[8]
                                i[4] = updatedFlowSecondLast                                                                # change second last node                                  
                                nodes_invert = changeFlowInNode(nodes_invert, toNode, updatedFlowSecondLast)                # change last node  
                                break 
                            else:  
                                # only three nodes 
                                if third[0] == path[0]:
                                    flowBefore = third[4] + third[8]                                                        # As wwtp has now flow
                                    if flowStartNode == second[4]:                                                          # Check if there is inflow from other nodes
                                        if newest[0] not in sewersNoC:
                                            notAgainstFlow = second[4]
                                        else:
                                            notAgainstFlow = second[4] - (newest[4] + newest[8])
                                    else:
                                        if newest[0] not in sewersNoC:
                                            notAgainstFlow = second[4]
                                        else:
                                            notAgainstFlow = second[4] - (newest[4] + newest[8])
                                else:
                                    flowBefore = third[4] + third[8]
                                    if newest[0] not in sewersNoC:
                                        notAgainstFlow = second[4]
                                    else:
                                        notAgainstFlow = second[4] - (newest[4] + newest[8])
                                updatedFlowSecondLast = notAgainstFlow + flowBefore
                                i[4] = updatedFlowSecondLast                                                                # change second last node
                                nodes_invert = changeFlowInNode(nodes_invert, toNode, updatedFlowSecondLast)                # change last node
                                break  
                        else:  
                            if secondNewest[4] == 0:                                                                        # no flow
                                flowBefore = third[4] + third[8]
                                updatedFlowSecondLast = flowBefore
                                flowForLast = second[8]                            
                                i[4] = updatedFlowSecondLast                                                                # change second last node
                                nodes_invert = changeFlowInNode(nodes_invert, toNode, updatedFlowSecondLast + flowForLast)  # change last node
                                break           
                            else: 
                                # first three nodes 
                                if third[0] == path[0]:
                                    flowBefore = third[4] + third[8]
                                            
                                    # Check if there is inflow from other nodes
                                    if flowStartNode == second[4] + second[8]:
                                        if newest[0] not in sewersNoC:
                                            notAgainstFlow = second[4]
                                        else:
                                            notAgainstFlow = second[4] - (newest[4] + newest[8])
                                    else:
                                        if newest[0] not in sewersNoC:  
                                            notAgainstFlow = second[4]  
                                        else:  
                                            notAgainstFlow = second[4] - (newest[4] + newest[8])
                                else:
                                    flowBefore = third[4] + third[8]
                                            
                                    if newest[0] not in sewersNoC:      # If neweswt not in sewers, then the new branch is reached. Take over all flow 
                                        notAgainstFlow = second[4]          
                                    else:
                                        notAgainstFlow = second[4] - (newest[4] + newest[8])       
                                flowForLast = secondNewest[8]
                                updatedFlowSecondLast = notAgainstFlow + flowBefore
                                i[4] = updatedFlowSecondLast                                                                # change second last node
                                nodes_invert = changeFlowInNode(nodes_invert, toNode, updatedFlowSecondLast + flowForLast)  # change last node
                                break 
                        finishReadingPath = 1
                        break
                    else:   # Not Last 3 nodes, More than three nodes
                        if secondNewest[8] == 0:
                            if secondNewest[4] == 0:                                                         # Connection to archPoint (unconnected point)
                                updatedFlowSecondLast = third[4] + third[8]                                  # change second last
                                i[4] = updatedFlowSecondLast                                                 # change second last node
                            else:  
                                if third[0] == path[0]:
                                    flowBefore = third[8]
                                    if flowStartNode == second[4]:                                           # Check if there is added flow
                                        if newest[0] not in sewersNoC:                                       # If neweswt not in sewers, then the new branch is reached. Take over all flow 
                                            notAgainstFlow = second[4]
                                            flowBefore = third[4] + third[8]  
                                        else:
                                            notAgainstFlow = second[4] - (newest[4] + newest[8])
                                    else: 
                                        if newest[0] not in sewersNoC:
                                            notAgainstFlow = second[4]
                                            flowBefore = third[4] + third[8] 
                                        else:  
                                            notAgainstFlow = second[4] - (newest[4] + newest[8])
                                            flowBefore = third[4] + third[8]  
                                else:                                                                        # not first 3 Points 
                                    flowBefore = third[4] + third[8]
                                    if newest[0] not in sewersNoC:                                           # If neweswt not in sewers, then the new branch is reached. Take over all flow 
                                        notAgainstFlow = second[4]
                                    else:
                                        notAgainstFlow = second[4] - (newest[4] + newest[8])   
                                updatedFlowSecondLast = notAgainstFlow + flowBefore 
                                i[4] = updatedFlowSecondLast                                                 # change second last node
                        else:                                                                                # is inhabited because has flow
                            if secondNewest[4] == 0:                                                         # inhabited point, no flow
                                if newest[0] not in sewersNoC:
                                    updatedFlowSecondLast = third[4] + third[8]                              # change second last
                                    i[4] = updatedFlowSecondLast                                             # change second last node
                                else:
                                    updatedFlowSecondLast = (third[4] - (newest[4] + newest[8])) + third[8]  # change second last
                                    i[4] = updatedFlowSecondLast                                             # change second last node 
                            else:                                                                            # There is flow
                                if third[0] == path[0]:                                                      # First 3 Points
                                    flowBefore = third[4] + third[8]
                                    if flowStartNode == second[4]:                                           # Check if there is added flow
                                        if newest[0] not in sewersNoC:
                                            notAgainstFlow = second[4]
                                        else:
                                            notAgainstFlow = second[4] - (newest[4] + newest[8])
                                    else:
                                        if newest[0] not in sewersNoC:
                                            notAgainstFlow = second[4]
                                        else:
                                            notAgainstFlow = second[4] - (newest[4] + newest[8])     
                                else:                                                                        # not first 3 Points
                                    flowBefore = third[4] + third[8]                                         # change second last node   
                                    if newest[0] not in sewersNoC:
                                        notAgainstFlow = second[4]
                                    else:
                                        notAgainstFlow = second[4] - (newest[4] + newest[8])
                                updatedFlowSecondLast = notAgainstFlow + flowBefore
                                i[4] = updatedFlowSecondLast                                                 # change second last node
                    break
    return nodes_invert

def updatePathShort(nodes_Invert, pathtonearestWTPInvert):
    """
    This function updateds the correct flow in the nodes on the path between two nodes.

    Input Arguments: 
    nodes_Invert              -    nodes inversed direction flow
    pathtonearestWTPInvert    -    Inverted path to nearest wwtp
    
    Output Arguments:
    nodes_Invert              -    Updated nodes
    """
    for s in nodes_Invert:
        if s[0] == pathtonearestWTPInvert[0]:
            flowFromNodeSwap = s[4] + s[8]
            break
    for d in nodes_Invert:
        if d[0] == pathtonearestWTPInvert[1]:
            d[4] = d[4] + flowFromNodeSwap
            break
    return nodes_Invert

def selectLowestNode(aggregatetPoints):
    """
    This function selects the lowest node of all aggregated nodes.

    Input Arguments: 
    aggregatetPoints     --    aggregated sources
    
    Output Arguments:
    lowestNode           --    Id of lowest node
    """
    height = 9999999999
    for i in aggregatetPoints:
        if i[8] < height:
            lowestNode, height = i[0], i[3]
    return lowestNode

def invertFlowToNearestWTP(pathtonearestWTPswap, sewers_A3, sewers):
    """
    This function selects the lowest node.

    Input Arguments: 
    pathtonearestWTPswap --    aggregated sources
    sewers_A3            --    Sewer network situation 3
    sewers               --    Sewer network 
    
    Output Arguments:
    sewers_A3            --    Updated Network
    """
    skip = 0
    for e in pathtonearestWTPswap:  
        if skip > 0:
            if sewers[e[0]] != ():  
                distEdgen = sewers_A3[e[0]][1]
                sewers_A3[oldEle] = (e[0], distEdgen)
        if skip == 0 or skip == 1:
            oldEle, skip = e[0], 1
    sewers_A3[pathtonearestWTPswap[-1][0]] = (), 0  # new origin
    return sewers_A3

def sumCostOfArchWWTPs(allPopNodesOntheWay, nodes, EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement):
    """
    This function sums the costs of all wwtp on the archpath. 

    Input Arguments: 
    allPopNodesOntheWay     --    All popluated nodes on the way
    nodes                   --    Nodes
    EW_Q                    --    Population equivalent in liters
    wwtpLifespan            --    Discounting years
    interestRate            --    Real interest rate
    
    Output Arguments:
    summedCostsWWPTS        --    total costs
    """
    summedCostsWWPTS = 0
    
    for decentralizedWTP in allPopNodesOntheWay:  # Last one is not considered as it is calculated just above
        for i in nodes:
            if i[0] == decentralizedWTP:
                flowWTPonTheWay = costWWTP((i[8] + i[4]), EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement)  # frher False: flowToNOde
                summedCostsWWPTS = summedCostsWWPTS + flowWTPonTheWay
                break
    return summedCostsWWPTS       

def getPathToWTPSubentwork(pathtonearestWTP, sewers):
    """
    This function gets the path to tbe closeset network

    Input Arguments: 
    pathtonearestWTP              --    Path to nearest WWTP
    sewers                        --    Network

    Output Arguments:
    pathSubNetworkToClosestWWTP   --    Path to closest wwtp
    """
    pathSubNetworkToClosestWWTP = []               
    for i in pathtonearestWTP:
        if i in sewers:
            pathSubNetworkToClosestWWTP.append(i) 
        else:
            continue
    return pathSubNetworkToClosestWWTP

def initialSetNetwork(sewers, toNode, fromNode, minDit): 
    """
    This function sets the initial configurations for the SNIP Algorithm

    Input Arguments: 
    sewers            --    Network
    toNode            --    Node 
    fromNode          --    Node
    minDit            --    Distance between the nodes

    Output Arguments:
    sewers           --    Initial network for SNIP
    """
    sewers[toNode] = fromNode, minDit
    return sewers  

def costComparison(totCostI, totCostII, totCostIII, resonableCosts, costConnection, ignoreRC, forceCentralConnection, finishedMerging):
    """
    This function compares different cost options by considering reasonable costs. 
    If reasonable costs are considered, a central connection is choosen even if
    decentral option would be cheaper (but central connection needs to be cheaper than reasonable costs).
    Reasonsable costs only matter if decentral conneciton is cheaper
    than central connection and central connection is cheaper than resonable costs. 
    
    Input Arguments: 
    totCostI, totCostII, totCostA3     --   Different costs for different systems
    resonableCosts                     --   reasonable costs
    costConnection                     --   Total costs of connection to check 
    forceCentralConnection             --   If a connection is forced, take cheapes central connection
    finishedMerging                    --   Only Force connection if Expansion Module is finished
    
    Output Arguments:
    wwtpSWAP                           --   Criteria wheter inversed flow is choosen.
    decCrit                            --   Criteria wheter a decentral option is choosen. 
    """

    # If a central connection is forced:
    if forceCentralConnection == True and finishedMerging == 1:
        if totCostI < totCostIII:
            decCrit, wwtpSWAP = 1, 0        # Selection central Option I
        else:
            wwtpSWAP, decCrit = 1, 1        # Selection central Option III
    else: 
        # if central connection costs are lower than reasonable costs, force a central connection.
        if ignoreRC != 1:
            if costConnection < resonableCosts: 
                # Connection Costs are lower than reasonable costs. Central connection is forced.
                if resonableCosts != 0:  
                    if totCostI < totCostIII:
                        wwtpSWAP, decCrit = 0, 1
                    else:
                        wwtpSWAP, decCrit = 1, 1
            else:
                if totCostIII >= totCostII: 
                    # Connection costs are higher than reasonable costs, decentral alternative can be considered
                    if totCostI < totCostII:                    # I is cheapest
                        decCrit, wwtpSWAP = 1, 0
                    else:
                        if totCostIII == totCostII:             # dez is same as swap --> go for swap
                            wwtpSWAP, decCrit = 1, 1
                        else:
                            if totCostI == totCostII:           # all options are the same --> Go for simple connect
                                decCrit, wwtpSWAP = 1, 0
                            else:                               # II is cheapest
                                wwtpSWAP, decCrit = 0, 0
                else:                                           # cost III are smaller than II
                    if totCostI < totCostIII:                   # I is cheapest
                        decCrit, wwtpSWAP = 1, 0
                    else:
                        if totCostIII == totCostI:              # regular is identical to swap --> go for regular
                            wwtpSWAP, decCrit = 0, 1
                        else:                                   # III is cheapest
                            wwtpSWAP, decCrit = 1, 1               
            return wwtpSWAP, decCrit
        else:                                                   # Don't consider reasonable costs
            if totCostIII >= totCostII: 
                if totCostI < totCostII:                        # I is cheapest
                    decCrit, wwtpSWAP = 1, 0
                else:
                    if totCostIII == totCostII:                 # dez is same as swap --> go for swap
                        wwtpSWAP, decCrit = 1, 1
                    else:
                        if totCostI == totCostII:               # all options are the same --> Go for simple connect
                            decCrit, wwtpSWAP = 1, 0
                        else:                                   # II is cheapest
                            wwtpSWAP, decCrit = 0, 0
            else:                                               # cost III are smaller than II
                if totCostI < totCostIII:                       # I is cheapest
                    decCrit, wwtpSWAP = 1, 0
                else:
                    if totCostIII == totCostI:                  # regular is identical to swap --> go for regular
                        wwtpSWAP, decCrit = 0, 1
                    else:                                       # III is cheapest
                        wwtpSWAP, decCrit = 1, 1
    return wwtpSWAP, decCrit

def delWWTP(WWTPs, ID):
    """
    This function deletes a wwtp.
    
    Input Arguments: 
    WWTPs         --   List with the wwtps
    ID            --   Id of WWTP to delete

    Output Arguments:
    WWTPs         --   Updated list of wwtps.
    """
    cnt = 0
    for i in WWTPs:
        if i[0] == ID:
            del WWTPs[cnt]
            return WWTPs
        cnt += 1
        
def calcNodesWithFlow(nodes):
    """
    This function counts the number of nodes which are a source.
    
    Input Arguments: 
    nodes                 --    nodes

    Output Arguments:
    anzBewNodes           --    Nr of sources
    """
    anzBewNodes = 0  # As the starting node does not need to be connected
    for i in nodes:
        if i[8] > 0:
            anzBewNodes += 1    
    return anzBewNodes
          
def loopTest(pathNearWTP, sewersNoC):
    """
    Test if in pathNearWTP are nodes which would create loops.
    
    Input Arguments: 
    pathNearWTP       --    path
    sewersNoC         --    Network before adding path.

    Output Arguments:
    pathNearWTP       --    Nr of sources
    loops             --    Criteria wheter loop was found
    """
    # Iterate path until last element is in existing network. From this node, calculate back path to closest wwtp
    position = -1
    for i in pathNearWTP:       
        if i[0] in sewersNoC:                   # if i[0] is not in path, Don't do anything and look if after leaving the network another node in the nwtorks is reached
            positionWhereInNetwork = position   # Last node in network
        position += 1
    
    positionWhereInNetwork += 1
    pathInNetwork = pathNearWTP[positionWhereInNetwork:]
    nodeInNet = pathInNetwork[0][0]
    startFindingPath = nodeInNet
    
    # Path in Network of sewersNoC (path with no loops to wwpt)
    pathWithoutLoops = []
    while 1:
        entry = sewersNoC[startFindingPath]
        if entry[0] == (): 
            pathWithoutLoops.append([startFindingPath, entry])
            break
        else:
            pathWithoutLoops.append([startFindingPath, entry])
            startFindingPath = entry[0]
    
    # Check if the path from this poin in sewersNoC is identical to pathNearWTP (which was calculated with added ArchPath)
    pathWithoutLoops = InvertandswapID(pathWithoutLoops)
    new = []
    
    for i in pathWithoutLoops:
        new.append(i[1])
    
    add = False
    for restNode in pathNearWTP:
        if add == True:
            new.append(restNode)
        if restNode[0] == nodeInNet:
            add = True
    it = range(len(new))
    
    # Compare paths and see if the path was changed
    for i in it:
        if new[i][0] == pathNearWTP[i][0]:
            continue
        else:  
            loops = True                    # Info: ArchPath has loops
            return new, loops
    loops = False
    return pathNearWTP, loops

def changeTD(nodes, edgesList, pumps, pathToNetwork, maxTD, minSlope, inflowNodes, sewers, minTD):
    """
    This function changes the trench depth. The function starts at the end of a path and then 
    and changes the trench depth in order that the maximum trench depth and minimum slope criteria 
    are fulfilled. In case there are inflowing branches, the trench depth is not changed. In case 
    the maximum trench depth is not fulfilled, a pump is added to the pump list.
    
    Input Arguments: 
    nodes             --    List with nodes
    edgesList         --    List with edges
    pumps             --    List with pumps
    pathToNetwork     --    Path 
    maxTD             --    Max Trench depth
    minSlope          --    Minimum slope criteria
    inflowNodes       --    Branch nodes
    sewers            --    Sewer Network
    minTD             --    Minimum Trench Depth
    
    Output Arguments:
    nodesCopy         --    Updated nodes
    pumps             --    List with updated pumps
    edgesList         --    List with edges where slope was recalulated is slope was laid (slope of pipes)
    """                                               
    nodesCopy = fastCopyNodes(nodes)                                                
    edgesIDCopy = fastCopy(edgesList)
    pathToNetwork = InvertandswapID(pathToNetwork)                               # inverse and swap path
    
    for i in pathToNetwork:                                                      # iterate over path and change trechdepth if needed up to maximum trench depth
        fromID, toID = i[0], i[1][0] 
        length = i[1][1]                                                         # Read length
        _, trenchFrom, hFrom = getTrenchDepth(nodesCopy, fromID)                 # Get trenchHeight 
        positionToNode, trenchTo, hTo = getTrenchDepth(nodesCopy, toID)          # Get trenchHeight 
        minRequiredHDiff = (minSlope * length) / 100                             # Minimum required hight difference for free flow

        if trenchFrom - minRequiredHDiff > hTo - maxTD:                          # Check if slope is steep enough and possible trench depth not too deep
            pumpFlow, lastID = False, pathToNetwork[-1][1][0]
            nodesCopy, pumps = checkTrenchLifting(lastID, sewers, pumps, nodesCopy, nodes, trenchFrom, trenchTo, minRequiredHDiff, positionToNode, inflowNodes, toID, fromID, minSlope, pumpFlow, hTo, hFrom, minTD, maxTD)
        else:
            pumpFlow, lastID = True, pathToNetwork[-1][1][0]
            nodesCopy, pumps = checkTrenchLifting(lastID, sewers, pumps, nodesCopy, nodes, trenchFrom, trenchTo, minRequiredHDiff, positionToNode, inflowNodes, toID, fromID, minSlope, pumpFlow, hTo, hFrom, minTD, maxTD)

        IDnew, cord_new, _, _, _, trenchTo = getPns(toID, nodesCopy)
        IDold, cord_old, _, _, _, trenchFrom = getPns(fromID, nodesCopy)
        newPipeSlope =  (trenchFrom - trenchTo) / length                        # New pipe slope
        
        # Update slope in edgesIDCopy which becomse the slope of the pipe
        for e in edgesIDCopy:
            if e[0][0] == fromID and e[1][0] == toID:
                e[3] = newPipeSlope * - 1                                       # Change slope with pipe installation
                break
                
            if e[1][0] == fromID and e[0][0] == toID:
                e[3] = newPipeSlope                                             # Change slope with pipe installation
                break
        
        edgesIDCopy = addToEdgeList(edgesIDCopy, length, newPipeSlope, IDnew, IDold, cord_new, cord_old)
    return nodesCopy, pumps, edgesIDCopy                     

def correctCoordinatesAfterClip(aggregatetPoints, streetVertices):
    '''
    Because after ArcGIS Clipping slightly different coordinate endings might occur, the coordinates
    in the ggregatetPoints file needs to be adapted.
    
    Input:
    aggregatetPoints    -    Aggregated nodes
    streetVertices      -    Nodes of street network
    
    Output:
    aggregatetPoints    -    Updated aggregated nodes
    '''
    # Read out new coordinates of clipped street
    copy_aggregatetPoints = []
    compareDigits = 1 # If same coordinate to two significatn digits, consither them as the same
    for i in aggregatetPoints:
        oldX, oldY = i[1], i[2] # Before splitting
        foundNewCoordinate = 0

        for e in streetVertices:
            if round(oldX, compareDigits) == round(e[1], compareDigits) and round(oldY, compareDigits) == round(e[2], compareDigits):
                newX = e[1]
                newY = e[2]
                foundNewCoordinate = 1
                break

        if foundNewCoordinate == 1:
            copy_aggregatetPoints.append([i[0], newX, newY, i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10]])
        else:
            copy_aggregatetPoints.append(i)
        
    aggregatetPoints = copy_aggregatetPoints
    return aggregatetPoints

def getTrenchDepth(liste, ID):
    """
    This function gets the trench depth of a node.       
    
    Input Arguments: 
    liste                --    Nodes
    ID                   --    ID of which to get trench depth.

    Output Arguments:
    toID                 --    Updated nodes
    positionToNode       --    Position in list
    trenchTo             --    Trench height
    height               --    Actual terrain height
    """
    positionToNode = 0
    for d in liste:
        if d[0] == ID:
            height, trenchTo = d[3], d[10]
            break  
        positionToNode += 1
    return positionToNode, trenchTo, height

def writeInEdges(pntsOnPath):
    """
    This function converts a path into edges.    
    
    Input Arguments: 
    pntsOnPath    --    path

    Output Arguments:
    edges         --    edges
    """
    edges, pos = [], 0

    for i in pntsOnPath:
        pntTo = i[0]
        if pos == 1:
            edges.append([pntFrom[0], [pntTo, pntFrom[1]]])
            pntFrom = i
        else:
            pntFrom = i
            pos = 1
    return edges

def costPump(pathBetweenWWTPs, pumps, pricekWh, pumpYears, interestRate):  # Calculate pump costs
    """
    This function calculates the pumping costs.
    
    Input Arguments: 
    pathBetweenWWTPs          --    path
    pumps                     --    list with pumpes
    pricekWh                  --    price per kWh
    pumpYears                 --    Lifespan of pumps
    interestRate              --    Real interest rate
    
    Output Arguments:
    pumpCosts                 --    Annual Pumping costs
    totalcostOverWholePeriod  --    Costs for the whole lifespan
    """
    pumpCosts, summingPumpCosts, totalCostOverWholePeriod = 0, 0, 0
    
    for i in pathBetweenWWTPs:
        for pmp in pumps:
            if pmp[0] == i:
                flow, heightDifference = pmp[1], pmp[2]
                summingPumpCosts, totalCostOverWholePeriod = getPumpCostsDependingOnFlow(flow, heightDifference, pricekWh, pumpYears, interestRate)   # pump is found on path
                pumpCosts += summingPumpCosts                                                                                                                                               # sum pumping costs
                break
    return pumpCosts, totalCostOverWholePeriod

def correctTD(nodes, path, minTD, listwtps, maxTD, minSlope, sewers):
    """
    This corrects trench depth of starting node and returns a list with all nodes which have inflowing edges
    
    Input Arguments: 
    nodes        --    nodes
    path         --    path
    minTD        --    Minium Trench Depth
    listwtps     --    List with wwtps
    maxTD        --    Max Trench Depth
    minSlope     --    Minium Slope
    sewers       --    Network
    
    Output Arguments:
    pnts         --    Update nodes
    inflowNodes  --    All inflowing nodes

    """
    inflowNodes = []
    pnts = fastCopyNodes(nodes)
    
    if len(path) == 0:
        return nodes, inflowNodes # Short path
    
    if len(path) == 1:
        pathlength = None
    if len(path) > 1:
        pathlength = True
    
    # Change trench depth in starting node
    for i in pnts:
        if i[0] == path[0]:
            flowupStream = i[4] + i[8]                                                              # initial flowupStream
            iterTrench = 9999999                                                                    # used for finding min possible trench depth
            trenchChange = 0                                                                        # Check if not a wwtp, the correct initial trench depth is changed
               
            # Check if the initial node is a wwpt with inflowing nodes.  If it is a wwtp with inflow from other edges, change the trenchdepth to minium possible depth
            if pathlength == True:                                                                  # If more than one node
                for toWTPflowwingNode in sewers:                                                    # Get all nodes flowing to this wwtp:
                    if sewers[toWTPflowwingNode][0] == path[0] and toWTPflowwingNode != path[1]:
                        inflowNodes.append(path[0]) 
                                
                        for s in nodes:
                            if s[0] == toWTPflowwingNode:
                                hFlowingToNode = s[10]
                                break
                        for s in nodes:
                            if s[0] == path[0]:
                                hNode, krit = s[10], i[3] - i[10]
                                break
                                
                        # Calculate minimal trench depth of inflowing nodes
                        dist = sewers[toWTPflowwingNode][1]
                        depth = (float(minSlope) / float(100.0)) * dist
                        neededTrenchDepth = round((hFlowingToNode - depth), 2)                       # Height of node which flows there    

                        if krit != minTD and hFlowingToNode > neededTrenchDepth and neededTrenchDepth < (hFlowingToNode - minTD):  # Check if toflowing node needs to be pumped. If so, don't changed heigth! (trenchdepth of toflowing node is deeper than trenchheightto
                            for n in nodes:
                                if n[0] == toWTPflowwingNode:
                                    if neededTrenchDepth < iterTrench:
                                        trenchChange = 1
                                        if neededTrenchDepth > hNode - minTD:
                                            newTD = i[3] - minTD
                                            iterTrench = newTD
                                        else:
                                            newTD = neededTrenchDepth
                                            iterTrench = newTD
                                    break   
            else:
                # Get all nodes flowing to this wwtp:
                for toWTPflowwingNode in sewers:
                    if sewers[toWTPflowwingNode][0] == path[0]:
                        inflowNodes.append(path[0])  
                        for s in nodes:
                            if s[0] == toWTPflowwingNode:
                                hFlowingToNode = s[10]
                                break
                        
                        for s in nodes:
                            if s[0] == path[0]:
                                hNode = s[10]
                                krit = i[3] - i[10]
                                break
                        
                        dist = sewers[toWTPflowwingNode][1]
                        depth = (float(minSlope) / float(100)) * dist
                        neededTrenchDepth = round((hFlowingToNode - depth), 2)  # Height of node which flows there    
                        
                        if krit != minTD and  hFlowingToNode > neededTrenchDepth and neededTrenchDepth < (hFlowingToNode - minTD):
                            for n in nodes: 
                                if n[0] == toWTPflowwingNode:
                                    if neededTrenchDepth < iterTrench:  # Calculate minimal trench depth of inflowing nodes
                                        trenchChange = 1
                                        if neededTrenchDepth > hNode - minTD:
                                            newTD = i[3] - minTD
                                            iterTrench = newTD
                                        else:
                                            newTD = neededTrenchDepth
                                            iterTrench = newTD
                                    break 
            if trenchChange == 1:
                i[10] = newTD                   # Don't change depth
            else:
                i[10] = i[3] - minTD            # Change depth
            break
        
    # Iterate path. If the flow chang es along the path, there is inflow. Save the node in inflow list. If a node with no inflow, correct the trench depth to minimum Trench depth.    
    for entry in path[1:]:
        for i in pnts:
            if i[0] == entry:
                if flowupStream != i[4]:        # If flow in current node is not equal to flow in node above
                    inflowNodes.append(entry)   # Add node to a list with all inflowing nodes. Don't change trench depth
                else:
                    i[10] = i[3] - minTD        # Change trench depth to minium trench depth
                flowupStream = i[4] + i[8]      # Flow upstream in order to compare flow in node down the stream
                break              
    return pnts, inflowNodes

def readPumpList(pumps, nodes):
    """
    This functions updates the list with the pumps.
    
    Input Arguments: 
    pumps            --    List with pumps
    nodes            --    Nodes
    
    Output Arguments:
    pumpsReadOut     --    Updated pump List
    """
    pumpsReadOut = []
    for i in pumps:
        for d in nodes:
            if d[0] == i[0]:
                z = [i[0], d[1], d[2], i[1]]  
                pumpsReadOut.append(z)
                break
    return pumpsReadOut

def mergeLists(part_A, part_B):
    """
    This functions merges to lists.
    
    Input Arguments: 
    part_A, part_B    --    Lists
    
    Output Arguments:
    mergedList        --    Merged Lists
    """
    mergedList = []
    for i in part_A:
        mergedList.append(i)
    for i in part_B:
        mergedList.append(i)
    return mergedList

def replaceLoopPath(sewers, loopFound, loopNode, part_A):
    """
    This functions replaces path in order that there are no loops.
    
    Input Arguments: 
    sewers            --    Sewer Network
    loopFound         --    Criteria if loop was found
    loopNode          --    Node where loop was found
    part_A            --    Parth of the path up to the looping node
    
    Output Arguments:
    part_A            --    Corrected path
    """
    initialloopNode = loopNode
    if loopFound == 1:
        srapList = []           # List with path from loopNode to wwtp
        
        # Iterate sewers up to fromNode
        while 1:
            wayNode = sewers[loopNode]
            if wayNode[0] == ():
                break
            srapList.append([loopNode, wayNode])
            loopNode = wayNode[0]

        srapList = InvertandswapID(srapList)  # Swapt ID Invert
        
        # Append points after loopnode
        copykrit = 0
        for i in part_A:
            if copykrit == 1:
                srapList.append(i)
            if i[0] == initialloopNode:
                srapList.append(i)  
                copykrit = 1
        part_A = srapList
    return part_A

def replaceCorretPath(sewers, path, loopNode, firstNodeInToNetwork):
    """
    This functions replaces path in order that there is no looped path.
    
    Input Arguments: 
    sewers                  --    Network
    path                    --    Path
    loopNode                --    Node where loop was found
    firstNodeInToNetwork    --    First node which is not in network anymore
    
    Output Arguments:
    part_B                  --    Correct Path without loops
    """
    part_B, cnt = [], 0

    if loopNode == path[-1][0]:  # If no path between the two points, 
        part_B.append(path[-1])
        return part_B
    
    while firstNodeInToNetwork[0] != ():
        try:
            newNode = sewers[firstNodeInToNetwork[0]]
            if newNode[0] == ():
                break
            if cnt == 1:
                part_B.append((firstNodeInToNetwork[0], newNode))
            cnt = 1
            firstNodeInToNetwork = newNode
        except KeyError:
            newNode = sewers[firstNodeInToNetwork[1][0]]  # End is reached
            break
    return part_B
                           
def checkTrenchLifting(lastID, sewers, pumps, nodesCopy, nodes, trenchFrom, trenchTo, minhDiffRequired, posToNode, inflowNodes, toID, fromID, minSlope, pumpFlow, hTo, hFrom, minTD, maxTD):
    """
    Check if trench depth can be changed or is needed to be low because of inflowwing nodes (Check if toFlow is InflowNode).
    
    Input Arguments: 
    pathToNetwork         --    Path
    sewers                --    Network
    pumps                 --    List of pumps
    nodesCopy             --    copy of nodes
    nodes                 --    nodes
    trenchFrom            --    Trench Depth From node
    minhDiffRequired      --    Minium required height difference
    posToNode             --    Position To Node
    inflowNodes           --    Inflowing nodes
    toID                  --    ID To Node
    fromID                --    ID From Node
    minSlope              --    Slope Criteria
    pumpFlow              --    Criteria wheter the flow is pumped or not
    hTo                   --    Height To node
    hFrom                 --    Height from node
    minTD                 --    Minimum trench depth
    maxTD                 --    Maximum trench depth
    
    Output Arguments:
    nodesCopy             --    Updated nodes
    pumps                 --    Updated list of pumps
    """
    pumps = removePumpCheck(pumps, toID)                    # Check if at toflowing node a pump was installed. If yes, delete the pump
    Zmin = hTo - minTD                                      # Intermediate Calculation
    ZToDepth = trenchFrom - minhDiffRequired                # Intermediate Calculation
        
    # If node has inflowing nodes or last node is reached
    if toID in inflowNodes or toID == lastID:               # Node has inflowing nodes
        isFreistehendeWWTP, liftTrenchDepth = True, 1       # If no connecting edge is found, the wwtp is freistehend and thus trenchlift is always possible.  # 0 means that trench can't be lifted
   
        # Get nodes flowing to inflow node
        for e in sewers:  
            if sewers[e][0] == toID:
                isPump = checkIfIsPump(pumps, e)                                                        # Check if not a pump. If water is alread 
                if isPump == False:                                                                     # Only calculate trench depth for inflowing edges if not pumped   
                    isFreistehendeWWTP = False                                                          # criteria if is a detached wwtp
                    lengthInflowNode = sewers[e][1]                                                     # length to node in intework
                    _, trenchInflowFrom, _ = getTrenchDepth(nodes, e)                                   # Get trench depth in inflowing nodes  
                    _, trenchTOID, _ = getTrenchDepth(nodes, toID)                                      # Get trench depth in inflowing nodes  
                    hDiffInflowNod = trenchInflowFrom - trenchTOID                                      # height difference new
                    newSlopToInflow = round(float(hDiffInflowNod) / float(lengthInflowNode) * 100, 3)   # Calc slope. If positive, flows downstream
                    if newSlopToInflow <= minSlope:                                                     # if new calculated slope is less steep, lift trench. (and not goes upwards) RIESENBAUSTELLE
                        liftTrenchDepth = 0                                                             # inflow node can be changed and is no problem
                break
            
        if isFreistehendeWWTP == True:                                                                  # If is a wwtp with no network
            liftTrenchDepth = 1                                                                         # lift trench

        if liftTrenchDepth == 1:
            if pumpFlow == True:                                                                        # If trench depth would be too deep
                nodesCopy[posToNode][10] = Zmin                                                         # Set new trench depth
                newHeightDiff = nodesCopy[posToNode][10] - trenchFrom                                   # Calculate height difference which needs to be pumped
                if newHeightDiff <= 0:
                    newHeightDiff = nodesCopy[posToNode][10] - ZToDepth          
                pumps = addPump(pumps, fromID, newHeightDiff, nodes)                                    # Add pump
            else:                                                                                       # free flow is possible and trench depth within limits
                if ZToDepth > Zmin and hTo < hFrom: 
                    nodesCopy[posToNode][10] = Zmin                                                     # Set new trench depth
                else:
                    nodesCopy[posToNode][10] = ZToDepth                                                 # Set new trench depth
                    if nodesCopy[posToNode][10] > Zmin:                                                 # If trench depth is above minium possible trench depth 
                        nodesCopy[posToNode][10] = trenchFrom - minTD                                   # Set new trench depth
        else:            
            if ZToDepth < hTo - maxTD:                                                                  # As trench lift is not possible, depth is not changed
                newHeightDiff = nodesCopy[posToNode][10] - trenchFrom                                   # Calculate height difference which needs to be pumped
                if newHeightDiff <= 0:                                                                  # trenchFrom is higher than nodesCopy[posToNode][10]. Might occur because similar heights
                    newHeightDiff = nodesCopy[posToNode][10] - ZToDepth                                 # Even though there is a slope, not steep enough. Pump only the height it would need for free flow
                pumps = addPump(pumps, fromID, newHeightDiff, nodes)                                    # Add pump
            else:                                                                                       # free flow is possible and trench depth within limits
                if trenchFrom == trenchTo:                                                              # Is the case in merging wwtps
                    nodesCopy[posToNode][10] = ZToDepth                         
                else:                                                                                   # Don't change trench depth                                                    
                    if ZToDepth < trenchTo:                                     
                        nodesCopy[posToNode][10] = ZToDepth                                             # Set new trench depth
    else:    
        if pumpFlow == True:                                                                            # has no inflowing node. Free flow is not possible
            nodesCopy[posToNode][10] = Zmin                                                             # Change trenchDepth of toNode to minTD as the water is pumped
            newHeightDiff = nodesCopy[posToNode][10] - trenchFrom                                       # Calculate height difference which needs to be
            if newHeightDiff <= 0:
                newHeightDiff = nodesCopy[posToNode][10] - ZToDepth 
            pumps = addPump(pumps, fromID, newHeightDiff, nodes)                                        # Add pump
        else:  
            if ZToDepth > Zmin and hTo < hFrom:                                                         # free flow is possible
                nodesCopy[posToNode][10] = Zmin                                                         # Set new trench depth
            else:                                                                                       # free flow is not possible
                nodesCopy[posToNode][10] = ZToDepth                                                     # Change trenchDepth of toNode
    return nodesCopy, pumps

def addPump(pumps, toID, pumpHeightDiference, nodes):
    """
    This function either updates a pump at a node or inserts a pump at a node
    
    Input Arguments: 
    pumps                  --    List of pumps
    toID                   --    ID To Node
    pumpHeightDifference   --    Height difference which needs to be pumped
    copypns                --    Nodes
    
    Output Arguments:
    pumps                  --    list with updated or new pumps    
    """
    newPump = False   
    
    # Calculate flow to be pumped
    for i in nodes:
        if i[0] == toID:
            pumpFlow = i[4] + i[8]
            break
    
    for pmp in pumps:                                                        # Check if pump at this node already exists  
        if pmp[0] == toID:                                                   # Pump is found
            pmp[1], pmp[2], newPump = pumpFlow, pumpHeightDiference, True    # replace flow and pumping height
            break   

    if newPump == False:  
        if pumpHeightDiference <= 0 or pumpFlow <= 0:
            raise Exception("ERROR ADDING PUMP:" + str(pumpHeightDiference) + "  "  + str(pumpFlow))
        pumps.append([toID, pumpFlow, pumpHeightDiference, 0])              # #Add new pump. Flow is not corrected and gets chaned later!
    return pumps 

def getNoLoopPath(sewers, path):
    """
    Test for loops in path. If found, replace path. 
    
    Input Arguments: 
    sewers                --    Sewer Network
    path                  --    Path
    
    Output Arguments:
    newPath               --    Path with no loops
    """
    if len(path) <= 1:                  
        return path                     # If no path but only one connecting edge. Loops are not checked
    else:                               # If first the a point in the TONETWORK Is reachd, calculate path on given Network!
        flowToId = path[-1][1][0]       # flow
        flowFromId = path[0][0]         # flow
        part_A = [path[0]]              # List with FROM-Network             
        ExitGivenNetwork = 0            # Used for loop check
        loopFound = 0                   # if there is a loop in the path of the fromNetwork
        loopNode = None                 # Node from where to loop is detected
        firstNodeInToNetwork = None                 
        cnter = 1
        b = path[1][0]
        
        # Test if first edge is not in network
        try:                                    
            if sewers[flowFromId][0] != b:
                if sewers[b][0] != flowFromId:
                    ExitGivenNetwork = 1                # Initial edge is not in network
        except KeyError:
            ExitGivenNetwork = 1
            
        # Iterate path
        for node in path[1:]:
            cnter += 1
            part_A.append(node)
            exitKrit, ide, ide1 = 0, node[0], node[0]
            
            if ExitGivenNetwork == 1:                   # If street was left already in first node, test for loops
                exKit = 0
                                                        # Test if network is entered again. 
                while exKit == 0:                       #Test if FROMID is reached. If yes, calculate path from here to wtpFrom
                    try:
                        new1 = sewers[ide1][0]
                        if new1 == ():
                            if ide1 == flowFromId:
                                correctTO = sewers[node[0]]  
                                firstNodeInToNetwork = [node[0], [correctTO[0], correctTO[1]]]  
                                loopFound = 1                                                   # Loop was found
                                loopNode = node[0]                                              # Node which links to destination wwtp
                                break
                            else:
                                break
                        ide1 = new1
                    except KeyError:
                        exKit = 1
            while 1:                                                                            # Test if toWWTP is reached. If yes, the closest toNetwork Point is found   
                try:                                                                            # If first edge is not in Network  
                    new = sewers[ide]
                    if new[0] == ():
                        if ide == flowToId:
                            correctTO = sewers[node[0]]
                            firstNodeInToNetwork = [node[0], [correctTO[0], correctTO[1]]]
                            del part_A[-1]
                            part_A.append(firstNodeInToNetwork)
                            exitKrit = 1
                            break
                        else:                                                                   # other ID reached --> Network in Between!
                            ExitGivenNetwork = 1
                            break
                    ide = new[0]
                                    
                except KeyError:
                    ExitGivenNetwork = 1
                    if cnter == len(path):                                                      # if last edge, check
                        if node[1][0] == flowToId:
                            firstNodeInToNetwork = [(), 0]
                            exitKrit = 1
                    break    
            if exitKrit == 1:
                break 
        part_A = replaceLoopPath(sewers, loopFound, loopNode, part_A)                           # If loop is found, replace path before loop node with path along existing network      
        if part_A[-1][1][0] == flowToId and part_A[0][0] == flowFromId:                         # If already end and starting node are reached
            return part_A  
        else:
            part_B = replaceCorretPath(sewers, path, loopNode, firstNodeInToNetwork)            # Calculate path to WWTO from firstNodeInToNetwork and replace street path in order to prevent loops
            newPath = mergeLists(part_A, part_B)                                                # Merge pathlists A & B     
            return newPath 

def primCalculations(nodes, PN, allPopNodesOntheWay): 
    """
    This functions calculateds the prim distances.
    
    Input Arguments: 
    nodes                 --    nodes
    PN                    --    PRIM Distances
    allPopNodesOntheWay   --    
    rasterSize            --    Raster Size
    rasterPoints          --    Raster Points

    Output Arguments:
    PN                    --    Updated PRIM Distances
    """              
    for a in allPopNodesOntheWay:
        idp1, p1, _, _, _, _ = getPns(a, nodes)
        for i in PN:
            aktuelleDistanz, idp0 = i[0], i[2]                               # weighted distance, id0
            if idp0 != idp1:                                                 # If not distance to itself
                distanz3d, _, _ = distanceCalc3d(p1, (i[3], i[4], i[5]))     # Calculate euclidian distance
                if distanz3d < aktuelleDistanz:                              # Check if EMST-distance to the new node is smaller. If so, replace the information
                    i[0], i[1] = distanz3d, idp1
    return PN

def getPointsNotInNetwork(pathBetweenWWTPs, nodesFromNetwork, nodesToNetwork):
    """
    This functions gets all points not in a network.
    
    Input Arguments: 
    pathBetweenWWTPs      --    path
    nodesFromNetwork      --    nodes in fromnetwork
    nodesToNetwork        --    nodes in tonetwork
    
    Output Arguments:
    notInaNetwork         --    Nodes not in a network
    """
    notInaNetwork = []
    for i in pathBetweenWWTPs[1:-1]:
        if i not in nodesFromNetwork and i not in nodesToNetwork:
            notInaNetwork.append(i)
    return notInaNetwork

def checkIfWWTPwereConnected(path, wwtpsIterate, WWTPs):
    """
    This functions gchecks if by connecting a wwtp another wwtp was connected or not. If yes, this wwtp gets deleted from list of wwtps.
    
    Input Arguments: 
    path                --    path between wwtps
    wwtpsIterate        --    wwtps to iterate
    WWTPs               --    list with WWTPs
    
    Output Arguments:
    WWTPs               --    Updatet WWTPs
    wwtpsIterate        --    Updatet wwtpsIterate
    """
    for entry in path[1:-1]:
        for wtp in wwtpsIterate:
            if wtp[0] == entry:
                WWTPs = delEntry(WWTPs, entry)
                wwtpsIterate = delEntry(wwtpsIterate, entry)
                break  
    return WWTPs, wwtpsIterate

def assignInitialFlow(nodes, path, initialFlow, position):
    """
    This functions assigns the initial flow to a wwtp 
    
    Input Arguments: 
    nodes                   --    nodes
    path                    --    path between wwtps
    initialFlow        --    initial flow of wwtp in from network
    position                --    Position in list
    
    Output Arguments:
    nodes                  --    Updated edges
    """
    for i in nodes:
        if i[0] == path[position]:
            i[4] = initialFlow - i[8]  # As wtp was included
            break 
    return nodes

def getPath(pathWithDistances):
    """
    This functions gets from a path with distances only the ID-Path.
    
    Input Arguments: 
    pathWithDistances     --    Path with ID and distances
    
    Output Arguments:
    liste                 --    Path with IDs
    """
    liste = []
    for i in pathWithDistances:
        liste.append(i[0])
    return liste                    

def getAggregatedNodesinListWWTP(WWTPs, pipeNetwork, aggregatedNodes):
    """
    This function makes a breath search for each wwtp to get nr of aggregated nodes to listWWTPs.
    
    Input Arguments: 
    WWTPs                             --    list with WWTPs
    pipeNetwork                       --    Sewer pipe network
    
    Output Arguments:
    listWWTPwithAggregatedNodes       --    List with wwtp where the nr of aggregated nodes is added
    """
    listWWTPwithAggregatedNodes = []  # Form: ID, total Flow, total nr of aggregated nodes
    for wwtp in WWTPs:
        allNodes = breathSearch(wwtp[0], pipeNetwork)
        nrOfNodes = 0           # List to store only inahbited nodes

        # Because of arch points which are uninhabited
        for node in allNodes:
            for i in aggregatedNodes:
                if i[0] == node:
                    nrOfNodes += 1
                    break
        listWWTPwithAggregatedNodes.append([wwtp[0], wwtp[1], nrOfNodes])

    return listWWTPwithAggregatedNodes

def checkifInWWTPs(WWTPs, ID):
    ''' 
    This functions check if there exists a WWTP with the iD
    Input: 
    WWTPs     -    List with wwtps
    ID        -    ID
    
    Ouput:
    criteria  -    Criteria wheter WWTP exists or not
    '''
    criteria = False
    for i in WWTPs:
        if i[0] == ID:
            criteria = True
            break 
    return criteria

def checkIfPathAreIdentical(tooLongPath, waytoclosestWWTP, network_ohneArchPath, network):
    ''' This function checks if two path are identical. If no, remove all not itentical nodes from network.
    
    Input: 
    tooLongPath                 -    List with wwtps
    waytoclosestWWTP  -    ID
    network_ohneArchPath        -    Network without added archPath
    network                     -    Network
    
    Ouput:
    network                     -    
    ''' 
    counter, notIdenticalNodes = 0, []

    if len(tooLongPath) > len(waytoclosestWWTP):
        weiterCheck = 1
    if len(tooLongPath) < len(waytoclosestWWTP):
        weiterCheck = 2
    if len(tooLongPath) == len(waytoclosestWWTP):
        weiterCheck = None
                   
    for i in tooLongPath:
        try:
            if i == waytoclosestWWTP[counter]:
                counter += 1
            else:
                if i not in network_ohneArchPath:
                    notIdenticalNodes.append(i)
        except IndexError:
            break
                  
    # If tooLongPath is longer, search remaining not check nodes where they were not in the network
    if weiterCheck == 1:
        for i in tooLongPath[len(waytoclosestWWTP):]:
            if i not in network_ohneArchPath:
                notIdenticalNodes.append(i)
                         
    # Remove all notIdenticalNodes as they are not connected (remove from network)
    for i in notIdenticalNodes:
        for z in network:
            if network[z] != ():
                if network[z][0] == i:
                    del network[z]
                    break
                          
    for i in notIdenticalNodes:
        try:
            del network[i]
        except KeyError:
            _ = 0  # already deleted    
    return network

def getFlowtoCalculateRC(allPopNodesOntheWay, flowFrom, flowTo, nodes):
    '''
    This function checks, wheter the FROM Network or the TO Network. Then is defineds the one with 
    the smaller flow as decentral.
    
    Next it iterates over all populated WWTPs along the connection path (except for the first node).
    
    Input:
    allPopNodesOntheWay    -    All populated nodes along a path
    fromFlowA              -    Flow of From Node
    fromFlowB              -    Flow of From Node
    
    fromToA                -    Flow of To Node
    fromToB                -    Flow of To Node
    
    Output:
    summedFlowDecentralWWTP    -     All flow of the decentral plants
    '''
    
    summedFlowDecentralWWTP = 0
                                
    if flowFrom < flowTo:                                   # From Network is smaller and considered decentral
        allPopNodesOntheWay = allPopNodesOntheWay[::-1]     # Invert
        for decentralizedWTP in allPopNodesOntheWay[:-1]:   # Last one is not considered as it is calculated just above
            for i in nodes:
                if i[0] == decentralizedWTP:
                    summedFlowDecentralWWTP += i[8] + i[4]
                    break 
    else:                                                   # To Network is smaller and considered decentral
        for decentralizedWTP in allPopNodesOntheWay[:-1]:   # Last one is not considered as it is calculated just above
            for i in nodes:
                if i[0] == decentralizedWTP:
                    summedFlowDecentralWWTP += i[8] + i[4]
                    break 

    return summedFlowDecentralWWTP
                            
                            
def getFullHypotheticalCosts(aggregatetPoints, WWTPs, sewers, EW_Q, wwtpLifespan, interestRate, pumps, pumpYears, pricekWh, edgeList, nodes, stricklerC, discountYearsSewers, operationCosts, f_SewerCost, fc_wwtpOperation, fc_wwtpReplacement):
    '''
    This function calculates total system systems and assumes decentralized solution for all not yet considered nodes.
    
    Input:
    aggregatetPoints    -    All nodes to consider
    WWTPs               -    Waste water treatment plants
    sewers              -    Sewer Network
    EW_Q                -    Waste water per EW
    wwtpLifespan        -    lifespan of wwtp
    interestRate        -    interst rate
    pumps               -    List with pumps
    pumpYears           -    Pumping year
    pricekWh            -    price per kWh
    edgeList            -    List with edges
    nodes               -    List with nodes
    stricklerC          -    Strickler coefficient
    discountYearsSewers -    Life span sewers
    operationCosts      -    operation costs sewers
    f_SewerCost         -    cost factor sewer
    fc_wwtpOperation    -    cost factor wwtp
    fc_wwtpReplacement  -    cost factor wwt
    
    Ouput:
    degCen              -    Z
    degCenWeighted      -    Z weighted
    fullCosts           -    Total costs
    listWWTPwithAggregatedNodes    -    All hypothetical WWTP with correct flow for z caluclations
    '''
    hypotheticalWWTP = fastCopy(WWTPs)
    connectedNodes = []                 # nodes to store in connected nodes 

    # Get ID of connected nodes
    for f in sewers:
        connectedNodes.append(f)
        if sewers[f][0] != () and sewers[f][0] not in connectedNodes:
            connectedNodes.append(sewers[f][0])

    # All not yet considered nodes are turned into a WWTP
    for i in aggregatetPoints:
        if i[0] not in connectedNodes:
            hypotheticalWWTP.append([i[0], i[8]])  # Add flow from itself
            
    # Classical Definition of Z
    sources, sinks = len(aggregatetPoints), len(hypotheticalWWTP)
    degCen = round((float(sources) - float(sinks)) / float(sources), 4)

    
    # Weighted Definition of Z
    listWWTPwithAggregatedNodes = getAggregatedNodesinListWWTP(hypotheticalWWTP, sewers, aggregatetPoints) # Get aggregated nodes
    sumFlow, weightedTerm = 0, 0
          
    for i in listWWTPwithAggregatedNodes:
        sumFlow += i[1]
        weightedTerm += float((i[1]/i[2]))

    degCenWeighted = (sumFlow - weightedTerm) / sumFlow
    
    # Read out only the network points
    _, flowPoints = readOnlyNetwork(nodes, sewers) 
    
    # Calculate final costs of whole system
    completePumpCosts, completeWWTPCosts, completePublicPipeCosts = calculatetotalAnnuities(hypotheticalWWTP, EW_Q, wwtpLifespan, interestRate, pumps, pumpYears, pricekWh, sewers, flowPoints, edgeList, nodes, stricklerC, discountYearsSewers, operationCosts, f_SewerCost, fc_wwtpOperation, fc_wwtpReplacement)     

    fullCosts = completePumpCosts + completeWWTPCosts + completePublicPipeCosts
    return degCen, degCenWeighted, fullCosts, completePumpCosts, completeWWTPCosts, completePublicPipeCosts, listWWTPwithAggregatedNodes
                  
def SNIP(OnlyExecuteMerge, outListFolder, runNr, nodes, anteilDaten, streetNetwork, startnode, edgeList, streetVertices, rasterSize, buildPoints, buildings, rasterPoints, inParameter, aggregatetPoints):
    """
    SNIP Algorithm
    
    Input Arguments
    OnlyExecuteMerge       -    Criteria wheter only the Expansion Module is executed
    outListFolder          -    Path to store results
    runNr                  -    Shows wheter expansion module is first time exectued (1) or used in merging module (2)
    nodes                  -    Nodes
    anteilDaten            -    Number of connected points
    streetNetwork          -    Intermediate Results
    startnode              -    List with WWTPs
    edgeList               -    Edges
    streetVertices         -    Street Network
    rasterSize             -    Raster Size
    buildPoints            -    Aggregated Sources
    buildings              -    Buildings
    rasterPoints           -    Raster Points
    inParameter            -    Parameters for SNIP
    writeOutList           -    Intermediate Results
    
    Output Arguments
    ExpansionTime, MergeTime                                                           -    Timers
    final_Network                                                                      -    Network
    flowPoints                                                                         -    Nodes with flow
    WWTPs                                                                              -    WWTPS
    final_wwtps                                                                        -    WWTPs
    final_Pumps                                                                        -    Pumps
    edgeList                                                                           -    list with costs
    completePumpCosts, completeWWTPCosts, completePublicPipeCosts, totalSystemCosts    -    costs
    """                             

    # Input parameters
    minTD = inParameter[0]                                  # [meter] Min trench depth
    maxTD = inParameter[1]                                  # [meter] Maximum trench depth
    minSlope = inParameter[2]                               # [%] Criteria of minimum slope without pumps needed
    f_merge = inParameter[3]                                # Factor do determine how the OST are reconnected. Puts the size in relation ot the distance. If large, the size gets more important
    resonableCostsPerEW = inParameter[4]                    # [CHF] Reasonable Costs per EW
    neighborhood = inParameter[5]                           # How large the neibhourhood is for searching paths
    streetFactor = inParameter[6]                           # How much the path follows the roads
    pricekWh = inParameter[7]                               # Electricity prices
    pumpYears = inParameter[8]                              # Number of years the pummp needs to pump [years]. Needed in order to compare lifespan of pipes to pumping costs
    discountYearsSewers = inParameter[9]                    # [year] Years the pips will survive
    interestRate = inParameter[10]                          # [%] Interest rate
    stricklerC = inParameter[11]                            # [m3/s] Strickler coefficient
    EW_Q = inParameter[12]                                  # [l] conversion factor from liter to EW 162 Liter 
    wwtpLifespan = inParameter[13]                          # [year] how long the wwtp operates
    operationCosts = inParameter[14]                        # [CHF / meter] Average operation costs per meter 
    
    arcpy.AddMessage(" ")
    arcpy.AddMessage("All Paramterers")
    arcpy.AddMessage("----------------")
    arcpy.AddMessage("minTD: " + str(minTD))
    arcpy.AddMessage("maxTD: " + str(maxTD))
    arcpy.AddMessage("f_merge: " + str(f_merge))
    arcpy.AddMessage("resonableCostsPerEW: " + str(resonableCostsPerEW))
    arcpy.AddMessage("neighborhood: " + str(neighborhood))
    arcpy.AddMessage("streetFactor: " + str(streetFactor))
    arcpy.AddMessage("pricekWh: " + str(pricekWh))
    arcpy.AddMessage("pumpYears: " + str(pumpYears))
    arcpy.AddMessage("discountYearsSewers: " + str(discountYearsSewers))
    arcpy.AddMessage("interestRate: " + str(interestRate))
    arcpy.AddMessage("stricklerC: " + str(stricklerC))
    arcpy.AddMessage("EW_Q: " + str(EW_Q))
    arcpy.AddMessage("wwtpLifespan: " + str(wwtpLifespan))
    arcpy.AddMessage("operationCosts: " + str(operationCosts))
    arcpy.AddMessage(" ")
    arcpy.AddMessage(" ")


    # Initialization of parameters
    firstIteration = 1                                      # initial parameter for first iteration
    initialPN = []                                          # used for PRIM
    sewers = {}                                             # Graph containing the network which is beeing built
    sewers_Current = []                                     # Graph containing all nodes of the current netowrk. The OST are as well included.
    PN = [(0, (), startnode, 0, 0, 0, 1)]                   # List for prims algorithm with all remaining nodes to be connected
    WWTPs = []                                              # List containing all wwtps. [(WTP_ID, summedFlow, fromNetworkPointbelow, length of edge fromNetworkPointbelow),()]
    pumps = []                                              # List with all pumps
    expansion = 1                                           # Abort criteria   
    totalSystemCosts = []                                   # List storing total system costs and Z
    hypoWWTPcorrectFlow = []                                # Initial
    
    #pumpInvestmentCosts = inParameter[15]                  # [CHF] Investment Costs pumps
    f_topo = inParameter[16]                                # [-] Weighting factor for the DEM graph creation of the topography
    f_SewerCost = inParameter[17]                           # [-] cost factor
    fc_wwtpOperation  = inParameter[18]                     # [-] cost factor
    fc_wwtpReplacement = inParameter[19]                    # [-] cost factor
                                 
    # Initial calculations
    hypoZWeighted = 0                                       # Initial weighted Z value
    hypoCostsOld = 0                                        # Initialization Hypothetical costs
    reActivationEM = 0                                      # If EM is first time loade == 0, otherwhise 1 (later for reconnection)
    firstMergeCrit = 1                                      # Criteria wheter in first or subsequent mering module execusion
    Zreached = 0                                            # Criteria wheter the requested Z is reached

    if OnlyExecuteMerge == 0:
        #arcpy.AddMessage("NODE: " + str(startnode))
        #arcpy.AddMessage("WWTPs: " + str(nodes))
        WWTPs = addOriginTolistWWTPs(nodes, WWTPs, startnode) # Add wwtp to listwwtps
        
    # Additional factors for aborting the SNIP at a certain Z-Value
    iterativeCostCalc = 0                                     # Wheter iteratiely the costs should be calculated   
    ignoreOnlyExecuteMerge = 0
    afterMM = False
    
    forceZ = False                                             # Criteria wheter SNIP is aborted at the optimum or articially forced to connect further
    ZtoReach = 0.99                                            # In case artifically connection is enforced SNIP terminates at this value (must be smaller than 1)

    if forceZ == True:  
        iterativeCostCalc = 1   
            
    # Create folder to store results 
    txtResultPath = outListFolder + "ResultAsTxt" + "/"     
    if not os.path.exists(txtResultPath):
        os.mkdir(txtResultPath)
  
    # Expansion Module is activated
    while expansion == 1:
        arcpy.AddMessage("Start expansion module...")

        
        if PN == []:                                          # Exit in case PN is empty 
            expansion = 0
            continue

        if OnlyExecuteMerge == 0 and ignoreOnlyExecuteMerge == 0: 
            while len(PN) > 0 or firstIteration == 1:                                                       
                if ZtoReach > hypoZWeighted:                                                                # Used to abort in case a certain z-value is reached
                    PN, FROMNODE, TONODE, weightFactorDijkstra, realDistance = getClosestNode(PN)           # Select next node to connect
                    
                    if TONODE not in sewers:   
                        sewerBeforeIteration = dict(sewers)                                                 # Used in order that only not connected nodes can be wwtps
             
                        # Calculate hypothetical costs: Calculate costs of current network, wwtps and if hypothetically all other wwtps wouldn't be connected and have a decentral conncetion. 
                        if iterativeCostCalc == 1 and runNr == 1:
                            hypoZ, hypoZWeighted, hypoCosts, hypoPumpCosts, hypoWWTPCosts, hypoPublicPipeCosts, hypoWWTPcorrectFlow = getFullHypotheticalCosts(aggregatetPoints, WWTPs, sewers, EW_Q, wwtpLifespan, interestRate, pumps, pumpYears, pricekWh, edgeList, nodes, stricklerC, discountYearsSewers, operationCosts, f_SewerCost, fc_wwtpOperation, fc_wwtpReplacement)
                            
                            if reActivationEM == 0:
                                if hypoCosts < hypoCostsOld:
                                    totalSystemCosts.append([hypoZ, hypoZWeighted, hypoCosts, hypoPumpCosts, hypoWWTPCosts, hypoPublicPipeCosts])   # Store costs in list
                                    
                            hypoCostsOld = hypoCosts
                            
                        if len(sewers) > 0:
                            firstIteration = 0
                        else:
                            sewers = initialSetNetwork(sewers, TONODE, FROMNODE, realDistance)                              # Initial Sewer Network
  
                        # Path finding module (PFM)
                        if firstIteration == 0:
                            idp0, p0, _, _, _, _ = getPns(FROMNODE, nodes)                                                     # get node information
                            idp1, p1, _, _, forceConnection, _ = getPns(TONODE, nodes)                                         # get node information
                            _, slopeMST, heightDiff = distanceCalc3d(p0, p1)                                                # calc slope of straight distance
                            try:                                                                                            # Try to find path on street with Dijkstra Algorithm
                                archPath, distStartEnd, _ = dijkstra(streetNetwork, idp0, idp1, heightDiff)                 # Djikstra   
                                Djikstradistancce = distStartEnd * weightFactorDijkstra                                     # weight distance
                                streetConnection = 1      
                                #arcpy.AddMessage("Path was found along the street..." + str(archPath))                                      # criteria whether the djikstra distance or MST distance was    
                            except KeyError:
                                #arcpy.AddMessage("no street was found. This could be because Djikstra failed or the buillding is too far away from the street netowrk.")
                                streetConnection = 0  
         
                            # Compare distance along the street and straight distance       
                            if streetConnection == 1:  # If a StreetNetwork exists
                                # Because the euclidian distance is calculated with the original heights (not considering the weighted from the aggregation) it might happen that the Djikstradistance is shorter than the Euclidian distance
                                if Djikstradistancce < realDistance:
                                    Djikstradistancce = realDistance
                                
                                # In case the MM has been executet one a* is not considered anymore (because there is always a connection along the street network)
                                if distStartEnd >= streetFactor * realDistance and afterMM != True:           # Fstreet Factor
                                    changeStreetGraph, streetConnection = 1, 0                                # Path along the terrain
                                else: 
                                    changeStreetGraph = 0                                                     # Path along the street

                            # If there is no street, the a* algorithm is used to find shortest path on DEM. If there is no such path, a MST is returned.
                            if streetConnection == 0: 
                                #arcpy.AddMessage("Try finding a path along the terrain (a*)...")
                                changeStreetGraph = 1                                           
                                archPathMST, boundingCandidates = aStar(rasterSize, rasterPoints, buildPoints, p0, p1, idp0, idp1, neighborhood, f_topo)
                                nodes = addDEMPntstoNodes(nodes, archPathMST, boundingCandidates, FROMNODE, TONODE, minTD)  # Add new DEM-Points to nodes
                                
                                if len(archPathMST) == 0:
                                    #arcpy.AddMessage("NO connection WAS FOUND" + str(streetConnection))                                                                  # If no A-Start Path was found, use MST
                                    if streetConnection == 1:
                                        archPathMST = archPath                                                              # But if there is a street connection, take the street distance
                                    else:
                                        archPathMST = [[idp0, [idp1, realDistance, slopeMST]]]    
    
                            # Change Street Graph and append edges of a* algorithm
                            if changeStreetGraph == 1:
                                wayProperty = 0 
                                edgeList, streetNetwork = addEdgesUpdateStreetNetwork(wayProperty, archPathMST, nodes, boundingCandidates, edgeList, streetNetwork)
                       
                        # Path calculations
                        if firstIteration == 0:
                            if changeStreetGraph == 0 :                                                                     # Insert archPath in network
                                sewers = appendToNetwork(sewers, archPath)                                                  # Append path from this node to closest WWTP from this point in P
                                tooLongPath = readPath(archPath, TONODE)                                                    # path from FROMNODE to TONODE  
                            else:
                                sewers = appendToNetwork(sewers, archPathMST)                                               # Append path from this node to closest WWTP from this point in P
                                tooLongPath = readPath(archPathMST, TONODE)                                                 # path from FROMNODE to TONODE
                                
                            pathtonearestWTP, closestARAtraditionell = getclosestWTP(TONODE, WWTPs, sewers)                 # Find Way to closest WWTP in Network P
                            pathtonearestWTPInvert = pathtonearestWTP[::-1]                                                 # Inverse Path to closest WWTP   
                            pathtonearestWTPswap = appendDistances(pathtonearestWTPInvert, sewers)                          # Find way to cloeset WTP in inverted Network 
                            pathNearWTP = appendDistances(pathtonearestWTP, sewers)                                         # Path to closest WWTP. The distances and slope get appended to the path OLD  
                            pathNearWTPInvert = pathNearWTP[::-1]                                                           # Invert pathOrigin to get correct flow to origin 
                                
                            # Check if the adding of the new path creates loops
                            if len(sewerBeforeIteration) > 1:  
                                pathNearWTPInvert, loopFound = loopTest(pathNearWTPInvert, sewerBeforeIteration) 
                                if loopFound == True:                                                                       # If loop is found
                                    pathNearWTP = pathNearWTPInvert[::-1]                                                   # Change pathNearWTP as a loop was found
                                    pathtonearestWTP = getPath(pathNearWTP)
                                    pathtonearestWTPInvert = pathtonearestWTP[::-1]
                                    pathtonearestWTPswap = appendDistances(pathtonearestWTPInvert, sewers)
                                    sewers = checkIfPathAreIdentical(tooLongPath, pathtonearestWTP, sewerBeforeIteration, sewers) # Compare IF pathtonearestWTP is identical to pathtonearestWTP. If not identical, the connection via an already existing pipe is better. (oterhweilse multiple flows)
                                
                            # Check if in archPath Nodes were added to sewer which are not necessary and remove them.
                            if changeStreetGraph == 0:
                                sewers = delnodesNotUsedInSewer(sewers, pathtonearestWTP, archPath, sewerBeforeIteration)
                            else:
                                sewers = delnodesNotUsedInSewer(sewers, pathtonearestWTP, archPathMST, sewerBeforeIteration)
                                
                            pathToNearestWTPwithDistances = writeInEdges(pathNearWTP)                                                # Path to closest WWTP with distances
                            pathToNearestWTPswapwithDistances = InvertandswapID(pathToNearestWTPwithDistances)                       # Inverse path with distances
                            allPopNodesOntheWay, flowAllArchPathWWTPs = getPopNodes(pathNearWTPInvert, nodes, sewerBeforeIteration)  # All source on the path to the wwtop 
                            pathSubNetworkToClosestWWTP = getPathToWTPSubentwork(pathtonearestWTP, sewerBeforeIteration)             # Get path to WWTP in Subnetwork
                            
                            # ====================================================
                            # Option  Module (OM) & Cost module (CM)
                            # ====================================================      
                            P_A3 = dict(sewers)
                            nodesA1, nodesA3 = fastCopyNodes(nodes), fastCopyNodes(nodes) 
                            swapCriteria, dezentralCriteria = 0, 0                                                                  # if swap takes places ((yes or no), if decentral (yes or no)
                            pumpsA1, pumpsA3 = fastCopy(pumps), fastCopy(pumps)
                            P_A3 = invertFlowToNearestWTP(pathtonearestWTPswap, P_A3, sewers)                                       # Invert flow along the way to the nearest WWTP

                            # --------
                            # Option 1
                            # --------
                            
                            # Option 1 - Sewer costs
                            nodesA1 = updateFlowA1(nodesA1, pathNearWTPInvert, allPopNodesOntheWay)                                                             # Correct flow along the path.
                            nodesA1, inflowNodesA1 = correctTD(nodesA1, pathtonearestWTP[:-1], minTD, WWTPs, maxTD, minSlope, sewers)                           # Set all trench Depth except inflow nodes to minimum trench depth              
                            nodesA1, pumpsA1, edgeListI = changeTD(nodesA1, edgeList, pumpsA1, pathToNearestWTPswapwithDistances, maxTD, minSlope, inflowNodesA1, sewers, minTD)
                            pipeCostA1, edgeListI = costToWTP(pathtonearestWTP, edgeListI, nodesA1, pumpsA1, minTD, TONODE, sewerBeforeIteration, discountYearsSewers, interestRate, stricklerC, operationCosts, f_SewerCost)  # Pipe costs

                            # Option 1 - WWTPs costs
                            flowWWFrom = getFlowWWTP(WWTPs, closestARAtraditionell)                                                                             # WWTP costs: Additional costs of building a larger WWTP. Flow of closestARAtraditionell
                            WWTPcostsA1 = costWWTP((flowWWFrom + flowAllArchPathWWTPs), EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement) # Costs of new WWTPs with original Flow, flow of connecting node and flow of all wwtps on archPaths
                            
                            # Option 1 - Pumping costs
                            pumpCostA1, pumpCostWholePeriodI = costPump(pathtonearestWTP, pumpsA1, pricekWh, pumpYears, interestRate)                           # Calculate annual pumping costs

                            # --------
                            # Option 2
                            # --------
                            
                            # Always the smaller network is considered to be the "decental" network used for calculating the reasonable costs 
                            for fl in nodes:
                                if fl[0] == pathNearWTPInvert[-1][0]:
                                    flowTo = fl[8] + fl[4]
                                    break
                            for fl in nodes:
                                if fl[0] == pathNearWTPInvert[0][0]:
                                    flowFrom = fl[8] + fl[4]
                                    break
                                
                            summedFlowDecentralWWTP = getFlowtoCalculateRC(allPopNodesOntheWay, flowFrom, flowTo, nodes)              

                            # Option 2 - WWTPs costs # TODO CHANGE FLOWTOADD (REMOVE)
                            costDecentralWWTPs = sumCostOfArchWWTPs(allPopNodesOntheWay[:-1], nodes, EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement)  # Calculate costs of new wwtps on added path
                            wtpCostClosestARA = costWWTP(flowWWFrom, EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement)                                             # cost of already build wwwt
                            WWTPcostsA2 = costDecentralWWTPs + wtpCostClosestARA        
                            
                            # Option 2 - Sewer costs
                            pipeCostA2, _ = costToWTP(pathSubNetworkToClosestWWTP, edgeList, nodes, pumps, minTD, TONODE, sewerBeforeIteration, discountYearsSewers, interestRate, stricklerC, operationCosts, f_SewerCost)  # Pipe costs
        
                            # Option 2 - Pumping costs
                            pumpCostA2, pumpCostWholePeriodA2 = costPump(pathSubNetworkToClosestWWTP, pumps, pricekWh, pumpYears, interestRate)  # Pumping costs II
                                     
                            # --------
                            # Option 3
                            # --------
                            # Option 3 - WWTPs costs
                            WWTPcostsA3 = WWTPcostsA1                                      
                                
                            # Option 3 - Sewer costs
                            if len(pathtonearestWTP) < 3:                                                                               # Single edge
                                nodesA3 = updatePathShort(nodesA3, pathtonearestWTPInvert)  
                            else:
                                nodesA3 = updatePathLong(nodes, nodesA3, pathtonearestWTPInvert, sewerBeforeIteration, TONODE)          # Several edges
                            nodesA3, inflowNodesA3 = correctTD(nodesA3, pathtonearestWTPInvert, minTD, WWTPs, maxTD, minSlope, P_A3)    # Set all trench depth except inflow nodes to minimum trench depth  
                            nodesA3, pumpsA3, edgeListIII = changeTD(nodesA3, edgeList, pumpsA3, pathToNearestWTPwithDistances, maxTD, minSlope, inflowNodesA3, P_A3, minTD)
                                 
                            # Option 3 - Pumping costs
                            pumpCostA3, pumpCostWholePeriodA3 = costPump(pathtonearestWTPInvert, pumpsA3, pricekWh, pumpYears, interestRate)  # Pump Costs III

                            # Calculate costs
                            pipeCostA3, edgeListIII = costToWTP(pathtonearestWTPInvert, edgeListIII, nodesA3, pumpsA3, minTD, closestARAtraditionell, sewerBeforeIteration, discountYearsSewers, interestRate, stricklerC, operationCosts, f_SewerCost)  # sum all, use negative slopes
                                
                            # Total option costs
                            totCostsA1 = pipeCostA1 + WWTPcostsA1 + pumpCostA1  # Option I        
                            totCostsA2 = pipeCostA2 + WWTPcostsA2 + pumpCostA2  # Option II
                            totCostsA3 = pipeCostA3 + WWTPcostsA3 + pumpCostA3  # Option III

                            # Calculate cost of cheapest central connection costs.
                            costConnectionA = calculateConnectionCosts(pipeCostA1, pumpCostWholePeriodI, pipeCostA2, pumpCostWholePeriodA2, pipeCostA3, pumpCostWholePeriodA3, WWTPcostsA1, interestRate, wwtpLifespan, discountYearsSewers)
                           
                            resonableCosts = summedFlowDecentralWWTP/EW_Q * resonableCostsPerEW  # Total reasonable costs # Calculate reasonable costs (If positive: central costs are more expensive, if negative: decentral is more expensive)                                                                                                       # 1: Ignore, 0: Don't Ignore 

                            swapCriteria, dezentralCriteria = costComparison(totCostsA1, totCostsA2, totCostsA3, resonableCosts, costConnectionA, 0, forceConnection, 1)  # Compare costs and select option, consider reasonable costs
                            
                            if dezentralCriteria == 0 and forceConnection != 1:
                                #arcpy.AddMessage("OPTION NO Connection")
                                if changeStreetGraph == 0:
                                    WWTPs, intermediateWWTPs = addToListWTPs(archPath, nodes, sewerBeforeIteration, WWTPs)
                                else:  
                                    WWTPs, intermediateWWTPs = addToListWTPs(archPathMST, nodes, sewerBeforeIteration, WWTPs)

                                sewers_Current = appendToSewers(sewers_Current, intermediateWWTPs)                          # Append to sewers_Current 
                                sewers = dict(sewerBeforeIteration)
                                sewers = appendListWTPsToNetwork(sewers, intermediateWWTPs)                                 # Append to network

                                # Add path to edgeList
                                for i in pathToNearestWTPwithDistances:                                                  
                                    fromID, toID = i[0], i[1][0] 
                                    length = i[1][1]                                                         # Read length
                                    IDnew, cord_new, _, _, _, trenchTo = getPns(toID, nodes)
                                    IDold, cord_old, _, _, _, trenchFrom = getPns(fromID, nodes)
                                    newPipeSlope =  (trenchFrom - trenchTo) / length                        # New pipe slope
                                    edgeList = addToEdgeList(edgeList, length, newPipeSlope, IDnew, IDold, cord_new, cord_old)
                            else:
                                if swapCriteria == 1:                                                                       # Switch WWTPs and connect
                                    #arcpy.AddMessage("OPTION SWAP")
                                    sewers = dict(P_A3)
                                    edgeList = fastCopy(edgeListIII)    # Replace nodes by new nodes with new flow  
                                    nodes = fastCopyNodes(nodesA3)
                                    WWTPs = delWWTP(WWTPs, closestARAtraditionell)                                          # Delete wwtps
                                    sewers_Current = appendToSewers(sewers_Current, pathtonearestWTP)                       # Append all newly connected nodes to pathtonearestWTP
                                    pumps = fastCopy(pumpsA3)                                                               # Replace pumps
                                    newAddFlow = getSummedFlow(nodes, TONODE)                                               # Get new flow
                                    WWTPs.append([TONODE, newAddFlow])                                                      # Add new wwtp with changed flow
                                    nodes = removeForceCriteria(nodes, TONODE)                                              # Remove force criteria as is connected
                                else:                                                                                       # Central Connection without WWTP switch
                                    #arcpy.AddMessage("OPTION NO SWAP")
                                    sewers_Current = appendToSewers(sewers_Current, pathtonearestWTP)                       # Append all newly connected nodes to pathtonearestWTP
                                    edgeList = fastCopy(edgeListI) # Replace nodes and pumps
                                    pumps = fastCopy(pumpsA1)
                                    nodes = fastCopyNodes(nodesA1)
                                    WWTPs = updateFlowInWWTP(WWTPs, nodes, closestARAtraditionell)                          # Update list with wwtps                                
                                    nodes = removeForceCriteria(nodes, closestARAtraditionell)

                        PN, initialPN, firstIteration = initialPrimCalc(firstIteration, startnode, nodes, PN, initialPN)    # Initia prim based calculation
                        
                        if len(sewers) > 1:
                            if dezentralCriteria == 0:
                                allPopNodesOntheWay = [TONODE]
                            PN = primCalculations(nodes, PN, allPopNodesOntheWay)                                           # Update prim distances
                else:                                                                                                       # Z-Value reached to abort
                    Zreached, expansion, firstIteration, PN = 1, 0, 0, []
            firstIteration = 0 

        # Merging modul after Expansion Module
        if Zreached == 0:                                                       # Go to Merging Module
            nodes, sewers, pumps, WWTPs, PN, edgeList, totalSystemCosts = mergingModule(OnlyExecuteMerge, forceZ, ZtoReach, hypoZWeighted, firstMergeCrit, aggregatetPoints, nodes, WWTPs, sewers_Current, edgeList, pumps, PN, sewers, streetNetwork, f_merge, minTD, maxTD, minSlope, discountYearsSewers, interestRate, stricklerC, operationCosts, pricekWh, pumpYears, wwtpLifespan, EW_Q, resonableCostsPerEW, f_SewerCost, fc_wwtpOperation, fc_wwtpReplacement, totalSystemCosts, iterativeCostCalc) 
            runNr, firstMergeCrit, reActivationEM = 0, 0, 1                     # Expansion module is finished, As from now on the EM is only reactivated    
            expansion = testExpansion(PN)                                       # Test if there is still expansion needed

    if reActivationEM == 0:                                                     # In case Expansion module is aborted because the wished Z-Values is reached, consider all not yet connected nodes as WWTPs
        final_wwtps = turnNodesIntoWWTP(WWTPs, nodes, aggregatetPoints, sewers) # Add all not yet connected nodes to WWTPs
    else:
        final_wwtps = WWTPsToDraw(WWTPs, nodes)                                 # Write all WWTPs
    final_Network, flowPoints = readOnlyNetwork(nodes, sewers)                  # Read out only the network points
    final_Pumps = readPumpList(pumps, nodes)                                    # Read out all pumps

    # Calculate costs of whole system
    completePumpCosts, completeWWTPCosts, completePublicPipeCosts = calculatetotalAnnuities(WWTPs, EW_Q, wwtpLifespan, interestRate, pumps, pumpYears, pricekWh, final_Network, flowPoints, edgeList, nodes, stricklerC, discountYearsSewers, operationCosts, f_SewerCost, fc_wwtpOperation, fc_wwtpReplacement)
    
    # Store the result as .txt files
    writeTotxt(txtResultPath, "inParameter", inParameter)
    writeToDoc(txtResultPath, "streetNetwork", streetNetwork)
    writeTotxt(txtResultPath, "aggregatetPoints", aggregatetPoints)
    writeTotxt(txtResultPath, "nodes", nodes)
    writeTotxt(txtResultPath, "WWTPs", WWTPs)
    writeTotxt(txtResultPath, "edgeList", edgeList)
    writeTotxt(txtResultPath, "pumps", pumps)
    writeTotxt(txtResultPath, "sewers_Current", sewers_Current)
    writeToDoc(txtResultPath, "sewers", sewers)
    writeTotxt(txtResultPath, "buildings", buildings)
    writeTotxt(txtResultPath, "buildPoints", buildPoints)
    
    # Calculate Z values
    calcZValue(aggregatetPoints, forceZ, iterativeCostCalc, hypoWWTPcorrectFlow, WWTPs)             # Z normal
    
    flowIfSameCheck(final_wwtps, aggregatetPoints)                                                  # Check if flow is lost
    arcpy.AddMessage("SNIP is successfully calculated.")
    arcpy.AddMessage("Costs per iteration:" + str(totalSystemCosts))
    MergeTime, ExpansionTime = 0, 0
                        
    return ExpansionTime, MergeTime, final_Network, flowPoints, WWTPs, final_wwtps, final_Pumps, edgeList, completePumpCosts, completeWWTPCosts, completePublicPipeCosts, totalSystemCosts, buildings, buildPoints, aggregatetPoints

def mergingModule(OnlyExecuteMerge, forceZ, ZtoReach, hypoZWeighted, firstMergeCrit, aggregatetPoints, nodes, WWTPs, sewers_Current, edgeList, pumps, PN, sewers, streetNetwork, f_merge, minTD, maxTD, minSlope, discountYearsSewers, interestRate, stricklerC, operationCosts, pricekWh, pumpYears, wwtpLifespan, EW_Q, resonableCostsPerEW, f_SewerCost, fc_wwtpOperation, fc_wwtpReplacement, totalSystemCosts, iterativeCostCalc):
    '''
    Merging Module
    
    Input Arguments:
    OnlyExecuteMerge     -    Shows wheter summary calculation with different folders: 1 --> yes
    forceZ               -    Criteria to reach full centralization
    ZtoReach             -    Abort Z criteria
    firstMergeCrit       -    Criteria wheter executed the first time or not
    nodes                -    Nodes
    WWTPs                -    WWTPs
    sewers_Current       -    sewer network from expansion module
    edgeList             -    List with edges
    pumps                -    Pumps
    PN                   -    Prim list with distances
    sewers               -    Sewer Network
    streetNetwork        -    Street Network
    f_merge              -    Merging Factor
    minTD                -    Minimum Trench depth
    maxTD                -    Maximum Trench Depth
    minSlope             -    Minimum slope
    discountYearsSewers  -    Discounting years
    interestRate         -    Inerest rate
    stricklerC           -    Strickler coefficient
    operationCosts       -    Operation costs
    pricekWh             -    energy costs
    pumpYears            -    Lifepspan pumps
    pumpInvestmentCosts  -    pump investment costs
    wwtpLifespan         -    wwtp lifespan
    EW_Q                 -    sewage per person per day
    resonableCostsPerEW  -    reasonable costs  
    totalSystemCosts     -    List conting Z and hypothethical costs
    iterativeCostCalc    -    
    
    Output Arguments:
    nodes:                -    Networ nodes
    sewers:               -    Sewers
    pumps:                -    Pumps
    WWTPs:                -    WWTPs
    PN:                   -    Prim distances
    edgeList              -    List with edges
    totalSystemCosts      -    Total System Costs
    '''
    sortedListWWTPs = fastCopy(WWTPs)
    sortedListWWTPs.sort(key=operator.itemgetter(1))    # Sort list along the size (second position in list)
    sortedListWWTPs = sortedListWWTPs[::-1]             # Invert sorting in oder that max values are first
    mergeCostStorage = []                               # Used to store cheapest connetio nin case Z == 1 wants to be achieved
    finishedExpansion = 0
    
    if iterativeCostCalc == 1:
        hypoZOld = totalSystemCosts[len(totalSystemCosts)-1][0] # Current Z value
    else:
        hypoZOld = 0
    
    # As long not merging is exited and there are wwtps to reconnect
    while len(sortedListWWTPs) > 1 and ZtoReach > hypoZWeighted and len(WWTPs) > 2: #last is new
        wwtpsIterate, expansionM = fastCopy(WWTPs), True

        # In case a connection is forced
        if forceZ == True and len(sortedListWWTPs) == 2 or finishedExpansion == 1: # Second Condition is new

            # NEW 19.05.2015
            if len(sortedListWWTPs) == 2:
                sortedListWWTPs = fastCopy(WWTPs)
                finishedExpansion = 1

            if ZtoReach > hypoZWeighted:    
                
                if len(sortedListWWTPs) == 2 or finishedExpansion == 1:
                    #arcpy.AddMessage("Forced furger merging..." + str(hypoZWeighted))
                    mergeCostStorage, cheapestM = getCheapestMerge(mergeCostStorage, WWTPs)                     # get cheapest connection
                    CurrentCheckedWWTP, checkBackConnectionID = cheapestM[1], cheapestM[2]

                else:
                    CurrentCheckedWWTP, checkBackConnectionID = sortedListWWTPs[0], sortedListWWTPs[0][0]       # Current WWTP to check for merge, Get largest wwtp to check
                    del sortedListWWTPs[0] 
            else:                                                                                               # Quit SNIP 
                sortedListWWTPs = []
                continue
        else:
            CurrentCheckedWWTP, checkBackConnectionID = sortedListWWTPs[0], sortedListWWTPs[0][0]               # Current WWTP to check for merge, Get largest wwtp to check
        
        if len(WWTPs) == 1:                                     # If only one WWTP left quit SNIP
            expansionM = False
            break
        
        while len(wwtpsIterate) > 1 and expansionM == True:    # Iterate until all WWTPs are checked   
            sewers_NoCon = dict(sewers)                        # Make copy 
            pumps_noCon = fastCopy(pumps)                      # Make copy   
            nodes_noCon = fastCopyNodes(nodes)                 # Make copy 
            WWTPS_noCon = fastCopy(WWTPs)                      # Make copy 
            PN_noCon = fastCopy(PN)                            # Make copy 
            sortedListWWTPs_noCon = fastCopy(sortedListWWTPs)  # Make copy 
            wwtpToWwtopConncetion = 0                          # Criteria to define wheter wwtps are connected or not                              
                
            while expansionM == True: 
                _, pZero, _, _, forceCentralConnection, _ = getPns(checkBackConnectionID, nodes)  # Select wwtp where node was connected and test if connection is allowed
                if forceCentralConnection == 1:  
                    wwtpsIterate = delEntry(wwtpsIterate, checkBackConnectionID) 
                    expansionM = True

                # Option Selection Module (OSM)   
                nodeIdOpt1, nodeOpt1Cor, nodeIdOpt3, nodeOpt3Cor = connectivityPotential(nodes, WWTPs, checkBackConnectionID, pZero, f_merge)    # Option 1 & Option 3
                nodeIdOpt2, nodeOpt2Cor, pathInClosestNetwork = getClosestNetworkWWTP(nodes, WWTPs, sewers_NoCon, pZero, checkBackConnectionID)  # Option 2
                mergeOptions = storeMergingOptions(nodeIdOpt1, nodeOpt1Cor, nodeIdOpt2, nodeOpt2Cor, nodeIdOpt3, nodeOpt3Cor)         # Store Merging Options in a list
                
                if mergeOptions == []:  # No more possibilities for merging
                    expansionM = False
                    break  

                # In case we want to reach Z == 1
                if ZtoReach > hypoZWeighted and finishedExpansion == 1: #
                    mergeOptions = [cheapestM[3]] 

                if expansionM == True: 
                    iterateMergeOptions = True  
                    while iterateMergeOptions == True:
                        #arcpy.AddMessage("Option to test merge: " + str(mergeOptions))
                        if len(mergeOptions) == 0:
                            iterateMergeOptions = False
                            expansionM = False
                            break
                        
                        for mOpt in mergeOptions:
                            try:                                                                                                      # Search path between WWTps with Djikstra
                                heightDiff = abs(pZero[2] - mOpt[1][2])                                                                          # heightdifference
                                archPathWWTP, distStartEnd, slopeDijkstra = dijkstra(streetNetwork, mOpt[0], checkBackConnectionID, heightDiff)  # Djkstra
                                edgeList = addToEdgeList(edgeList, distStartEnd, slopeDijkstra, mOpt[0], checkBackConnectionID, mOpt[1], pZero)  # If the new streetDistance is not already added in edgeList, add to edgeList

                            except KeyError:
                                raise Exception("There can be no sewer path be determined between the two wwtps.")

                            if nodeIdOpt2 == mOpt[0]:
                                archPathWWTP = mergePathClosestNetwork(pathInClosestNetwork, archPathWWTP, edgeList)     # If closest network, check path again.
                            
                            archPathWWTP = getNoLoopPath(sewers, archPathWWTP)                              # Path is changed in order that there are no loops.      
                            sewers = appendToNetwork(sewers, archPathWWTP)                                  # Append no loop path to sewers
                            inversearchPathWWTP = invertArchPath(archPathWWTP)      
                                
                            # Get parts of path which are already part of a the starting and ending network
                            nodesToNetwork = findInNetworkPath(archPathWWTP, sewers_NoCon, sewers)          # finds part of the From Network which need to be considered in costs calculations
                            nodesFromNetwork = findInNetworkPath(inversearchPathWWTP, sewers_NoCon, sewers) # finds part of the To Network which need to be considered in costs calculations
                            WWTPTO, WWTPFROM = nodesToNetwork[0], nodesFromNetwork[0]                       # WWTPs
                                                                    
                            # Update
                            PathListToPotentialWWTP = InvertandswapID(inversearchPathWWTP)                  # swap ID and invert path
                            pathBetweenWWTPs = readOnlyNodesID(archPathWWTP)                                # Only ID
                            pathBetweenWWTPsInvert = pathBetweenWWTPs[::-1]                                 # Invert abd only ID
                            archPathWWTPInvert = InvertandswapID(archPathWWTP)                              # swap ID and invert path
    
                            # Add ArchPath to edges # new that full path is entered into edgeList
                            for edge in archPathWWTP:
                                newNode, oldNode = edge[0], edge[1][0]
                                IDnew, cord_new, _, _, _, _ = getPns(newNode, nodes)
                                IDold, cord_old, _, _, _, _ = getPns(oldNode, nodes)
                                distance, slope, _ = distanceCalc3d(cord_new, cord_old)
                                edgeList = addToEdgeList(edgeList, distance, slope, IDnew, IDold, cord_new, cord_old)
    
                            #=========================================
                            # Check if networks need to be removedRemove networks
                            # If there are intercrossed networks to remove, remove network nodes from sewers, clear nodes and add again to PN
                            #=========================================                                                                                                              
                            netWorkToRemove, needsNetworkRemoving = findNetworkToRemove(archPathWWTP, nodesToNetwork, nodesFromNetwork, sewers_NoCon, WWTPs)    # Find sub networks to remove
        
                            if needsNetworkRemoving == 1:
                                #arcpy.AddMessage("Network needs to be removed")
                                allNodesToAddToPN = []

                                # Iterate over all wwts which needs to be deleted
                                for wtpToDelet in netWorkToRemove:   
                                    allNodesToDelet = breathSearch(wtpToDelet, sewers)          
                                    nodes = clearFlow(nodes, allNodesToDelet)                   # Clear copypnts flow and pressure
                                    pumps = removePumpsWhereNoNetwork(pumps, allNodesToDelet)   # Remove if there are pumps in between 
                                    sewers = removeInSewers(sewers, allNodesToDelet)            # Update the Ps and delete all network nodes from sewers
                                    WWTPs = delEntry(WWTPs, wtpToDelet)                         # delet WWTP
                                    allNodesToAddToPN.append((wtpToDelet, allNodesToDelet))     # Append PN distances
                                    sortedListWWTPs = delEntry(sortedListWWTPs, wtpToDelet)     # Delete WWTP
                                
                                sewers = appendToNetwork(sewers, archPathWWTP)                  # Add the path between the wwtps to the sewer network
                                                 
                            sewers_B1, sewers_B3 = dict(sewers), dict(sewers)                   # copy of sewers
                            sewers_B1 = insertPathDirection(sewers_B1, PathListToPotentialWWTP) # Add all nodes between wwtps to copy of sewers
                            sewers_B3 = insertPathDirection(sewers_B3, inversearchPathWWTP)     # Add all nodes between wwtps to copy of sewers
                                                
                            # Initialisation wwwtp reconnection
                            nodes_BI = fastCopyNodes(nodes)
                            nodes_B3 = fastCopyNodes(nodes)
                            pumpWWTPListB1, pumpWWTPB3 = fastCopy(pumps), fastCopy(pumps)       # copy list with pumps
                            flowWWFrom = getFlowWWTP(WWTPs, WWTPFROM)                           # Get flow of largestWWTPID (same as idwWWTPZero)
                            flowWWTO = getFlowWWTP(WWTPs, WWTPTO)                               # Get flow of WWTPTO 
                                                               
                            notInaNetworkBeforeConnection = getNotInNetwork(pathBetweenWWTPs, sewers_NoCon)                                                                 # Get all nodes not in a network before interconnectin path is added
                            notInaNetwork = getPointsNotInNetwork(pathBetweenWWTPs, nodesFromNetwork, nodesToNetwork)                                                       # Get all nodes not in a network after innterconnection path is added
                            flowInitial_from, flowInitial_to = getFlowInitial(nodesFromNetwork, nodesToNetwork, nodes, flowWWFrom, flowWWTO, notInaNetworkBeforeConnection) # Calculate initial Flow:  if next node is in Network, substract flow from flowWWFrom 
                                    
                            # -----------
                            # Cost Module
                            # -----------

                            # --------
                            # Option 1
                            # --------
                            
                            # Sewer Costs
                            if len(nodesFromNetwork) > 1:
                                nodes_BI = getInflowingNodes(nodes_BI, nodesFromNetwork, flowInitial_from, notInaNetworkBeforeConnection)       # Only the inFlow from other nodes on the path gets calculated
                            else:
                                nodes_BI = assignInitialFlow(nodes_BI, pathBetweenWWTPs, flowInitial_from, -1)                                  # Change flow in nodes of first element
        
                            if len(nodesToNetwork) > 1:
                                nodes_BI = getInflowingNodes(nodes_BI, nodesToNetwork, flowInitial_to, notInaNetworkBeforeConnection)           # Only the inFlow from other nodes on the path gets calculated
                            else:
                                nodes_BI = assignInitialFlow(nodes_BI, pathBetweenWWTPs, flowInitial_to, 0)                                     # Change flow in nodes of last element

                            nodes_BI = clearFlow(nodes_BI, notInaNetwork)                                                                       # Clear flow where not in starting or ending network
                            nodes_BI = changeFlowAlongPath(nodes_BI, pathBetweenWWTPs)                                                          # Change Flow along the path in nodes_BI
                            nodes_BI, inflowNodesWTPB1 = correctTD(nodes_BI, archPathWWTPInvert[:-1], minTD, WWTPs, maxTD, minSlope, sewers)    # Set all trench depth except inflow nodes to minimum trench depth
                            nodes_BI, pumpWWTPListB1, edgeList1 = changeTD(nodes_BI, edgeList, pumpWWTPListB1, archPathWWTPInvert, maxTD, minSlope, inflowNodesWTPB1, sewers, minTD)
                            pipeCostsB1, edgeList1 = costsBetweenWWTPs(pathBetweenWWTPs, edgeList1, nodes_BI, 0, pumpWWTPListB1, minTD, discountYearsSewers, interestRate, stricklerC, operationCosts, f_SewerCost)  # first zero: regular flow, second zero: proportional costs
                                        
                            # Pumping Costs
                            pumpCostB1, _ = costPump(pathBetweenWWTPs, pumpWWTPListB1, pricekWh, pumpYears, interestRate)  # Calculate pump costs
                                        
                            # Option 1 - WWTPs costs
                            wtpCostB1 = costWWTP((flowWWFrom + flowWWTO), EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement)                 # WWTP-costs: Additional costs of building a larger WWTP
                                            
                            # --------
                            # Option 2
                            # --------
                                                        
                            # Sewer Costs       
                            if len(nodesFromNetwork) == 0 and len(nodesToNetwork) == 0:
                                pipeCostsB2 = 0                                                                                                                   # As there are none pipes on the way
                            else:
                                # The costs of the pipes on the way between the wwtps needs to be calculated
                                pipeCostB2a, _ = costsBetweenWWTPs(nodesFromNetwork, edgeList, nodes, 0, pumps, minTD, discountYearsSewers, interestRate, stricklerC, operationCosts, f_SewerCost)      # first zero: regular flow, second zero: proportional costs
                                pipeCostB2b, _ = costsBetweenWWTPs(nodesToNetwork, edgeList, nodes, 0, pumps, minTD, discountYearsSewers, interestRate, stricklerC, operationCosts, f_SewerCost)        # first zero: regular flow, second zero: proportional costs
                                pipeCostsB2 = pipeCostB2a + pipeCostB2b                                                                                                                                 # sum costs of the two networks
                                pumpCostB2A, _ = costPump(nodesFromNetwork, pumps, pricekWh, pumpYears, interestRate)                                         # [CHF] Calculate pump costs
                                pumpCostB2B, _ = costPump(nodesToNetwork, pumps, pricekWh, pumpYears, interestRate)                                           # [CHF] Calculate pump costs
                                pumpCostB2 = pumpCostB2A + pumpCostB2B
                                        
                            # WWTP costs   
                            wtpCostB2a = costWWTP(flowWWFrom, EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement)                            
                            wtpCostB2b = costWWTP(flowWWTO, EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement)
                                
                            if needsNetworkRemoving == 1:           # Calculate costs of crossed WWTPs
                                sumCostcrossedWWTP = getCostsOfCrossedWWTPs(allNodesToAddToPN, pathBetweenWWTPs, WWTPS_noCon, sewers_NoCon, nodes_noCon, EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement)
                                wtpCostB2 = wtpCostB2a + wtpCostB2b + sumCostcrossedWWTP     # Total wwtp costs
                            else:
                                wtpCostB2 = wtpCostB2a + wtpCostB2b     # Total wwtp costs
                                
                            # --------
                            # Option 3
                            # --------
                            
                            # Sewer costs
                            if len(nodesFromNetwork) > 1:  # Change nodes because of againstFlow (if needed) how it would look like with connection
                                nodes_B3 = getInflowingNodes(nodes_B3, nodesFromNetwork, flowInitial_from, notInaNetworkBeforeConnection)
                            else:
                                nodes_B3 = assignInitialFlow(nodes_B3, pathBetweenWWTPs, flowInitial_from, -1)  # Change initial flow in nodes   
                    
                            if len(nodesToNetwork) > 1:
                                nodes_B3 = getInflowingNodes(nodes_B3, nodesToNetwork, flowInitial_to, notInaNetworkBeforeConnection)
                            else:
                                nodes_B3 = assignInitialFlow(nodes_B3, pathBetweenWWTPs, flowInitial_to, 0)     # Change initial flow in nodes
                            
                            nodes_B3 = clearFlow(nodes_B3, notInaNetwork)                                       # Clear flow wheter there is no flow
                            nodes_B3 = changeFlowAlongPath(nodes_B3, pathBetweenWWTPsInvert)                    # Change flow along the path
                            nodes_B3, inflowNodesWTPB3 = correctTD(nodes_B3, pathBetweenWWTPsInvert[:-1], minTD, WWTPs, maxTD, minSlope, sewers)  # Set all trench depth except inflow nodes to minimum trench depth     
                            nodes_B3, pumpWWTPB3, edgeList3 = changeTD(nodes_B3, edgeList, pumpWWTPB3, archPathWWTP, maxTD, minSlope, inflowNodesWTPB3, sewers, minTD)
                                        
                            # Pumping Costs
                            pumpCostB3, _ = costPump(pathBetweenWWTPsInvert, pumpWWTPB3, pricekWh, pumpYears, interestRate)  # Calculate pump costs
                                        
                            # Sewer costs
                            pipeCostB3, edgeList3 = costsBetweenWWTPs(pathBetweenWWTPsInvert, edgeList3, nodes_B3, flowInitial_to, pumpWWTPB3, minTD, discountYearsSewers, interestRate, stricklerC, operationCosts, f_SewerCost)  # first zero: regular flow, second zero: proportional costs
                                    
                            # WWTP costs            
                            wtpCostA3 = wtpCostB1                             # Total wwtp costs
                            totCostBI = pipeCostsB1 + wtpCostB1 + pumpCostB1  # Total Costs Option 1        
                            totCostB2 = pipeCostsB2 + wtpCostB2 + pumpCostB2  # Total Costs Option 2
                            totCostB3 = pipeCostB3 + wtpCostA3 + pumpCostB3   # Total Costs Option 3
        
                            #arcpy.AddMessage("-------------Cost Comparison----------") 
                            #arcpy.AddMessage("totCostBI: " + str(totCostBI))
                            #arcpy.AddMessage("totCostB2: " + str(totCostB2))
                            #arcpy.AddMessage("totCostB2: " + str(totCostB3))
                            
                            wwtpSWAP, wwtpToWwtopConncetion = costComparison(totCostBI, totCostB2, totCostB3, 0, 0, 1, forceZ, finishedExpansion) # Compare costs, don't consider reasonable  costs
                                
                            if totCostBI - totCostB2 < totCostB3- totCostB2:
                                mergeCosts = totCostBI - totCostB2
                            else:
                                mergeCosts = totCostB3- totCostB2
                                
                            # Update if wwtp reconnection took place
                            if wwtpToWwtopConncetion == 1:
                                WWTPs, wwtpsIterate = checkIfWWTPwereConnected(pathBetweenWWTPs, wwtpsIterate, WWTPs)  # Check if by a wwtp connection another wwpt was connected as well If yes, delete these
                
                                # Changes depending on options
                                if wwtpSWAP == 1:
                                    wwtpsIterate = delEntry(wwtpsIterate, WWTPFROM)                            # Femove found wwtp from copylistWTP
                                    WWTPs = delEntry(WWTPs, WWTPFROM)                                          # Delete wwtp
                                    edgeList = fastCopy(edgeList3)                                             # Replace list with edges
                                    nodes = fastCopyNodes(nodes_B3)                                            # Replace nodes by new nodes with new flow
                                    WWTPs = updateFlowInWWTP(WWTPs, nodes, WWTPTO)                             # Update flow in WWTPs
                                    pumps = fastCopy(pumpWWTPB3)                                               # Replace pumps
                                    sewers = dict(sewers_B3)                                                   # Add new connection to sewers 
                                    expansionM = False
                                    sortedListWWTPs = delEntry(sortedListWWTPs, WWTPFROM)                      # Delete in wwtps to check for merging
                                    iterateMergeOptions = False                                                
                                        
                                    # If you want to merge until all but one WWTP are connected
                                    if ZtoReach > hypoZWeighted and finishedExpansion == 1:
                                        copyList = []
                                        newFlow = getFlowWWTP(WWTPs, WWTPTO)
                                        _, corrWWTPTO, _, _, _, _ = getPns(WWTPTO, nodes)
                                        for e in mergeCostStorage:
                                            if e[1][0] != WWTPFROM and e[3][0] != WWTPFROM:
                                                if e[1][0] == WWTPFROM or e[3][0] == WWTPFROM:
                                                    if e[1][0] == WWTPFROM:
                                                        eNew = [e[0], [WWTPTO, newFlow], WWTPTO, e[3]]
                                                        copyList.append(eNew)
                                                    else:
                                                        eNew = [e[0], e[1], e[2], [WWTPTO, corrWWTPTO]]
                                                        copyList.append(eNew)
                                                else:
                                                    copyList.append(e)
                                        mergeCostStorage = copyList

                                    # Calculate costs of current network, wwtps and if hypothetically all other wwtps wouldn't be connected but each have a decentral conncetion
                                    if iterativeCostCalc == 1 and ZtoReach > hypoZWeighted:
                                        hypoZ, hypoZWeighted, hypoCosts, hypoPumpCosts, hypoWWTPCosts, hypoPublicPipeCosts, _ = getFullHypotheticalCosts(aggregatetPoints, WWTPs, sewers, EW_Q, wwtpLifespan, interestRate, pumps, pumpYears, pricekWh, edgeList, nodes, stricklerC, discountYearsSewers, operationCosts, f_SewerCost, fc_wwtpOperation, fc_wwtpReplacement)
                                            
                                        if hypoZ > hypoZOld:
                                            totalSystemCosts.append([hypoZ, hypoZWeighted, hypoCosts, hypoPumpCosts, hypoWWTPCosts, hypoPublicPipeCosts ])
                                        hypoZOld = hypoZ
        
                                    # Create new PN and force connetion in nodes and restart the expansion-module if there was an intermediate network in between
                                    if needsNetworkRemoving == 1:
                                        PN, nodes, sortedListWWTPs = updatePrimEdges(allNodesToAddToPN, sewers, nodes, sortedListWWTPs)

                                        if PN == []: # ALL deleted nodes were on path, don't quit expnasion module
                                            break 
                                        else:
                                            nodes = removeForceCriteria(nodes, WWTPFROM)
                                            nodes = adForceCriteria(nodes, allNodesToAddToPN)                   # Add a force connection criteria that all deleted nodes get connected in the EM
                                            return nodes, sewers, pumps, WWTPs, PN, edgeList, totalSystemCosts  # Return to expansion-module
                                    break 
                                else:
                                    #arcpy.AddMessage("INFO: WWTP MERGE, NO SWAP. Delete WWTP: " + str(WWTPTO) + "  " + str(CurrentCheckedWWTP))  
                                    sortedListWWTPs = delEntry(sortedListWWTPs, WWTPTO)                         # Delete in wwtps to check for merging

                                    edgeList = fastCopy(edgeList1)                                              # Prim Edges

                                    wwtpsIterate = delEntry(wwtpsIterate, WWTPTO)                               # remove found wwtp from wwtpsIterate
                                    WWTPs = delEntry(WWTPs, WWTPTO)                                             # Delete wwtp
                                    nodes = fastCopyNodes(nodes_BI)                                             # Replace nodes by new nodes with new flow
                                    WWTPs = updateFlowInWWTP(WWTPs, nodes, WWTPFROM)                            # Update flow in WWTPs
                                    pumps = fastCopy(pumpWWTPListB1)                                            # Replace pumps list
                                    sewers_Current = appendToSewers(sewers_Current, pathBetweenWWTPs[1:-1])     # Add new connection to sewers sewers_Current                                                                # Used for next checking if connection is wortwhile
                                    sewers = dict(sewers_B1)                                                    # Add new connection to sewers    
                                    expansionM = False
                                    iterateMergeOptions = False                                                 
                                               
                                    # If Z = 1 wants to be reached:
                                    if ZtoReach > hypoZWeighted:
                                        copyList = []
                                        newFlow = getFlowWWTP(WWTPs, WWTPFROM)
                                        _, corWWTOFROM, _, _, _, _ = getPns(WWTPFROM, nodes)
                                        for e in mergeCostStorage:
                                            if e[1][0] != WWTPTO and e[3][0] != WWTPTO:
                                                # Replace mergecost elements
                                                if e[1][0] == WWTPTO:
                                                    eNew = [e[0], [WWTPFROM, newFlow], WWTPFROM, e[3]]
                                                    copyList.append(eNew)
                                                if e[3][0] == WWTPTO:
                                                    eNew = [e[0], e[1], e[2], [WWTPFROM, corWWTOFROM]]
                                                    copyList.append(eNew)
                                                else:
                                                    copyList.append(e)
                                        mergeCostStorage = copyList
                                                
                                    # Calculate costs of current network, wwtps and if hypothetically all other wwtps wouldn't be connected but each have a decentral conncetion
                                    if iterativeCostCalc == 1 and ZtoReach > hypoZWeighted:    
                                        hypoZ, hypoZWeighted, hypoCosts, hypoPumpCosts, hypoWWTPCosts, hypoPublicPipeCosts, _ = getFullHypotheticalCosts(aggregatetPoints, WWTPs, sewers, EW_Q, wwtpLifespan, interestRate, pumps, pumpYears, pricekWh, edgeList, nodes, stricklerC, discountYearsSewers, operationCosts, f_SewerCost, fc_wwtpOperation, fc_wwtpReplacement)
                                        if hypoZ > hypoZOld:
                                            totalSystemCosts.append([hypoZ, hypoZWeighted, hypoCosts, hypoPumpCosts, hypoWWTPCosts, hypoPublicPipeCosts])
                                        hypoZOld = hypoZ
                                        
                                    if needsNetworkRemoving == 1:
                                        PN, nodes, sortedListWWTPs = updatePrimEdges(allNodesToAddToPN, sewers, nodes, sortedListWWTPs)   
                                            
                                        if PN == []: # ALL deleted nodes were on path, therefore don't quit expansion module
                                            break 
                                        else:
                                            #arcpy.AddMessage("Return to Expansion Module...")
                                            nodes = removeForceCriteria(nodes, WWTPFROM)
                                            nodes = adForceCriteria(nodes, allNodesToAddToPN)                   # Add a force connection criteria that all deleted nodes get connected in the EM
                                            return nodes, sewers, pumps, WWTPs, PN, edgeList, totalSystemCosts  # Return to expansion-module

                                    # Create new PN and force connetion in nodes and restart the expansion-module
                                    if needsNetworkRemoving == 1:
                                        if PN == [] and OnlyExecuteMerge != 1:                                  # ALL deleted nodes were on path, don't quit expnasion module
                                            break 
                                        else:
                                            nodes = removeForceCriteria(nodes, WWTPTO)                          # Remove force conneciton criteria
                                            return nodes, sewers, pumps, WWTPs, PN, edgeList, totalSystemCosts  # Return to expansion-module
                                    break                                                                       # Stop iteration merge options
                            else:                                                                               # If no merging took place Stop iteration
                                #arcpy.AddMessage("Do not make any changes...")
                                wwtpsIterate = delEntry(wwtpsIterate, WWTPFROM)                                         # delete wwtps in iteration list of wwtps
                                sortedListWWTPs = fastCopy(sortedListWWTPs_noCon)                                       # Restore as it was before because not connection took place
                                sewers = dict(sewers_NoCon)                                                             # Restore as it was before because not connection took place
                                pumps = fastCopy(pumps_noCon)                                                           # Restore as it was before because not connection took place
                                nodes = fastCopyNodes(nodes_noCon)                                                      # Restore as it was before because not connection took place
                                WWTPs = fastCopy(WWTPS_noCon)                                                           # Restore as it was before because not connection took place
                                PN = fastCopy(PN_noCon)                                                                 # Restore as it was before because not connection took place
                                expansionM = False                                                                      # If no reconnection, stop checking
                                iterateMergeOptions = False                                                             # stop iterating options
                                mergeCostStorage.append([mergeCosts, CurrentCheckedWWTP, checkBackConnectionID, mOpt])  # In case no merge store costs of additional merging. This is used in case we want to force higher Zw values

                    # If merge took place and no swap took place, add the WWTP to sortedListWWTPs to check again. 
                    if wwtpToWwtopConncetion == 1:
                        if wwtpSWAP == 1:
                            #In case the WWTPTO is not in sortedListWWTPs, insert at top position
                            addWWTPagain = 1
                            for i in sortedListWWTPs:
                                if i[0] == WWTPTO: 
                                    addWWTPagain = 0
                                    break
                            if addWWTPagain == 1:
                                for wwtp in WWTPs:
                                    if wwtp[0] == WWTPTO:
                                        sortedListWWTPs.insert(0, [wwtp[0], wwtp[1]])
                    else:                                                           # If all three options were not merged
                        del sortedListWWTPs[0]                                      # Remove wwtp to check     
                break                              

    return nodes, sewers, pumps, WWTPs, PN, edgeList, totalSystemCosts

def getCheapestMerge(mergeCostStorage, WWTPs):
    '''
    This function gets from all merges the cheapest. This function is used in case
    SNIP is artifically continued to reach higher Z values.
    
    Input:
    mergeCostStorage    -    List storing costs of all merges
    WWTPs               -    List with WWTPs
    
    Outpt:
    mergeCostStorage    -    All merges withhout cheapest merge
    '''
    cheapestMerge = 999999999999
    
    # Iterate all merges and get cheapest
    for e in mergeCostStorage:
        if e[0] < cheapestMerge:
            ok, ok2 = checkifInWWTPs(WWTPs, e[3][0]), checkifInWWTPs(WWTPs, e[2])
            if ok == True and ok2 == True and e[3][0] != e[2]:
                cheapestM = e
    # Remove element
    copyList = []
    for e in mergeCostStorage:
        ok, ok2 = checkifInWWTPs(WWTPs, e[3][0]), checkifInWWTPs(WWTPs, e[2])
        if ok == True and ok2 == True and e[3][0] != e[2]:    
            if e[0] != cheapestMerge: 
                copyList.append(e) 
    mergeCostStorage = copyList
    return mergeCostStorage, cheapestM

def updatePrimEdges(allNodesToAddToPN, network, nodes, sortedListWWTPs):
    ''' 
    Reconnect intermediate network nodes. Change list with prim distances (PN) in order that expansion module can get started again
    
    Input:
    allNodesToAddToPN  -    List with all nodes to add to PN
    network            -    Sewer Network
    nodes              -    All nodes
    sortedListWWTPs    -    List with wwpt to check for a merge
    
    Output:
    PN_new             -    Updatet Prim Edges
    nodes              -    Updated nodes
    sortedListWWTPs    -    updated list to check for a merge
    '''
    PN_new, delNodes, nodesUpdateCalc = [], [], []  # List to store from each node the shortest edge to the network
                          
    # Get all nodes of network inbetween which need to be deleted
    for i in allNodesToAddToPN:
        for e in i[1]:
            delNodes.append(e) # All Nodes which need recalculation

    # Correct flow and only read out populated nodes
    for i in delNodes:
        for e in nodes:
            if e[0] == i:
                if e[8] != 0:
                    nodesUpdateCalc.append(i)
                break
    
    # Iterate all not connected populated nodes and get shortest edge to Network
    for nodeToCon in nodesUpdateCalc:
        intermediatePN = []
        idp1, p1, isPop, _, _, _ = getPns(nodeToCon, nodes)  # Get unconnected node

        if isPop > 0: 
            if idp1 in network:

                # Remove nodes already from list with wwtps to check for a merge.
                cntt, RemoveFromMerge = 0, False
                for e in sortedListWWTPs:
                    if e[0] == idp1:
                        delPosition, RemoveFromMerge = cntt, True
                        break
                    cntt += 1
                if RemoveFromMerge == True:
                    del sortedListWWTPs[delPosition]
            else:              
                # Calculate distance from origin to all nodes
                for i in nodes:
                    if i[0] not in nodesUpdateCalc and i[0] in network: 
                        idp0, p0, flowCriteria = i[0], (i[1], i[2], i[3]), i[8]  
        
                        # if not point to itself, point not already connected and no ArchPoint
                        if flowCriteria > 0:                                                                    
                            distanz3d, _, _ = distanceCalc3d(p1, p0)  
                            intermediatePN.append([distanz3d, idp0, idp1, p1[0], p1[1], p1[2], 1])
            
                # Get closest PN to network
                closestDist = 99999999999
                for f in intermediatePN:
                    if f[0] < closestDist:
                        closestDist, closestEdge = f[0], f
                PN_new.append(closestEdge)
    
    # Change criteria in order that connection is forced
    for i in PN_new:
        for e in nodes:
            if e[0] == i[1]:
                e[5] = 1    # Make that connection is forced
                break
    return PN_new, nodes, sortedListWWTPs

def storeMergingOptions(nodeIdOpt1, nodeOpt1Cor, nodeIdOpt2, nodeOpt2Cor, nodeIdOpt3, nodeOpt3Cor):
    ''' 
    This functions stores different merging options in a list.
    
    Input:
    nodeIdOpt1             -    ID option 1
    nodeOpt1Cor            -    coordinate option 1
    nodeIdOpt2             -    ID option 2
    nodeOpt2Cor            -    coordinate option 2
    nodeIdOpt3             -    ID option 3
    nodeOpt4Cor            -    coordinate option 4

    Output:
    mergeOptions           -    List with merging options
    '''
    mergeOptions = []
    if nodeIdOpt1 != None:
        mergeOptions = [[nodeIdOpt1, nodeOpt1Cor]]
                        
    if nodeIdOpt2 != None:
        insert = False
        for i in mergeOptions:
            if i[0] != nodeIdOpt2:
                insert = True
                break
        if insert == True or len(mergeOptions) == 0:
            mergeOptions.insert(0, ([nodeIdOpt2, nodeOpt2Cor]))

    if nodeIdOpt3 != None:
        insert = False
        for i in mergeOptions:
            if i[0] != nodeIdOpt3:
                insert = True
                break
        if insert == True or len(mergeOptions) == 0:
            mergeOptions.insert(0, ([nodeIdOpt3, nodeOpt3Cor]))                  
    return mergeOptions

def relax(StreetNetwork, u, v, D, P):
    """
    Relaxing function of Djikstra Algorithm

    Input Arguments: 
    StreetNetwork    --    StreetNetwork
    u,v,D            --    Djikstra-related varaibles
    P                --    Djikstra-List with all distances to each node
    """
    inf = float('inf')
    d = D.get(u, inf) + StreetNetwork[u][v] # Possible shortcut estimate
    if d < D.get(v, inf):                   # Is it really a shortcut?
        D[v], P[v] = d, u                   # Update estimate and p

def dijkstraAlgorithm(StreetNetwork, s):
    """
    Djikstra Algorithm to find shortest rout on street network.

    Input Arguments: 
    StreetNetwork    --    StreetNetwork
    s                --    Starting node

    Output Arguments:
    D                --    Distances to every node from s
    P                --    Djikstra-List
    """
    D, P, Q, S = {s:0}, {}, [(0, s)], set()  

    while Q:  # All nodes are checked?
        smallesvalue, cnt = 99999999, 0
        for i in Q:
            if i[0] < smallesvalue:
                smallesvalue, u = i[0], i[1]
                smallesPosition = cnt
            cnt += 1
        del Q[smallesPosition]

        if u in S: continue                     # Already visited? Skip it
        S.add(u)                                # We've visited it now
        for v in StreetNetwork[u]:              # Go through all its neighbors
            relax(StreetNetwork, u, v, D, P)    # Relax the out-edges
            Q.append((D[v], v))                 # visited
    return D, P  

def getClosestNode(PN):
    """
    This function sorts a list with all nodes left to connect and deletes the closest node. 

    Input Arguments: 
    PN                    --    Distances to all nodes
    
    Output Arguments:
    PN                    --    Updates distances to all nodes
    minDit                --    Smallest distance
    fromNode              --    Node FROM
    toNode                --    Node TO
    factorDistanz         --    Distance factor
    euclidianDistance     --    Euclidian distance
    """
    minDit, zahler = 9999999999, 0 
    for i in PN:
        if i[0] < minDit:
            minDit, fromNode, toNode, factorDistanz = i[0], i[1], i[2], i[6]   
            deletPosition = zahler
        zahler += 1
    del PN[deletPosition]

    # calculate real Distance
    if minDit > 0:
        euclidianDistance = float(minDit) / float(factorDistanz)                # Euclidian distance converted back because was weighted
    else:
        euclidianDistance = 0
    return PN, fromNode, toNode, factorDistanz, euclidianDistance

def dijkstra(streetNetwork, idp0, idp1, heightDiff):
    """
    This function gets the path from the Djikstra list.

    Input Arguments: 
    streetNetwork             --    Distances to all nodes
    idp0                      --    Start node
    idp1                      --    Endnode
    heightDiff                --    Height Difference
    
    Output Arguments:
    archPathList              --    Updates distances to all nodes.
    distStartEnd              --    Distance between the two nodes.
    slopeDijkstra             --    Slope between the two nodes on Djikstra distance.
    """                                              
    distances, listDijkstra = dijkstraAlgorithm(streetNetwork, idp0)                # calculate djikstra distances
    scrapPathDjika, distStartEnd = writePath(listDijkstra, distances, idp0, idp1)                 
    archPathList = archPath(scrapPathDjika, distances)                              # the achPathList contains all intermediary, not pouplated nodes of the street network. List gets afterwards appended to P
    slopeDijkstra = (float(heightDiff) / float(distStartEnd)) * 100                 # Slope in % of the djikstra-street distance                    
    return archPathList, distStartEnd, slopeDijkstra

def createStreetGraph(in_FC):
    """
    This function feads out the coordinates from a line shapefile

    Input Arguments: 
    in_FC    --    Path to Shapefile

    Output Arguments:
    edges    --    list with all edges edges = ([(X_FROM, Y_FROM), (X_TO, Y_TO)], ...)
    """
    rows = arcpy.SearchCursor(in_FC)  # the rows are read in the shapefile
    fields = arcpy.ListFields(in_FC)  # the fields are read in the shapefile
    edges = []
    
    for row in rows:                                # Iterate row
        for field in fields:                        # Iterate fields
            if field.name == "X_START":
                X_FROM = float(row.getValue(field.name))
            if field.name == "X_END":
                X_TO = float(row.getValue(field.name))
            if field.name == "Y_START":
                Y_FROM = float(row.getValue(field.name))
            if field.name == "Y_END":
                Y_TO = float(row.getValue(field.name))
        z = [(X_FROM, Y_FROM), (X_TO, Y_TO)]  
        edges.append(z)
    return edges

def createPolyPoint(pntsFlow, outListStep_point):
    """
    This function creates polypoints for ArcGIS.

    Input Arguments: 
    in_FC    --    Path to Shapefile

    Output Arguments:
    edges    --    list with all edges edges = ([(X_FROM, Y_FROM), (X_TO, Y_TO)], ...)
    """
    punkte, pointList = arcpy.Point(), []
    for entry in pntsFlow:
        punkte.X, punkte.Y = entry[1], entry[2]
        pointList.append(arcpy.PointGeometry(punkte))

    if len(pointList) == 0: # Could not create PolyPoint as no points to draw are found
        draw = 0
        return draw
    
    arcpy.CopyFeatures_management(pointList, outListStep_point)
    draw = 1
    return draw

def createPolyPointWWTP(pntsFlow, outListStep_point):
    """
    This function creates polypoints for ArcGIS.

    Input Arguments: 
    in_FC    --    Path to Shapefile

    Output Arguments:
    edges    --    list with all edges edges = ([(X_FROM, Y_FROM), (X_TO, Y_TO)], ...)
    """
    punkte, pointList = arcpy.Point(), []
    for entry in pntsFlow:
        punkte.X, punkte.Y = entry[2], entry[3]
        pointList.append(arcpy.PointGeometry(punkte))

    if len(pointList) == 0:
        arcpy.AddMessage("Could not create PolyPoint as no points to draw are found")
        draw = 0
        return draw
    
    arcpy.CopyFeatures_management(pointList, outListStep_point)
    draw = 1
    return draw

def createPolyPointPump(pntsFlow, outListStep_point):
    """
    This function create shapefiles for ArcGIS.
    
    Input Arguments: 
    pntsFlow             --    ist with point coordinates
    outListStep_point    --    Path to Output Folder
    
    Output Arguments:
    draw                 --    Criteria wheter points are drawn or not (0: no, 1: yes)
    """
    punkte, pointList, draw = arcpy.Point(), [], 1

    if len(pntsFlow) == 0:
        arcpy.AddMessage("INFO: No pumps are drawn ")
    
    for entry in pntsFlow:
        punkte.X, punkte.Y = entry[1], entry[2]
        pointList.append(arcpy.PointGeometry(punkte))

    if len(pointList) == 0: # Could not create PolyPoint as no points to draw are found
        draw = 0
        return draw
    arcpy.CopyFeatures_management(pointList, outListStep_point)
    return draw

def writeWWTPs(shapefile_path, pointListNear):
    """
    This function create shapefiles for ArcGIS.
    
    Input Arguments: 
    shapefile_path       --    Path to Shapefile
    pointListNear        --    List with fields
    
    """
    from arcpy import env
    env.workspace = shapefile_path
    arcpy.AddField_management(shapefile_path, "ID_scratch", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "Id")              # Delete automatically created field
    arcpy.AddField_management(shapefile_path, "ID", "FLOAT")
    arcpy.AddField_management(shapefile_path, "FLOW", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "ID_scratch")      # Delete created field
    
    rows = arcpy.UpdateCursor(shapefile_path)
    cnt = 0
      
    for row in rows:
        row.setValue("ID", pointListNear[cnt][0])
        row.setValue("FLOW", pointListNear[cnt][1])
        rows.updateRow(row)
        cnt += 1
    return

def writeStartnode(shapefile_path, startnodeID):
    """
    This function create shapefiles for ArcGIS.
    
    Input Arguments: 
    shapefile_path       --    Path to Shapefile
    pointListNear        --    List with fields
    """
    from arcpy import env
    env.workspace = shapefile_path
    arcpy.AddField_management(shapefile_path, "ID_scratch", "FLOAT")
      
    arcpy.DeleteField_management(shapefile_path, "Id")  # Delete automatically created field
    arcpy.AddField_management(shapefile_path, "ID", "FLOAT")
    arcpy.AddField_management(shapefile_path, "FLOW", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "ID_scratch")  # Delete created field
    
    rows = arcpy.UpdateCursor(shapefile_path)
      
    for row in rows:
        row.setValue("ID", startnodeID)
        rows.updateRow(row)
    return

def readOutAllStreetVerticesAfterAggregation(nodes, rasterPoints, rasterSize):
    """
    This function reads out the aggregated street inlets
    
    Input Arguments: 
    nodes               --    nodes
    rasterPoints        --    Raster Points
    rasterSize          --    Raster Size
    
    Output Arguments:
    streetVert          --    Vertices of Street Network
    """
    rows = arcpy.SearchCursor(nodes)  # define rows
    fields = arcpy.ListFields(nodes)  # define fields
    streetVert, IDNEW = [], 100000

    # Iterate shapefile 
    for row in rows:
        for field in fields:
            if field.name == "X_START":
                X_start = row.getValue(field.name)
            if field.name == "X_END":
                X_end = row.getValue(field.name)
            if field.name == "Y_START":
                Y_start = row.getValue(field.name)
            if field.name == "Y_END":
                Y_end = row.getValue(field.name)

        # Assign height from closest rasterPoint (Takes a lot of time)
        fertig1, fertig2 = 0, 0

        for i in rasterPoints: 
            if (i[1] <= X_start + rasterSize and i[1] >= X_start - rasterSize) and (i[2] <= Y_start + rasterSize and i[2] >= Y_start - rasterSize):
                heightSTART, fertig1 = i[3], 1
            if i[1] <= X_end + rasterSize and i[1] >= X_end - rasterSize and i[2] <= Y_end + rasterSize and i[2] >= Y_end - rasterSize:
                heightEND, fertig2 = i[3], 1
            if fertig1 == 1 and fertig2 == 1:
                #arcpy.AddMessage("Height was found of start and end node of street arc")
                break
        
        if fertig1 == 0 and fertig2 == 0:
            arcpy.AddMessage(" ERROR: Node of streetwork could not get assigned heigh of raster cell as no raster call was found." + str(fertig1) + " " + str(fertig2))
            #continue # In case a street has no height
               
        copyZ, copyZ2 = 1, 1

        # Prevent that dublicate points are added
        for pnt in streetVert:
            if pnt[1] == X_start and pnt[2] == Y_start:
                copyZ = 0
            if pnt[1] == X_end and pnt[2] == Y_end:
                copyZ2 = 0

        # If due to rounding the start and end becomes the same
        if X_start == X_end and Y_start == Y_end:
            copyZ2 = 0

        if copyZ == 1:
            IDNEW += 1
            newZ = [IDNEW, X_start, Y_start, heightSTART]
            streetVert.append(newZ)
                
        if copyZ2 == 1:
            IDNEW += 1
            newz2 = [IDNEW, X_end, Y_end, heightEND]
            streetVert.append(newz2)
    return streetVert

def mergePathClosestNetwork(pathInClosestNetwork, archPathWWTP, edgeList):
    """
    This function appends to the path from a WWTO to the node of a closest network the path in the closest network to the
    wwtp in the closest network. The edges are iterated to get the length.
    
    In case there are "triangle-Distances" (distances which are walked twice) they are removed.
    
    Input Arguments: 
    pathInClosestNetwork    --    Path in closeset network to wwtp
    archPathWWTP            --    Path from from wwtp to node of closest network
    
    Output Arguments:
    archPathWWTP            --    Updated path
    """
    cnt = 0
    for e in pathInClosestNetwork:
        currentEntry = e
        if cnt == 1:
            for s in edgeList:
                if s[0][0] == currentEntry and s[1][0] == entryBefore or s[0][0] == entryBefore and s[1][0] == currentEntry:
                    archPathWWTP.insert(0, [currentEntry, [entryBefore, s[2]]])
                    break
            entryBefore = e
        if cnt == 0:
            entryBefore, cnt = e, 1
    
    # Find Triangle- Double Path (path which is double because of djikstra-street search)
    noTraingleArchPath, TriangleNetworktoRemove, shortCutNode = [], [], None      
    cnt = 1
    for e in archPathWWTP:
        nodeToTest, t = e[0], []
        
        # Iterate rest of archPath to check if double entry
        for a in archPathWWTP[cnt:]:
            if a[0] != nodeToTest:
                t.append(a[0])
            else:
                shortCutNode = a[0]
                TriangleNetworktoRemove = t
                break
        
        if TriangleNetworktoRemove != []:   
            break           # abort    
        cnt += 1    

    firstShortcutNotAdded = 0
    
    if shortCutNode == None:
        return archPathWWTP
    else:
        for e in archPathWWTP:
            if e[0] not in TriangleNetworktoRemove:
                if firstShortcutNotAdded == 0:
                    if e[0] == shortCutNode:
                        firstShortcutNotAdded = 1
                    else:
                        noTraingleArchPath.append(e)  
                else:
                    noTraingleArchPath.append(e)
            else:           # Remove node from TriangleNetworktoRemove
                cnt = 0
                for i in TriangleNetworktoRemove:
                    if i == e[0]:
                        del TriangleNetworktoRemove[cnt]
                        break
                    cnt += 1
        return noTraingleArchPath

def delnodesNotUsedInSewer(sewers, pathtonearestWTP, archPath, sewerBefore):
    '''
    This Function check if there were added nodes from the archPath to the sewers and then not used. In case there
    are such points, they get removed.
    
    Input:
    pathtonearestWTP    -    Path between nodes
    archPath            -    Path added in iteration
    sewerBefore         -    State of sewer before iteration

    Output:
    sewers              -    Updated sewers
    '''                
    nodesInSewerNotUsed = []
    for e in archPath[1:]:
        isInPath = False
        for i in pathtonearestWTP:
            if e[0] == i:
                isInPath = True
                break
        
        if isInPath == False:
            try:    
                _ = sewerBefore[e[0]]  # Test if not in sewer
            except:
                nodesInSewerNotUsed.append(e[0])
   
    if nodesInSewerNotUsed != []:
        sewers = delFromSewers(sewers, nodesInSewerNotUsed)    # Delete all nodes not used in sewers    
    return sewers

def readStreetPoints(nodes):
    """
    This function reads in the street vertices (which need to be optained from line-shape file)
    
    Input Arguments: 
    nodes             --   Point file with Fields "X","Y", "StreetID" 

    Output Arguments:
    streetVert        --    Nodes of street network
    """
    rows = arcpy.SearchCursor(nodes)  # Rows of the shapefile
    fields = arcpy.ListFields(nodes)  # Fields of the shapefile
    streetVert = []

    # Iterate shapefile 
    for row in rows:
        for field in fields:
            if field.name == "StreetID":
                StreetID = int(row.getValue(field.name))
            if field.name == "X":
                X = float(row.getValue(field.name))
            if field.name == "Y":
                Y = float(row.getValue(field.name))
        z = (StreetID, X, Y)
        streetVert.append(z)  
    return streetVert

def addedgesID(edges, streetVertices):
    """
    This function assing Ids and dicstance to edges
    
    Input Arguments: 
    edges                 --   list with coordinates of edges
    streetVertices        --   list with coordinates of streetVertices
    
    Output Arguments:
    newEdges              --   List with edges and distance: (((id, p1x, p1y), (id, p1x, p1y), distance, slope), .....)
    """
    pipeDiameter =  0.0 # Initially it is zero
    newEdges = []
    for edge in edges:
        pnt1, pnt2 = edge[0], edge[1]
        p1x, p1y, p2x, p2y = pnt1[0], pnt1[1], pnt2[0], pnt2[1]

        for i in streetVertices:
            if i[1] == p1x and i[2] == p1y:
                p0 = [i[0], p1x, p1y, i[3]]
            if i[1] == p2x and i[2] == p2y:
                p1 = [i[0], p2x, p2y, i[3]]
                
        distanz, slope, _ = distanceCalc3d((p0[1], p0[2], p0[3]), (p1[1], p1[2], p1[3]))
        
        if distanz != 0:
            wayProperty = 1 
            newEdges.append([p0, p1, distanz, slope, pipeDiameter, wayProperty])
    return newEdges

def appendStreetIDandCreateGraph(edges):
    """
    This function appens the ID to the read out coordinates of the line-file. Writes it in graph form
    
    Input Arguments: 
    edges             --    All edges
    
    Output Arguments:
    graph             --    Street Network Graph
    """
    graph = {}  # empty new graph

    for entry in edges:
        to_pt, from_pt, distanz = entry[0], entry[1], entry[2]
        ID_from, ID_to = from_pt[0], to_pt[0]
        
        # Allow that more than one entry is in dictionary
        if ID_from not in graph:
            graph[ID_from] = {ID_to: distanz}

        if ID_from in graph:
            z = graph[ID_from]
            z[ID_to] = distanz
            graph[ID_from] = z
                
        if ID_to not in graph:
            graph[ID_to] = {ID_from: distanz}
        
        if ID_to in graph:
            z = graph[ID_to]
            z[ID_from] = distanz
            graph[ID_to] = z
    return graph

def writePath(dijkstraList, distancesDijka, start, end):
    """
    This function writes out the direct path from list generated with djikstra algorithm. 
    Start node needs to be the same from where the dijkstra list is calculcated.
    
    Input Arguments: 
    dijkstraList            --    All edges
    distancesDijka          --    List with dijkstra distances 
    start                   --    start node
    end                     --    end node
    
    Output Arguments:
    path                    --    Path
    distStartEnd            --    Total distance from start to end
    """
    path, up, distStartEnd = [end], None, distancesDijka[end]
    
    while up != start and start != end: # As long the startin node is not found
        up = dijkstraList[end]
        path.append(up)
        end = up
    return path, distStartEnd

def readBuildingPoints(in_FC):
    """
    This function reads out the buildling points from a shapefile.
    
    Input Arguments: 
    in_FC             --   Path to point shapefile
    
    Output Arguments:
    pointListNear     --    List with buildlings
    """
    pointListNear, rows, fields = [], arcpy.SearchCursor(in_FC), arcpy.ListFields(in_FC) 
   
    for row in rows:  # Iterate shapefile
        for field in fields:
            if field.name == "FID":
                ID = row.getValue(field.name)
            if field.name == "POINT_X":
                X = row.getValue(field.name)  
            if field.name == "POINT_Y":
                Y = row.getValue(field.name)
            if field.name == "Q":  # flow amount
                quantity = row.getValue(field.name)
        pointListNear.append((ID, X, Y, 0, quantity))  # ID, x, y, pop_orig
    return pointListNear

def readClosestPointsAggregate(in_FC):
    """
    This function reads out the closest Points for each building after a near calculation was performed
    
    Input Arguments: 
    in_FC             --   Path to point shapefile.
    
    Output Arguments:
    pointListNear     --     Form: [[ID, x, y, z, flowInNode (is zero because no connetion), flow in node, joker position, flow in Node], ...]
    """
    rows = arcpy.SearchCursor(in_FC) 
    fields = arcpy.ListFields(in_FC) 
    pointListNear = []
    
    # Iterate shapefile
    for row in rows:
        for field in fields:
            if field.name == "FID": 
                ID = row.getValue(field.name)
            if field.name == "NEAR_X":
                NEAR_X = row.getValue(field.name)
            if field.name == "NEAR_Y":
                NEAR_Y = row.getValue(field.name)
            if field.name == "Q":  
                Q = row.getValue(field.name)
                
            # NEW 12.03.2015
            if field.name == "POINT_X":
                X = row.getValue(field.name)  
            if field.name == "POINT_Y":
                Y = row.getValue(field.name)  
                
        # Check if Near was not calculated. If Near was too far away because of distancecriteria, ArcGIS per default returns the value -1
        if NEAR_X == -1 or NEAR_Y == 1:
            z = [ID, X, Y, 0, 0, 0, 0, 0, Q]
        else:
            z = [ID, NEAR_X, NEAR_Y, 0, 0, 0, 0, 0, Q]
        pointListNear.append(z)
    return pointListNear

def aggregate(pointListNear, AggregateKritStreet, outListStep_point, minTD):
    """
    This function aggregates the buildings on the streets according to given input parameter AggregateKritStreet (sewer inlets).
    If a sewer inlet is within a defined distance to a next sewer inlet, the sewer inlets are merged and aggregated.
    
    Input Arguments: 
    pointListNear             -   Nearest Points to buildlings on network. Aggregated.
    AggregateKritStreet       -   Criteria to aggregate streets.
    outListStep_point         -   Path to point shapefile. 
    minTD                     -   Minimum Trench Depth
    
    ID, x, y, pop, flow, nearDist, costExist, costDezentral, height, Quan (This is the amount of Watewater from each Building
    Output Arguments:
    aggregatetStreetInlets    -    Sewer inlets
    buildingList              -    list with buildings
    """
    pointList, buildingList, aggregatetStreetInlets, punkte = [], [], [], arcpy.Point()
    
    # Aggregate points. The criteria AggregateKritStreet decides wheter a new point to the street is added or a aggregated
    for i in pointListNear:
        IdGeb, p0, StartDist = i[0], (i[1], i[2], i[3]), 10000000000

        
        # Calculate wheter there is a point closer than the defined distance
        if len(aggregatetStreetInlets) > 0:     
            for e in aggregatetStreetInlets:
                p1 = (e[1], e[2], e[3])
                distanz = distanceCalc2d(p0, p1)
    
                # Search closest closestStreetInlet to get distance
                if distanz < StartDist:  
                    closestX, closestY, closest_height, closest_quan, StartDist = e[1], e[2], e[3], e[8], distanz

            # If cloest point is further away than closest StreetInlet, make new StreetInlet
            if StartDist > AggregateKritStreet:
                gebListe = [IdGeb]  # Add geb to gebListe
                aggregatetStreetInlets.append([i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], gebListe, i[3] - minTD])
                buildingList.append((i[1], i[2], [IdGeb]))  # add information which buildling waters to which streetinlet
    
            # if closest StreetInlet is closer than AggregateKritStreet, merge StreetInlet
            if StartDist < AggregateKritStreet:      
                for e in aggregatetStreetInlets:
                    if e[1] == closestX and e[2] == closestY:
                        
                        # Add BuildingID
                        for building in buildingList:
                            if closestX == building[0] and closestY == building[1]:
                                building[2].append(IdGeb)
        
                        e[9].append(IdGeb)  # Add Buildings

                        # Sum quan
                        existingQuan = i[8]
                        newQuan = (existingQuan + closest_quan)
                        e[1], e[2], e[8] = closestX, closestY, newQuan  # Change nearest point in building
        
                        # Sum the heigths
                        existingHeight = i[3]
                        summedHeight = (existingHeight + closest_height)
                        e[3] = summedHeight
                        break
                    
        # add initial point                
        if len(aggregatetStreetInlets) == 0:
            buildingList.append((i[1], i[2], [IdGeb]))
            aggregatetStreetInlets.append([i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], [IdGeb], i[3] - minTD])

    # Assign unique ID to aggregated nodes
    aggregID = 0
    for i in aggregatetStreetInlets:
        i[0] = aggregID
        aggregID += 1

    # Divide the summed height by the number of aggregated nodes
    for i in aggregatetStreetInlets:
        summe, anzBuildinging = i[3], len(i[9])
        i[3] = float(summe) / float(anzBuildinging)

    # Draw nodes in ArcGIS
    for entry in aggregatetStreetInlets:
        punkte.X, punkte.Y = entry[1], entry[2]
        pointList.append(arcpy.PointGeometry(punkte))
    
    arcpy.CopyFeatures_management(pointList, outListStep_point)
    return aggregatetStreetInlets, buildingList

def removeForceCriteria(nodes, ID):
    '''
    This function removes the force criteria.
    Input:
    nodes    -    All nodes
    ID       -    ID to remove force criteria
    
    Output:
    nodes    -    All nodes with updated force criteria
    '''
    for i in nodes:
        if i[0] == ID:
            i[5] = 0
            break
    return nodes

def adForceCriteria(nodes, allNodesToAddToPN):
    '''
    This function adds the force criteria.
    Input:
    nodes      -    All nodes
    allNodesToAddToPN   -    nodeList with IDs to remove force criteria
    
    Output:
    nodes    -    All nodes with updated force criteria
    
    '''
    for f in allNodesToAddToPN:
        for z in f[1]:
            for i in nodes:
                if z == i[0]:
                    i[5] = 1        # Make that a connection is enforced (later in the EM)
                    break
    return nodes

def drawAllNodes(nodes, allNodesPath):
    """
    This function create points of all nodes in a list
    
    Input Arguments: 
    nodes              --    list with nodes
    allNodesPath       --    path
    """
    punkte, pointList = arcpy.Point(), []

    # Draw Points
    for entry in nodes:
        punkte.X, punkte.Y = entry[1], entry[2]
        pointList.append(arcpy.PointGeometry(punkte))
    arcpy.CopyFeatures_management(pointList, allNodesPath)  # Create a copy of the Polyline objects
    return

def splitStreetwithInlets(in_FC, outListStep_point, aggregatetStreetFile):
    """
    This function splits a ArcGIS line with given input Points (in this case the street inlets)     
    
    Input Arguments: 
    in_FC                   --   Path to shapefile to split
    outListStep_point       --   Shapefile Path
    aggregatetStreetFile    --   Shapefile Path
    """ 
    arcpy.SplitLineAtPoint_management(in_FC, outListStep_point, aggregatetStreetFile, "2 Meters")  # The Search radius would be two Meters.
    return

def updatefieldsPoints(shapefile_path):
    """
    This function creates fields in a shapefile.
    
    Input Arguments: 
    shapefile_path    --   Path to point shapefile.
    """
    from arcpy import env
    env.workspace = shapefile_path
    arcpy.AddField_management(shapefile_path, "ID_scratch", "FLOAT")
      
    arcpy.DeleteField_management(shapefile_path, "Id")              
    arcpy.AddField_management(shapefile_path, "ID", "FLOAT")
    arcpy.AddField_management(shapefile_path, "StreetID", "FLOAT")
    arcpy.AddField_management(shapefile_path, "NEAR_DIST", "FLOAT")
    arcpy.AddField_management(shapefile_path, "LENGTH", "DOUBLE")
    arcpy.AddField_management(shapefile_path, "X_START", "DOUBLE")
    arcpy.AddField_management(shapefile_path, "X_END", "DOUBLE")
    arcpy.AddField_management(shapefile_path, "Y_START", "DOUBLE")
    arcpy.AddField_management(shapefile_path, "Y_END", "DOUBLE")
    arcpy.DeleteField_management(shapefile_path, "ID_scratch")                                      # Delete automatically created field 
    arcpy.CalculateField_management(shapefile_path, "LENGTH", "float(!SHAPE.LENGTH!)", "PYTHON")    # Update field in ArcGIS
    
    xStart = 'float(!shape.firstpoint!.split() [0])'
    xEnd = 'float(!shape.lastpoint!.split() [0])'
    yStart = 'float(!shape.firstpoint!.split() [1])'
    yEnd = 'float(!shape.lastpoint!.split() [1])'
      
    arcpy.CalculateField_management(shapefile_path, "X_START", xStart, "PYTHON")
    arcpy.CalculateField_management(shapefile_path, "X_END", xEnd, "PYTHON")
    arcpy.CalculateField_management(shapefile_path, "Y_START", yStart, "PYTHON")
    arcpy.CalculateField_management(shapefile_path, "Y_END", yEnd, "PYTHON")
    return

def updateFieldNode(shapefile_path):
    """
    This function creates fields in a shapefile.
    
    Input Arguments: 
    shapefile_path             --   Path to point shapefile.
    """
    from arcpy import env
    env.workspace = shapefile_path
    arcpy.AddField_management(shapefile_path, "ID_scratch", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "Id")          # Delete automatically created field
    arcpy.AddField_management(shapefile_path, "ID", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "ID_scratch")  # Delete automatically created field
    return

def updateFieldStreetInlets(shapefile_path):
    """
    This function creates fields in a shapefile.
    
    Input Arguments: 
    shapefile_path             --   Path to point shapefile.
    """
    from arcpy import env
    env.workspace = shapefile_path
      
    arcpy.AddField_management(shapefile_path, "ID_scratch", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "Id")  # Delete automatically created field
    arcpy.AddField_management(shapefile_path, "SourceFlow", "FLOAT")
    arcpy.AddField_management(shapefile_path, "ID_AGGREG", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "ID_scratch")
    return

def writefieldsAllNodes(shapefile_path, pointListNear):
    """
    This function creates fields in a shapefile and adds a value from a list
    
    Input Arguments: 
    shapefile_path             --   Path to point shapefile.
    pointListNear    --    list with values
    """
    rows = arcpy.UpdateCursor(shapefile_path)
    counter = 0
    for row in rows:
        row.setValue("ID", pointListNear[counter][0])
        rows.updateRow(row)
        counter += 1
    return

def writefieldsStreetInlets(shapefile_path, pointListNear):
    """
    This function creates fields in a shapefile and adds a field value from a list.
    
    Input Arguments: 
    shapefile_path   --   Path to point shapefile.
    pointListNear    --    list with values
    """
    cnt, rows = 0, arcpy.UpdateCursor(shapefile_path)

    for row in rows:       
        row.setValue("SourceFlow", pointListNear[cnt][8])
        row.setValue("ID_AGGREG", pointListNear[cnt][0])
        rows.updateRow(row)
        cnt += 1
    return

def writefieldsNodes(shapefile_path, pointListNear):
    """
    This function creates fields in a shapefile and adds a value from a list
    
    Input Arguments: 
    shapefile_path           --   Path to point shapefile.
    pointListNear            --    list with values
    """
    from arcpy import env
    env.workspace = shapefile_path

    arcpy.AddField_management(shapefile_path, "ID_scratch", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "Id")      
    arcpy.AddField_management(shapefile_path, "ID", "FLOAT")
    arcpy.AddField_management(shapefile_path, "POP" , "FLOAT")
    arcpy.AddField_management(shapefile_path, "FLOW", "FLOAT")
    arcpy.AddField_management(shapefile_path, "TRENCHD", "FLOAT")
    arcpy.AddField_management(shapefile_path, "HEIGHT", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "ID_scratch")  
    
    counter, rows = 0, arcpy.UpdateCursor(shapefile_path)
      
    for row in rows:
        row.setValue("ID", pointListNear[counter][0])  
        row.setValue("HEIGHT", pointListNear[counter][3])
        row.setValue("FLOW", pointListNear[counter][4])
        row.setValue("POP", pointListNear[counter][6])
        row.setValue("TRENCHD", pointListNear[counter][7])
        rows.updateRow(row)
        counter += 1
    return

def writeFieldNodesPUMPS(shapefile_path, pointListNear):
    """
    This function creates fields in a shapefile and adds a value from a list
    
    Input Arguments: 
    shapefile_path       --   Path to point shapefile.
    pointListNear        --    list with values
    """
    from arcpy import env
    env.workspace = shapefile_path
    
    arcpy.AddField_management(shapefile_path, "ID_scratch", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "Id")           
    arcpy.AddField_management(shapefile_path, "ID", "FLOAT")
    arcpy.AddField_management(shapefile_path, "FLOW", "FLOAT")
    arcpy.DeleteField_management(shapefile_path, "ID_scratch")    
    
    cnt, rows = 0, arcpy.UpdateCursor(shapefile_path)
    
    for row in rows:
        row.setValue("ID", pointListNear[cnt][0])
        row.setValue("FLOW", pointListNear[cnt][3])
        rows.updateRow(row)
        cnt += 1
    return

def initialPrimCalc(firstIteration, origin, nodes, PN, initialPN):
    """
    This functions calculateds the initial prim distances.
    
    Input Arguments: 
    firstIteration         -    Criteria wheter first Iteration or not
    origin                 -    Starting node
    nodes                  -    nodes
    PN                     -    PRIM distances
    initialPN              -    initial PRIM distances
    
    Output Arguments:
    PN                     -    Updated PRIM Distances
    """
    
    if firstIteration == 1:
        idp0, p0, _, _, _, _ = getPns(origin, nodes)  
        
        # Calculate distance from origin to all nodes
        for i in nodes: 
            idp1, p1, flowCriteria = i[0], (i[1], i[2], i[3]), i[8]  
            
            # if not point to itself, point not already connected and no ArchPoint 
            if flowCriteria > 0 and idp0 != idp1:  
                distanceinclSlope, _, _ = distanceCalc3d(p1, p0)    
                distFak = 1                                                         
                PN.append([distanceinclSlope, idp0, idp1, i[1], i[2], i[3], distFak])

        initialPN = fastCopy(PN)  # Make copy 
        firstIteration = 0
        return PN, initialPN, firstIteration
    else:
        return PN, initialPN, firstIteration

def primResultGISList(DrawHouses, VerticGraph, pntsFlowMin, streetVertices, buildings, buildPoints, rasterPoints, edgesList):
    """
    This function reads out the coordinates for calculating GIS polylines.
    
    Input Arguments: 
    DrawHouses              --    Criteria if house connections should be drawn
    VerticGraph             --    Graph from prim Algorithm
    pntsFlowMin             --    Points from prim Algorithm
    streetVertices          --    street vertices from original vertices
    buildings               --    nodes of buildlings
    buildPoints             --    list with buildinges
    rasterPoints            --    raster points
    OnlyExecuteMerge        --    Criteria whter one big merge

    Output Arguments:
    listforCalc             --    [([x,y],[x2,y2]), ...]
    """
    listforcalc = []
    
    # Draw House Connections
    if DrawHouses == 0:
        for node in buildings:
            pt_to1_X, pt_to1_Y, gebListe = node[0], node[1], node[2]

            for house in gebListe:     
                for geb in buildPoints:  # Buildling coordinates
                    if geb[0] == house:
                        ID, pt_from1_X, pt_from1_Y, flow = geb[0], geb[1], geb[2], geb[4]
                        break
                        
                z1 = (ID, [pt_from1_X, pt_from1_Y], [pt_to1_X, pt_to1_Y], flow, 0, 0)  # Add house connection
                listforcalc.append(z1)
          
    # Draw the graph in ArcGIS
    for i in VerticGraph:                       # Get flow, pt_from and pt_to and append to z_list
        pt_from, pt_to, foundFROM, foundTO = i, VerticGraph[i][0], 0, 0

        if pt_from != () and pt_to != ():
            for entry in pntsFlowMin:           # Recalculate position in list and get flow
                if entry[0] == pt_from:
                    flow = entry[4] + entry[6]  # Current Flow + flow from node
                    trenchDepthFrom = entry[7]  # Trench
                    break
                
            for entry in pntsFlowMin: 
                if entry[0] == pt_to:
                    trenchDepthTo = entry[7]    # Trench
                    break

            averageTD = (abs(trenchDepthFrom) + abs(trenchDepthTo)) / 2
               
            # Get pipe diameter
            for entry in edgesList: 
                if (entry[0][0] == pt_from and entry[1][0] == pt_to) or (entry[0][0] == pt_to and entry[1][0] == pt_from):
                    pipeDiameter = entry[4]
                    break

            # Recalculate position in original vertices list and get coordinates 
            for entry in pntsFlowMin:
                if entry[0] == pt_from:
                    foundFROM, korX_from_pt, korY_from_pt = 1, entry[1], entry[2]
                    break
                    
            for entry in pntsFlowMin:
                if entry[0] == pt_to:
                    foundTO, korX_to_pt, korY_to_pt = 1, entry[1], entry[2]
                    break
    
            # If the points in VerticGraph are no StreetVertices but DEM-points, search coordinates in DEM-List
            if foundFROM == 0:
                for DEMPunkt in rasterPoints:
                    if DEMPunkt[0] == pt_from:
                        korX_from_pt, korY_from_pt = DEMPunkt[1], DEMPunkt[2]
                        break
                        
            if foundTO == 0:
                for DEMPunkt in rasterPoints:
                    if DEMPunkt[0] == pt_to:
                        korX_to_pt, korY_to_pt = DEMPunkt[1], DEMPunkt[2]
                        break
                        
            z = (0, [korX_from_pt, korY_from_pt], [korX_to_pt, korY_to_pt], flow, pipeDiameter, averageTD)
            listforcalc.append(z) 
    return listforcalc

def createVirtualDEM(rasterSize, border, buildPoints):
    """
    This function lays a virtual flat DEM over the boudning box of the points plus a border distance.
    
    Input Arguments:
    rasterSize       --    List with aggregated points
    border           --    Size of a tile [m]
    buildPoints      --    List with buildlings
    
    Output Arguments:
    startnode        --    Return top left ID coordinate of tile with highest density
    """
    tupleTopLeft, tupleBottomRight = getTopFleftBottomRightTuple(buildPoints)  # Get top right coordinate
    
    # Add Border
    x_tupleTopLeft, y_tupleTopLeft = (tupleTopLeft[0] - border), (tupleTopLeft[1] + border)
    x_tupleBottomRight, y_tupleBottomRight = (tupleBottomRight[0] + border), (tupleBottomRight[1] - border)
    
    tupleTopLeft = [x_tupleTopLeft, y_tupleTopLeft]
    x_tupleBottomRight = [x_tupleBottomRight, y_tupleBottomRight]
    
    rasterPoints = []

    # add raster points
    width = tupleBottomRight[0] - tupleTopLeft[0]
    height = tupleTopLeft[1] - tupleBottomRight[1]
    nrOfRasterCellsWidht = int(round((width / rasterSize + 1), 0))
    nrOfRasterCellsHeight = int(round((height / rasterSize + 1), 0))
    lineCnt, ID = 0, 1000000  
    
    # Add Raster points until all lines are rastered
    while lineCnt <= nrOfRasterCellsHeight:
        for i in range(nrOfRasterCellsWidht):
            x_ = i * rasterSize + tupleTopLeft[0]
            y_ = tupleTopLeft[1] - lineCnt * rasterSize
            z = 0
            rasterPoints.append((ID, x_, y_, z))  
            ID += 1
        lineCnt += 1
    return rasterPoints

def createPolyLine(polyline_calc, outListStep):
    """
    This function creates polylines for ArcGIS.
    
    Input Arguments: 
    polyline_calc             --   Path to point shapefile.
    outListStep    --    list with values
    """
    from arcpy import env
    env.workspace = outListStep

    point, array, featureList = arcpy.Point(), arcpy.Array(), []  # A list that will hold each of the Polyline objects   

    for feature in polyline_calc:
        # For each point, set the x,y properties and add to the array object.
        coordPair1 = feature[1]
        point.X, point.Y = coordPair1[0], coordPair1[1]
        array.add(point)

        coordPair2 = feature[2]
        point.X, point.Y = coordPair2[0], coordPair2[1]
        array.add(point)   

        polyline = arcpy.Polyline(array)  # Create a Polyline object based on the array of points
        array.removeAll()  # Clear the array for future use
        featureList.append(polyline)  # Append to the list of Polyline objects

    arcpy.CopyFeatures_management(featureList, outListStep)  # Create a copy of the Polyline objects.
    
    arcpy.AddField_management(outListStep, "FLOW", "FLOAT")  # Add Field
    rows = arcpy.UpdateCursor(outListStep)
    arcpy.AddField_management(outListStep, "DIAMETER", "FLOAT")  # Add Field
    rows = arcpy.UpdateCursor(outListStep)
    arcpy.AddField_management(outListStep, "TRENCH", "FLOAT")  # Add Field
    rows = arcpy.UpdateCursor(outListStep)

    cnt = 0
    for row in rows:
        row.setValue("FLOW", polyline_calc[cnt][3])
        row.setValue("DIAMETER", polyline_calc[cnt][4])
        row.setValue("TRENCH", polyline_calc[cnt][5])
        rows.updateRow(row)
        cnt += 1
    return

def assignStreetVertAggregationMode(aggregatetPoints, streetVertices, minTD):
    """
    This function assigns the ID of the street network to the aggregated nodes.
    
    Input Arguments: 
    aggregatetPoints               --    List with coordinates
    streetVertices                 --    Defined interval...
    minTD                          --    Defined classes...

    Output Arguments:
    forSNIP                        --    Dictionary of accumulated Streetvertices with unique ID
    aggregatetPoints               --    aggregatetPoints with corrected ID
    """
    vertexDict = {}  
    
    # Add all streetVertices into vertexDict
    for vertex in streetVertices:
        xVertex, yVertex, zVertex = vertex[1], vertex[2], vertex[3]
        trenchDepth = zVertex - minTD
        vertexDict[vertex[0]] = [vertex[0], xVertex, yVertex, zVertex, 0, 0, 0, 0, 0, [], trenchDepth]

    # Add all agreggatet points into vertexDict with scrapID and replace if point already existing
    for point in aggregatetPoints:
        x1, y1 = point[1], point[2]
        # Check wheter the added point was in streetVertices
        for i in vertexDict:
            streetID, xNear, yNear, zNear = vertexDict[i][0], vertexDict[i][1], vertexDict[i][2], vertexDict[i][3]
            
            # Assign correct ID and FLOW
            if xNear == x1 and yNear == y1:
                point[0] = streetID
                vertexDict[streetID] = point
                vertexDict[streetID][3] = zNear  # Assign height
                vertexDict[streetID][10] = zNear - minTD
    forSNIP = dictToList(vertexDict)  # Dictionary to list
    return forSNIP, aggregatetPoints
                
def dictToList(dictionary):
    """
    This function converts a ditionary to a list
    
    Input Arguments: 
    dictionary    --    Dictionary

    Output Arguments:
    newList       --    List
    """
    newList = []
    for i in dictionary:
        newList.append(dictionary[i])
    return newList 
       
def archPath(dijkstra, distances):
    """
    This functions writes out the edges from the djikstra algorithm of the new path to be added (archPath).
    
    Input Arguments: 
    dijkstra          --    Path of djikstra
    distances         --    Dijkstra-Distances
    
    Output Arguments:
    archPoints        --    List with path
    """
    archPoints, cnter = [], -1
    dijkstra = dijkstra[::-1]  # Invert pathOrigin to get correct flow to origin
    
    # ArchPath: Calculate Distances on the path and the whole path to archPoints
    for i in dijkstra:
        cnter += 1 
        if cnter > 0:
            nextElement = i
            interDist = distances[nextElement] - distances[oldElement]      
            archPoints.append((oldElement, (nextElement, interDist)))  # Form: [(from, (to, dist)), (....)]
        oldElement = i
    return archPoints

def writeToDoc(outListStep_point, name, inList):
    """
    This functions writes to a .txt file.
    
    Input Arguments: 
    outListStep_point    --    Path to folder where .txt file is stored
    name                 --    Name of .txt file
    inList               --    List with entries to write to .txt file    
    """
    outfile = outListStep_point + name + ".txt"
    zList = []
    myDocument = open(outfile, 'w')
    for i in inList:
        z = (str(i) + "\t" + str(inList[i]))
        zList.append(z)

    for item in zList:
        myDocument.write(item)
        myDocument.write("\n")
    myDocument.close()
    return

def writeToDocLine(outListStep_point, name, inList):
    """
    This functions writes to a .txt file.
    
    Input Arguments: 
    outListStep_point    --    Path to folder where .txt file is stored
    name                 --    Name of .txt file
    inList               --    List with entries to write to .txt file    
    """
    outfile = outListStep_point + name + ".txt"
    zList = []
    myDocument = open(outfile, 'w')
    for i in inList:
        z = ({i: inList[i]})
        zList.append(z)
    for item in zList:
        myDocument.write(str(item))
        myDocument.write("\n")
    myDocument.close()
    return

def writeTotxt(outListStep_point, name, inList):
    """
    This functions writes to a .txt file.
    
    Input Arguments: 
    outListStep_point    --    Path to folder where .txt file is stored
    name                 --    Name of .txt file
    inList               --    List with entries to write to .txt file    
    """
    outfile = outListStep_point + name + ".txt"
    myDocument = open(outfile, 'w')
    for item in inList:
        myDocument.write(str(item))
        myDocument.write("\n")
    myDocument.close()
    return

def writeOutPipes(outListStep_point, name, inList):
    """
    This functions writes to a .txt file.
    
    Input Arguments: 
    outListStep_point    --    Path to folder where .txt file is stored
    name                 --    Name of .txt file
    inList               --    List with entries to write to .txt file    
    """    
    outfile = outListStep_point + name + ".txt"
    zList = []
    myDocument = open(outfile, 'w')
    for i in inList:
        zList.append(str(i))
    for item in zList:
        myDocument.write(item)
        myDocument.write("\n")
    myDocument.close()
    return
    
def writeLogFile(outListStep_point, name, inList):
    """
    This function writes out list to doc file in order to import it easy into excel
    
    Input Arguments: 
    outListStep_point    --    Path to folder where .txt file is stored
    name                 --    Name of .txt file
    inList               --    List with entries to write to .txt file    
    """
    outfile = outListStep_point + name + ".txt"
    myDocument = open(outfile, 'w')
    for item in inList:
        myDocument.write(item)
        myDocument.write("\n")
    myDocument.close()
    return

def writeToDocWTPs(outListStep_point, name, WWTPs): 
    """
    This function writes out list to doc file in order to import it easy into excel
    
    Input Arguments: 
    outListStep_point    --    Path to folder where .txt file is stored
    name                 --    Name of .txt file
    WWTPs             --    List with entries to write to .txt file    
    """
    outfile = outListStep_point + "/" + name + ".txt"
    zList = []
    myDocument = open(outfile, 'w')
    zList.append(str("ID") + "\t" + str("PEOPLE"))

    for entry in WWTPs:
        z = (str(entry[0]) + "\t" + str(entry[1]))
        zList.append(z)

    for item in zList:
        myDocument.write(item)
        myDocument.write("\n")
    myDocument.close()
    return

def readRasterPoints(in_FC, anzForID):
    """
    This function reads out all raster points and calculates the raster cell size. The DEM must be 
    right-angled for correct reading out the resolution.

    Input Arguments: 
    in_FC              --    Path to raster points.
    anzForID           --    Number of aggregated nodes to consider in SNIP.
 
    Output Arguments:
    pointListNear      --    list containing all raster points with unique ID, r
    rasterCellSize     --    raster Cell Size
    """
    rows = arcpy.SearchCursor(in_FC) 
    fields = arcpy.ListFields(in_FC) 
    pointListNear, ID = [], 1000000
    for row in rows:
        for field in fields:
            if field.name == "POINT_X":
                X = row.getValue(field.name)
            if field.name == "POINT_Y":
                Y = row.getValue(field.name)
            if field.name == "GRID_CODE":  # Z-Value
                Z = row.getValue(field.name)
        
        pointListNear.append((ID, X, Y, Z))
        ID += 1

    # Assign new ID. This is needed in order that DEM-Points don't get confused with conneciton points
    anzForID = anzForID + 1000000
    pointListNearNewID = []
    
    for i in pointListNear:
        pointListNearNewID.append((anzForID, i[1], i[2], i[3]))
        anzForID = anzForID + 1
    pointListNear = pointListNearNewID
    
    # Calculate rasterCellSize by selecting the two first cell and calculate x-difference
    count = 0
    for cellPoint in pointListNear[:2]:  
        secondX = cellPoint[1]
        secondY = cellPoint[2]
        if count == 1:
            rasterCellSize = math.fabs(secondX - firstX)
            
            if rasterCellSize == 0: # Use X for calculation
                rasterCellSize = math.fabs(secondY - firstY)
            
        firstX = cellPoint[1] 
        firstY = cellPoint[2] 
        count = 1
        
    if rasterCellSize == 0:
        raise Exception("Error: Raster dimension is zero.")
    return pointListNear, rasterCellSize

def readLines(pathInFile):
    """
    This functions reads out lines of a .txt file
    
    Input Arguments: 
    pathInFile       --    Path to .txt file
    
    Output Arguments:
    readLines        --      Statistics
    """
    inputfile = open(pathInFile, 'r')   # Set Path to existing .txt-File with R results
    lineArray = inputfile.readlines()   # Read in single result lines from txtfile into List
    readLines = []                      # Create empty list to fill in whole txt File with results
    position = 0                        
    while position < len(lineArray):    # Iterate through list with line-results by position
        entry = lineArray[position]     # Get line at position
        readLines.append(entry)         # Append line at position to empty List
        position += 1                   # For the loop
    inputfile.close()                   # Close result .txt file 
    return readLines

def readInFastValues(pathInFile):
    """
    This functions reads out lines of a .txt file and creates a dictionary with unique ID
    Input Arguments:
    pathInFile    --    Path to .txt file

    Output Arguments:
    readLines        --      Statistics
    """
    inputfile = open(pathInFile, 'r')   # Set Path to existing .txt-File with R results
    lineArray = inputfile.readlines()   # Read in single result lines from txtfile into List
    readLines = {}                      # Create empty list to fill in whole txt File with results
    position = 0
    while position < len(lineArray):    # Iterate through list with line-results by position
        entry = lineArray[position]     # Get line at position
        entry = entry.replace("\n", "")
        if len(entry) > 0:
            single = entry.split(",")
            anzEntries = len(single)
            newLine = []
            for i in range(anzEntries):
                newLine.append(float(single[i]))
            readLines[position] = newLine  
            position += 1               
        else:
            break
    inputfile.close()                   # Close result .txt file
    return readLines

def readInRasterSize(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    rasterSize        --      Raster Size
    """
    inputfile = open(pathInFile, 'r')       # Set Path to existing .txt-File with R results
    lineArray = inputfile.readlines()       # Read in single result lines from txtfile into List
    txtFileList, position = [], 0           # Create empty list to fill in whole txt File with results
    while position < len(lineArray):        # Iterate through list with line-results by position
        entry = lineArray[position]         # Get line at position
        txtFileList.append(entry)           # Append line at position to empty List
        position += 1                   
    inputfile.close()                       # Close result .txt file 
    rasterSize = float(txtFileList[0][:-1]) # Read out only string
    return rasterSize

def readInStartnode(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile            --    Path to .txt file
    
    Output Arguments:
    rasterSize            --      Raster Size
    """
    inputfile = open(pathInFile, 'r')   # Set Path to existing .txt-File with R results
    lineArray = inputfile.readlines()   # Read in single result lines from txtfile into List
    txtFileList, position = [], 0       # Create empty list to fill in whole txt File with results
    while position < len(lineArray):    # Iterate through list with line-results by position
        entry = lineArray[position]     # Get line at position
        txtFileList.append(entry)       # Append line at position to empty List
        position += 1                   # For the loop
    for i in txtFileList:
        s = i
    f = s.split(",",)
    idStartnode = int(f[0][2:-1])
    x, y = float(f[1][:-1]), float(f[2][:-1])
    startnode = [[idStartnode, x, y]]
    inputfile.close()                   # Close result .txt file 
    return startnode

def readInedgesID(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    liste                --      List with edges
    """
    txtFileList = readLines(pathInFile)
    liste = []
        
    # Read out only string
    for i in txtFileList:
        splittedLine = i.split("],",)

        # add first bracket
        z = []
        firstBracket = splittedLine[0].split()
        
        # Add first element of first bracket
        fstEle = firstBracket[0][2:-1] 
        z.append(int(fstEle))

        for entry in firstBracket[1:-1]:
            og = float(entry[:-1])
            z.append(og) 
            
        lastEle = firstBracket[-1]  # Add last element of first bracket
        z.append(float(lastEle))
        
        # add second bracket
        z1 = []
        firstBracket = splittedLine[1].split()
        fstEle = firstBracket[0][1:-1]      # Add first element of first bracket
        z1.append(int(fstEle))

        for entry in firstBracket[1:-1]:
            og = float(entry[:-1])
            z1.append(og) 
            
        lastEle = firstBracket[-1]          # Add last element of first bracket
        z1.append(float(lastEle))
        
        # add all in between
        for entry in splittedLine[2:-1]:
            z1.append(float(entry[:-1]))

        # add last two entries
        lastEntries = splittedLine[2].split()
        forList = [z, z1, float(lastEntries[0][:-1]), float(lastEntries[1][:-1]), float(lastEntries[2][:-1]), float(lastEntries[3][:-1])]        
        liste.append(forList)               # append who line
    return liste

def readInbuildings(openPath, cluster):
    """
    This functions reads in a .txt file
    
    Input Arguments: 
    openPath         --    Path to .txt file
    
    Output Arguments:
    buildings        --    List with read in lines
    """
    txtFileList = readLines(openPath)
    buildings = []
    
    for i in txtFileList:
        liste = []
        _ = i.split(",", 2)
        if cluster == 1:
            lastElemente = _[2][2:-4].split()
        else:
            lastElemente = _[2][2:-3].split()
        if len(lastElemente) == 1:
            liste.append(float(lastElemente[-1]))
        else:
            for b in lastElemente[:-1]:
                liste.append(float(b[:-1]))
            liste.append(float(lastElemente[-1]))
        buildings.append([float(_[0][1:]), float(_[1]), liste])
    return buildings

def readInWWTPs(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    liste                --   List with edges
    """
    txtFileList = readLines(pathInFile)
    liste = []
    # Read out only string
    for i in txtFileList:
        _ = i.split(",", 2)
        ID = int(_[0][1:])
        flow = float(_[1][1:-2])
        liste.append([ID, flow])               # append who line
    return liste

def readInSewersCurrent(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    liste                --    List with edges
    """
    txtFileList = readLines(pathInFile)
    liste = []
        
    # Read out only string
    for i in txtFileList:
        liste.append(int(i))               # append who line
    return liste

def readInPumps(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    liste                --    List with edges
    """
    txtFileList = readLines(pathInFile)
    liste = []
    
    # Read out only string
    for i in txtFileList:       
        _ = i.split(",", 4)        
        ID = int(_[0][1:])
        flow = float(_[1][1:])
        hDiff = float(_[2][1:])
        lastElement = float(_[3][1:-2])
        liste.append([ID, flow, hDiff, lastElement])               # append who line
    return liste

def readInDictionary(pathInFile):
    """
    This functions reads a .txt file into a dictionary.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    outDictionary        --    Dictionary
    """
    txtFileList = readLines(pathInFile)
    outDictionary = {}  # empty Dictionary
    
    for i in txtFileList:
        spl = i.split(None, 1)  # get index
        index = int(spl[0])
        entries = spl[1].split(",",)
        subDict = {}

        # First Entry
        firstEntry = entries[0].split()
        subDict[int(firstEntry[0][1:-1])] = float(firstEntry[1][:-1])

        # entries in between
        if len(entries) >= 2:
            for entry in entries[1:-1]:
                splitEntry = entry.split(None,)
                subDict[int(splitEntry[0][:-1])] = float(splitEntry[1])
            # Last Entry
            lastEntry = entries[-1].split()
            subDict[int(lastEntry[0][:-1])] = float(lastEntry[1][:-1])
        else:  
            lastEntry = entries[-1].split()
            subDict[int(lastEntry[0][1:-1])] = float(lastEntry[1][:-1])
        outDictionary[index] = subDict
    return outDictionary

def readInSewers(pathInFile):
    """
    This functions reads a .txt file into a dictionary.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    outDictionary        --   Dictionary
    """
    txtFileList = readLines(pathInFile)
    outDictionary = {}  # empty Dictionary
    
    for i in txtFileList:
        spl = i.split(None, 1)  # get index
        index = int(spl[0])
        entries = spl[1].split(",",)
    
        # First Entry
        firstEntry = entries[0].split()
        secondEntry = entries[1].split()

        if firstEntry[0][1:] == '()':
            outDictionary[index] = ((), float(secondEntry[0][:-1]))
        else:
            outDictionary[index] = (int(firstEntry[0][1:]), float(secondEntry[0][:-1]))
    return outDictionary

def readInAggreatedPoints(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    liste                --    List witz aggregated Points
    """
    txtFileList = readLines(pathInFile)
    liste = []  # Read out only string
    
    for i in txtFileList:
        z = []
        splittedLine = i.split(" ", 10)
        z.append(float(splittedLine[0][1:-1]))  # add first element
        
        # add all in between
        for entry in splittedLine[1:-1]:
            z.append(float(entry[:-1]))

        # last element
        lstEle = splittedLine[-1]
        lstEle = lstEle[1:-3]
        lstEle = lstEle.replace(",", "")
        lstEle = lstEle.split()
        
        subList = []
        for i in lstEle:
            subList.append(float(i))
        
        z.append(subList)  # add last elemen
        liste.append(z)
    return liste

def getClosest(PN):
    """
    This functions gets the closest node and deletes the according PRIM-Distance in PN.

    Input Arguments: 
    p0, p1    --    Two points
    
    Output Arguments:
    distance              --    Distance between the points
    slope                 --    Slope between the two points
    heightDiff            --    Height difference betwen the two points
    """
    zahler, minDit = 0, 999999999999 
    for i in PN:
        if i[0] < minDit:
            minDit = i[0]  # shortest distance
            deletPosition = zahler
        zahler += 1
    del PN[deletPosition]  
    return PN, deletPosition

def readInforSNIP(pathRasterPoints):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathRasterPoints           --    Path to .txt file
    
    Output Arguments:
    lines                      --    List with aggregated Points
    """
    txtFileList = readLines(pathRasterPoints)
    lines = []
    for i in txtFileList:
        lineElements = i.split(None, 10)
        allElements = i.split("],",)
        
        # First eight elements
        z = [int(lineElements[0][1:-1]), float(lineElements[1][:-1]), float(lineElements[2][:-1]), float(lineElements[3][:-1]), float(lineElements[4][:-1]), float(lineElements[5][:-1]), float(lineElements[6][:-1]), float(lineElements[7][:-1]), float(lineElements[8][:-1])]
 
        # Copy buildling list
        houses = lineElements[9].split("]",)
        lastList = str(houses[0][1:])

        if len(lastList) > 2:
            lastList = lastList.replace(",", "")
            lastList = lastList.split()
            ele = []
            for i in lastList:
                ele.append(float(i))
        else:
            ele = []
        z.append(ele)
        
        #Add last Elements
        nachList = allElements[1].split()  # Copy entries after list
        z.append(float(nachList[0][:-1]))
        lines.append(z)
    return lines

def readInRasterPoints(pathRasterPoints):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathRasterPoints       --    Path to .txt file
    
    Output Arguments:
    demPnts                --    DEM Points
    """
    txtFileList = readLines(pathRasterPoints)  
    demPnts = []
    for i in txtFileList:
        lineElements = i.split()
        demPnts.append([int(lineElements[0][1:-1]), float(lineElements[1][:-1]), float(lineElements[2][:-1]), float(lineElements[3][0:-1])])
        
        
        # Calculate rasterCellSize by selecting the two first cell and calculate x-difference
    count = 0
    for cellPoint in demPnts[:2]:  
        secondX = cellPoint[1]
        if count == 1:
            rasterCellSize = math.fabs(secondX - firstX) #/ 2 # NEU 12.03.2015
        else:
            firstX = cellPoint[1]
            count = 1
        
    return demPnts, rasterCellSize

def readInbuildPoints(openPath):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    openPath           --    Path to .txt file
    
    Output Arguments:
    buildPoints        --    Buildling points
    """
    txtFileList = readLines(openPath)
    buildPoints = []
    for i in txtFileList:
        lineElements = i.split()
        buildPoints.append([float(lineElements[0][1:-1]), float(lineElements[1][:-1]), float(lineElements[2][:-1]), float(lineElements[3][:-1]), float(lineElements[4][:-1])]) #, float(lineElements[5][:-1])])
    return buildPoints

def readInParameters(pathInputParameterDesign):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInputParameterDesign    --    Path to .txt file
    
    Output Arguments:
    parameters                  --    Parameters
    """
    txtFileList = readLines(pathInputParameterDesign)
    parameters = []
    for i in txtFileList:
        parameters.append(float(i))
    return parameters

def readResultsresult_VerticGraph(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    outList              --    Results SNIP
    """
    txtFileList = readLines(pathInFile)
    
    # Convert Dictionary
    outDict = {}
    for entry in txtFileList:
        i = entry.split()
        if str(i[1][1:-1]) == str(()):
            outDict[int(i[0])] = ((), float(i[2][:-1]))
        else:
            outDict[int(i[0])] = (int(i[1][1:-1]), float(i[2][:-1])) 
    return outDict

def writeTotxtInterResults(outListStep_point, name, inList):
    """
    This function opens a .txt file and write values in it.
    
    Input Arguments: 
    outListStep_point   --    Path to point shapefile.
    name                --    Name (String Value)
    inList              --    list with values 
    """
    outfile = outListStep_point + name + ".txt"
    myDocument = open(outfile, 'w')
    for item in inList:
        myDocument.write(str(item))
        myDocument.write("\n")
    myDocument.close()
    return

def readInstreetVertices(pathstreetVertices):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathstreetVertices    --    Path to .txt file
    
    Output Arguments:
    StreetVertices        --    Street Vertices
    """
    txtFileList = readLines(pathstreetVertices)
    StreetVertices = []
    for i in txtFileList:
        lineElements = i.split()
        StreetVertices.append([int(lineElements[0][1:-1]), float(lineElements[1][:-1]), float(lineElements[2][:-1]), float(lineElements[3][0:-1])])
    return StreetVertices

def readResultspointsPrim(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    outList              --    Results SNIP
    """
    txtFileList = readLines(pathInFile)
    
    # Convert Dictionary
    outList = []
    for entry in txtFileList:
        i = entry.split()
        z = [int(i[0][1:-1]), float(i[1][:-1]), float(i[2][:-1]), float(i[3][:-1]), float(i[4][:-1]), float(i[5][:-1]), float(i[6][:-1]), float(i[7][:-1])] #, float(i[8][:-1]), float(i[9][:-1]), float(i[10][:-1])]
        outList.append(z)
    return outList

def readResultslistWTPs(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    outList              --    Results SNIP
    """
    txtFileList = readLines(pathInFile)
    
    # Convert Dictionary
    outList = []
    for entry in txtFileList:
        i = entry.split()
        z = [int(i[0][1:-1]), float(i[1][:-1])]
        outList.append(z)  
    return outList

def readResultslistPumps(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    outList              --    Results SNIP
    """
    txtFileList = readLines(pathInFile)
    
    # Convert Dictionary
    outList = []
    for entry in txtFileList:
        i = entry.split()
        z = [int(i[0][1:-1]), float(i[1][:-1]), float(i[2][:-1])]
        outList.append(z)
    return outList

def readResultswtpstodraw(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    outList              --    Results SNIP
    """
    txtFileList = readLines(pathInFile)
    
    # Convert Dictionary
    outList = []
    for entry in txtFileList:
        i = entry.split()
        z = [int(i[0][1:-1]), float(i[1][:-1]), float(i[2][:-1]), float(i[3][:-1])]
        outList.append(z) 
    return outList

def readInLine(pathInFile):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    outDiction           --    Ditionary
    """
    txtFileList = readLines(pathInFile)
    outDiction = {}
    
    # Read out list with dictionary
    for i in txtFileList:  # read out linde
        splitted = i.split(None, 1)             # split line into two elements
        index = int(splitted[0][1:-1])          # get index of main dictionary
        splitTwo = splitted[1].split(",")       # split second entry
        subDict = {}                            # create subDictionary
        for sub in splitTwo:
            subSplit = sub.split()              # split second index
            firstChar = subSplit[0][0]          # get first character
            lastChar = subSplit[-1][-1]         # get last character
            
            if firstChar == "{":            
                subIndex = int(subSplit[0][1:-1])
            else:
                subIndex = int(subSplit[0][:-1])
            
            if lastChar == "}":
                subcontent = float(subSplit[1][1:-2])
            else:
                subcontent = float(subSplit[1][:-1])
            
            subDict[subIndex] = subcontent
        outDiction[index] = subDict  # copy single line into graph
    return outDiction


def getStatistics(startnode, sewers, pointsPrim, aggregatetPoints, WWTPs, edgeList, nrOfNeighboursDensity, EWQuantity, buildings, buildPoints):
    '''
    This function calculates statistical information.
    
    Input Arguments: 
    startnode                      --    Start node
    sewers                         --    Dictionary of SNIP
    aggregatetPoints               --    Aggregated Nodes
    WWTPs                          --    WWTPs
    nrOfNeighboursDensity          --    Density
    EWQuantity                     --    Sewage flow per person
    edgeList                       --    List with edges
    
    Output Arguments:
    listWWTPwithAggregatedNodes    --    List with wwtp where the nr of aggregated nodes is added
    '''
    arcpy.AddMessage(" ")
    arcpy.AddMessage("Statistics")
    arcpy.AddMessage("----------")
    totalPublicPipeLength, totalPrivatePipeLenth = 0.0, 0.0

    for i in sewers:
        if sewers[i] != ((), 0):
            totalPublicPipeLength += sewers[i][1]    # sum public pipe length

    # Draw House Connections
    for node in buildings:
        pt_to1_X, pt_to1_Y, gebListe = node[0], node[1], node[2]

        for house in gebListe:     
            for geb in buildPoints:                         # Buildling coordinates
                if geb[0] == house:
                    distance = distanceCalc2d((geb[1], geb[2]), (pt_to1_X, pt_to1_Y))
                    totalPrivatePipeLenth += distance       # sum private pipe length
                    break
    
    arcpy.AddMessage(" ")
    arcpy.AddMessage("Total public pipe Length:  " + str(totalPublicPipeLength))
    arcpy.AddMessage("Total private pipe Length: " + str(totalPrivatePipeLenth))
    arcpy.AddMessage(" ")
        
    # Calculate degree of centralization
    sources = len(aggregatetPoints)
    sinks = len(WWTPs)
    degCen = round((float(sources) - float(sinks)) / float(sources), 4)
    arcpy.AddMessage("Non-weighted Degree of centrality: " + str(degCen))                       # Classical Definition of Z according to Ambros (2009)

    # Weighted Definition of Z
    listWWTPwithAggregatedNodes = getAggregatedNodesinListWWTP(WWTPs, sewers, aggregatetPoints) # Get aggregated nodes
    sumFlow, weightedTerm = 0, 0
        
    for i in listWWTPwithAggregatedNodes:
        sumFlow += i[1]
        weightedTerm += float((i[1]/i[2]))
              
    degCenWeighted = (sumFlow - weightedTerm) / sumFlow
    arcpy.AddMessage("weighted degree of centrality: " + str(degCenWeighted))
    arcpy.AddMessage(" ")
        
    # Calculate average trench depth
    cnt = 0.0
    avreageTrenchDepth = 0
    for i in pointsPrim:
        if i[7] != 0:               # only public sewer
            avreageTrenchDepth += i[7]
            cnt += 1
        
    avreageTrenchDepth = avreageTrenchDepth / cnt
    arcpy.AddMessage("Average Trench Depth: " + str(avreageTrenchDepth))
        
    # Calculate height of pipe network
    cnt = 0.0
    averageHeight = 0
    for i in pointsPrim:
        if i[7] != 0:               # only public sewer
            averageHeight += i[3]
            cnt += 1
        
    averageHeight = averageHeight / cnt
    arcpy.AddMessage("Average height of sewer network: " + str(averageHeight))
           
    # Average wwtp size
    averageWWTPSize, totFlow = 0, 0
        
    cnt = 0.0
    for i in WWTPs:
        averageWWTPSize += i[1]
        totFlow += i[1]
        cnt += 1
    averageWWTPSize = averageWWTPSize / cnt
        
    # Statistics
    statistics = [startnode, sources, sinks, degCen, degCenWeighted, nrOfNeighboursDensity, totalPublicPipeLength, totalPrivatePipeLenth, avreageTrenchDepth, averageHeight, averageWWTPSize]
    return statistics
