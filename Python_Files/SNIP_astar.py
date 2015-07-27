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

## You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#    
# The Djikstra and a* algorithm are adapted from Hetland (2010).
# The algorithm is developed for Python 2.7 and ArcGIS 10.2
#
# Literature
# ----------
# Eggimann Sven, Truffer Bernhard, Maurer Max (2015): To connect or not to connect? 
# Modelling the optimal degree of centralisation for wastewater infrastructures.   
# Water Research, XY, .....
#
# Hetland M.L. (2010): Python Algorithms. Mastering Basic Algorithms in the Python Language. apress.
#
# Contact:   sven.eggimann@eawag.ch
# Version    1.0
# Date:      1.12.2015
# Autor:     Eggimann Sven
# ======================================================================================
import math
from SNIP_functions import *

def distanceCalc2d(p0, p1):
    """
    This functions calculates the euclidian distance in 2d space.

    Input Arguments: 
    p0, p1                --    Points with coordinates. Form: (ID, X, Y, Z)
    
    Output Arguments:
    distance              --    Distance between the points not taking in account the height difference.
    """
    distanz = math.hypot(p0[0] - p1[0], p0[1] - p1[1])         
    return distanz

def topographicFactor(slope, heightDifferencePRIM, f_topo):
    """
    This functions pushes the node away depending on the topography (Topographic Weighting). 
    
    Input Arguments: 
    slope                 --    Slope
    heightDifferencePRIM  --    Height difference [m]
    f_topo                --    exponent
    
    Output Arguments:
    factor                --    Factor
    """
    if heightDifferencePRIM == 0:
        factor = 1
    else:
        factor = float(abs(heightDifferencePRIM)**f_topo)  # Weight distance
    return factor

def distanceDistDEMPath(p0, p1):
    """
    This functions calculates the euclidian distance in 2d space.

    Input Arguments: 
    p0, p1               --    Points with coordinates. Form: (ID, X, Y, Z)
    
    Output Arguments:
    distance             --    Distance between the points taking in account the height difference.
    slope                --    Slope between the two points.
    """
    distancePlanar = math.hypot(p0[0] - p1[0], p0[1] - p1[1])               # Planar distance          
    distanz = math.sqrt(pow(distancePlanar, 2) + pow((p1[2] - p0[2]), 2))   # Distance
    slope = (float((p1[2] - p0[2]))/float(distancePlanar))*100              # Slope
    return distanz, slope

def aStar(rasterSize, rasterPoints, buildPoints, start, end, idp0, idp1, neighborhood, f_topo):
    """
    This function calculates the a* path based on a Input Raster. In case the a* search takes too long because of the many DEM-points,
    ,the function is aborted.

    Input Arguments: 
    rasterSize              --    Size of the raster cells
    rasterPoints            --    Raster points. Form: [[ID, X, Y],...]
    buildPoints             --    List with all building coordinates.
    start, end              --    Id of start node, Id of end Node
    idp0, idp1              --    Coordinates of start and end node
    neighborhood            --    How large the search window is for the a* algorithm    [m]
    
    Output Arguments:
    aStarPath               --    A* Path
    boundingCandidates      --    All rasterpoints in bounding box of a* path
    """
    maxNumberOfDEMPoints = 2000          # Maximal Number of points the a* Algorithm runs on
    shortListFROM, shortListTO = [], []
    
    # Assign FROMNODE and TONODE the closest DEM-cell Points
    for i in rasterPoints:
        if (i[1] - rasterSize < start[0] and start[0] < (i[1] + rasterSize)) and (i[2] - rasterSize < start[1] and start[1] <(i[2] + rasterSize)):
            shortListFROM.append(i)
        if (i[1] - rasterSize < end[0] and end[0] < (i[1] + rasterSize)) and (i[2] - rasterSize < end[1] and end[1] <(i[2] + rasterSize)): 
            shortListTO.append(i)
       
    # Search closest point in shortList of startNode
    fr, dist = [start[0], start[1]], 9999999999

    for i in shortListFROM:
        to = [i[1], i[2]]
        disttoCellPoint = distanceCalc2d(fr, to)

        if disttoCellPoint < dist:
            startCell, startX, startY = i[0], i[1], i[2]
            dist = disttoCellPoint

    # Search closest point in shortList of pointTO
    fr, dist = [end[0], end[1]], 9999999999
    
    for i in shortListTO:
        to = [i[1], i[2]]
        disttoCellPoint = distanceCalc2d(fr, to)

        if disttoCellPoint < dist:
            goalCell, endX, endY = i[0], i[1], i[2]
            dist = disttoCellPoint

    # Create graph
    candidateListWithCosts, _, boundingCandidates = createDEMGraph(rasterSize, rasterPoints, buildPoints, startX, startY, endX, endY, neighborhood, f_topo)

    # If path searching would take too much time, quit without search
    if len(candidateListWithCosts) > maxNumberOfDEMPoints:
        aStarPath = []
        boundingCandidates = []
        return aStarPath, boundingCandidates

    # Calculate a*-path
    path = aStarAlgorithm(candidateListWithCosts, startCell, goalCell, startX, startY, endX, endY, rasterPoints)
    path = path[::-1]

    #Write path out with coordinates
    if len(path) > 2:                                                #if no path was found, ArchPath list is equal to MST
        pathCoordinates = []
        for i in path:
            for pnt in boundingCandidates:
                if i == pnt[0]:
                    pathCoordinates.append([pnt[0], pnt[1], pnt[2], pnt[3]]) #ID, X, Y, Z
                    break
    
        #Correct first and last point in path as this path is still a cell-coordinate
        pathCoordinates[0] = [idp0, start[0], start[1], start[2]]       # Coordinate of pointFROM
        pathCoordinates[-1] = [idp1, end[0], end[1], end[2]]            # Coordinate of pointFROM

    if len(path) == 2:  # If no path is found, and empty archPath is returned
        return [], [] 
    
    #Convert point archPathList into edges
    aStarPath, iterator = [], 0
    for i in pathCoordinates:
        if iterator == 1:
            pt_aktuell = i
            distan, slope = distanceDistDEMPath((pt_vorher[1],pt_vorher[2], pt_vorher[3]), (pt_aktuell[1], pt_aktuell[2], pt_aktuell[3]))   # Calculate distance
            aStarPath.append([pt_vorher[0], [pt_aktuell[0], distan, slope], pt_vorher[3], pt_aktuell[3]])            
            pt_vorher = i
        else:
            pt_vorher, iterator= i, 1
    return aStarPath, boundingCandidates #Return a DEM-ArchPathList

def readValuesInBoundingBox(rasterSize, ID_Point, x, y, pointList, buildingscells):
    """
    This function collects only the raster points in the bounding box of two points which have no buildling on it.

    Input Arguments: 
    rasterSize         --    Size of raster cell
    ID_Point           --    Point ID
    x,y                --    Coordinates of ID_Point
    pointList          --    List with all points in the bounding box 
    buildlingscells    --    All cells where there is a buildling located
 
    Output Arguments:
    candiateList       --    List with candidates [[id, x,y, z],...]
    """
    candiateList, X_northWest, Y_northWest, X_southEast, Y_southEast = [], x - rasterSize, y + rasterSize, x + rasterSize, y - rasterSize

    for point in pointList:
        idp, x, y = point[0], point[1], point[2]

        if ID_Point != idp:
            if x >= X_northWest and x <= X_southEast and y <= Y_northWest and y >= Y_southEast:
                canCopy = 1
                
                for i in buildingscells:    #Check if the cell to which the connect would take place is a buildlingCell
                    if i == idp:
                        canCopy = 0
                        break
                    
                if canCopy == 1:
                    candiateList.append(point)   
    return candiateList

def checkQuadrant(startX, startY, endX, endY):
    """
    This function only collects the cell points in the bounding box of the two points. 

    Input Arguments: 
    startX, startY        --    X & Y Coordinates of startnode  
    endX, endY            --    X & Y Coordinates of endnode  
    
    Output Arguments:
    candiateList          --    List with candidates [[id, x,y, z],...]
    """
    #Only select those points which are within the bounding box of the start and end
    if startX <= endX and startY <= endY: #Quadrant 1
        quadrantSituation = 1
        return quadrantSituation
                
    if startX <= endX and startY >= endY: #Quadrant 2
        quadrantSituation = 2
        return quadrantSituation
                
    if startX >= endX and startY >= endY: #Quadrant 2
        quadrantSituation = 3
        return quadrantSituation

    if startX >= endX and startY <= endY: #Quadrant 4
        quadrantSituation = 4
        return quadrantSituation

def createDEMGraph(rasterSize, rasterPoints, buildingPoints, startX, startY, endX, endY, neighborhood, f_topo):
    """
    This creates a graph out of DEM points.

    Input Arguments: 
    rasterSize           --    Size of raster cell
    rasterPoints         --    List with raster points
    buildingPoints       --    Coordinates of buildings
    startX, startY       --    Start Coordinate 
    endX, endY           --    End Coordinate
    neighborhood         --    Neighborhood Criteria
 
    Output Arguments:
    pnts                 --    updated nodes
    dictionaryGraph      --    Graph made out of DEM Points
    buildingcells        --    Raster cells where there is a buildling on it
    pointListBoundingBox --    All pnts within a bounding box
    """
    dictionaryGraph, buildingcells, pointListBoundingBox = {}, [], []
    quadrant = checkQuadrant(startX, startY, endX, endY)                # iterate rasterPoints to read out only these points which lie in the bounding box
     
    #Variable to define neighborhood which is looked at as well and thus DEM-Points selected
    if quadrant == 1:
        for i in rasterPoints:
            if i[1] >= startX - neighborhood and i[1] <= endX + neighborhood and i[2] >= startY - neighborhood and i[2] <= endY + neighborhood:
                pointListBoundingBox.append(i)
    elif quadrant == 2:
        for i in rasterPoints:
            if i[1] >= startX - neighborhood and i[1] <= endX + neighborhood and i[2] <= startY + neighborhood and i[2] >= endY - neighborhood:
                pointListBoundingBox.append(i)
    elif quadrant == 3:
        for i in rasterPoints:
            if i[1] <= startX + neighborhood and i[1] >= endX - neighborhood and i[2] <= startY + neighborhood and i[2] >= endY - neighborhood:
                pointListBoundingBox.append(i)
    elif quadrant == 4:
        for i in rasterPoints:
            if i[1] <= startX + neighborhood and i[1] >= endX - neighborhood and i[2] >= startY - neighborhood and i[2] <= endY + neighborhood:
                pointListBoundingBox.append(i)
        
    #Check which cells are house-cells and store in list
    for i in pointListBoundingBox:
        currID, currX, currY, currZ = i[0], i[1], i[2], i[3]

        #Check if building is within a cell
        for geb in buildingPoints:
            xGeb, yGeb = geb[1], geb[2]
            if currX - 0.5*rasterSize <= xGeb <= currX + 0.5*rasterSize and currY - 0.5*rasterSize <= yGeb <= currY + 0.5*rasterSize:   # Check if there is a buildling within raster cell
                if currX == startX and currY == startY or currX == endX and currY == endY:                                              # Prevent that the starting or end node are classified as buildings
                    _ = 0
                else:
                    buildingcells.append(currID)
                    break

    #create graph       
    for i in pointListBoundingBox:
        listWithCosts, isBuildingCell, currID, currX, currY, currZ = [], 0, i[0], i[1], i[2], i[3]

        #Check if cell is a building-cell
        for cellID in buildingcells:
            if cellID == currID:
                isBuildingCell = 1
                break

        if isBuildingCell == 0:
            candidateList = readValuesInBoundingBox(rasterSize, currID, currX, currY, pointListBoundingBox, buildingcells)

            #append distance., slope and costs to neighbours
            for candidate in candidateList:
                toID, toX, toY, toZ = candidate[0], candidate[1], candidate[2], candidate[3]
                
                if currX != toX and currY != toY:# new 7.01.2015
                    heightDiff = toZ - currZ
                    distanz2d = math.hypot(currX - toX, currY - toY)    
                    distanz3d = math.sqrt(pow(distanz2d, 2) + pow(heightDiff, 2))
                    slope = float(heightDiff) / float(distanz3d)
                    factor = topographicFactor(slope, heightDiff, f_topo)   
                    weightedDistance = distanz3d * factor                   
                    listWithCosts.append([toID, toX, toY, toZ, distanz3d, slope, weightedDistance])        

        #create graph
        toDictionary = {}
        for i in listWithCosts:
            toDictionary[i[0]] = i[6]           # store costs in dictionary
        dictionaryGraph[currID] = toDictionary           
    return dictionaryGraph, buildingcells, pointListBoundingBox

def aStarAlgorithm(G, start, goal, startX, startY, endX, endY, listWithCoordinates):
    """
    A-Star Algorithm to find shortest path between two points.

    Input Arguments: 
    G                   --    DEM Graph
    start               --    Starting node
    goal                --    Ending nod
    startX, startY      --    Start Coordinate 
    endX, endY          --    End Coordinate
    listWithCoordinates --    List with nodes
 
    Output Arguments:
    backPath            --    Shortest Path
    """
    backPath = []
    
    # Define aStarHeuristic function. A straight line aStar-Heuristic is choosen.
    def aStarHeuristic(from_node):
        for i in listWithCoordinates: 
            if i[0] == from_node:
                dist2d = math.hypot(i[1] - endX, i[2] - endY)
                return dist2d
    
    P, Q = {}, [(aStarHeuristic(start), None, start)]
    
    while Q:                                                      # Still unprocessed nodes?
        smallesvalue, cnt = 99999999, 0
        for i in Q:
            if i[0] < smallesvalue:
                smallesvalue = i[0]
                d, p, u = i[0], i[1], i[2]
                smallesPosition = cnt
            cnt += 1
        del Q[smallesPosition]
        
        if u in P: continue                                       # Already visited? Skip it
        P[u] = p                                                  # Set path predecessor
        if u == goal:                                             # read path back
            if u == start:                                        # If the start and end point are within the same cell
                break                                             # exit loop
            backPath.append(u)
            entry = None
            while entry != start:
                entry = P[u]
                u = entry
                backPath.append(entry)
            return backPath                                       # Arrived! Path was found by a*-algorithm
        
        for v in G[u]:                                            # Go through all neighbors 
            w = G[u][v] - aStarHeuristic(u) + aStarHeuristic(v)   # Modify weight of aStarHeuristic
            Q.append((d + w, u, v))

    backPath = [start, goal]                                      # If no path is found, create shortest path
    return backPath                                               # No Path was found a*
