# ======================================================================================
# Copyright 2014  Swiss Federal Institute of Aquatic Science and Technology
#
# This file is part of SNIP (Sustainable Network Infrastructure Planning)
# SNIP is used for determining the optimal degree of centralization for waste
# water infrastructures.
# 
# Where find detailed information about SNIP:
#
# Eggimann Sven, Truffer Bernhard, Maurer Max (2015): To connect or not to connect? 
# Modelling the optimal degree of centralisation for wastewater infrastructures.   
# Water Research, XY, .....
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
    
# Contact: sven.eggimann@eawag.ch
# Version 1.0
# Date: 1.12.2015
# 
# ======================================================================================

import sys, os
from SNIP_functions import *


print("Imported SNIP_FUNCTIONsI.py")


'''def getAggregatedNodesinListWWTP(WWTPs, pipeNetwork, aggregatedNodes):
    """
    This function makes a breath search for each wwtp to get nr of aggregated nodes to listWWTPs.
    
    Input Arguments: 
    WWTPs                             --    list with WWTPs
    pipeNetwork                          --    Sewer pipe network
    
    Output Arguments:
    listWWTPwithAggregatedNodes          --    List with wwtp where the nr of aggregated nodes is added
    """
    listWWTPwithAggregatedNodes = []                    # Form: ID, total Flow, total nr of aggregated nodes
    
    testFUNCTION("ACHTUNGGGGGGGGGGGGG")

    for wwtp in WWTPs:
        allNodes = breathSearch(wwtp[0], pipeNetwork)
        onlyInhabitedNodes, nrOfNodes = [], 0           # List to store only inahbited nodes

        # Because of arch points which are uninhabited
        for node in allNodes:
            for i in aggregatedNodes:
                if i[0] == node:
                    nrOfNodes += 1
                    break
        z = [wwtp[0], wwtp[1], nrOfNodes]
        listWWTPwithAggregatedNodes.append(z)

    return listWWTPwithAggregatedNodes
'''
   
def getStatistics(startnode, VerticGRAPH, pointsPrim, aggregatetPoints, WWTPs, edgeList, nrOfNeighboursDensity, EWQuantity, buildings, buildPoints):
    '''
    This function calculates statistical information.
    
    Input Arguments: 
    startnode                      --    Start node
    VerticGRAPH                    --    Dictionary of SNIP
    aggregatetPoints               --    Aggregated Nodes
    WWTPs                          --    WWTPs
    nrOfNeighboursDensity          --    Density
    EWQuantity                     --    Sewage flow per person
    edgeList                       --    List with edges
    
    Output Arguments:
    listWWTPwithAggregatedNodes    --    List with wwtp where the nr of aggregated nodes is added
    '''
     # Calculate total length of pipes
    # -------------------------------
    totalPublicPipeLength, totalPrivatePipeLenth = 0.0, 0.0
    
    # Calculate public Sewer Length
    for i in VerticGRAPH:
        if VerticGRAPH[i] != ((), 0):
            totalPublicPipeLength += VerticGRAPH[i][1]    # sum public pipe length

    # Draw House Connections
    for node in buildings:
        pt_to1_X, pt_to1_Y, gebListe = node[0], node[1], node[2]

        for house in gebListe:     
            for geb in buildPoints:  # Buildling coordinates
                if geb[0] == house:
                    # TwoD Distance as approximation
                    distance = distanceCalc2d((geb[1], geb[2]), (pt_to1_X, pt_to1_Y))
                    totalPrivatePipeLenth += distance    # sum private pipe length
                    break

            
    print("Total public pipe Length:  " + str(totalPublicPipeLength))
    print("Total private pipe Length: " + str(totalPublicPipeLength))
        
    # ----------------------------------
    # Calculate degree of centralization
    # ----------------------------------
    
    # Classical Definition of Z
    sources = len(aggregatetPoints)
    sinks = len(WWTPs)
    degCen = round((float(sources) - float(sinks)) / float(sources), 4)
        
    print("Non-weighted Degree of centrality: " + str(degCen))
    print("Nr of Sources: " + str(sources))
    print("Nr of sinks: " + str(sinks))
        
    # Weighted Definition of Z
    listWWTPwithAggregatedNodes = getAggregatedNodesinListWWTP(WWTPs, VerticGRAPH, aggregatetPoints) # Get aggregated nodes
    sumFlow, weightedTerm = 0, 0
    test = 0
        
    for i in listWWTPwithAggregatedNodes:
        test += i[2]
        sumFlow += i[1]
        weightedTerm += float((i[1]/i[2]))
              
    degCenWeighted = (sumFlow - weightedTerm) / sumFlow
    print("weighted degree of centrality: " + str(degCenWeighted))
        
    # Calculate average trench depth
    cnt = 0.0
    avreageTrenchDepth = 0
    for i in pointsPrim:
        if i[7] != 0:               # only public sewer
            avreageTrenchDepth += i[7]
            cnt += 1
        
    avreageTrenchDepth = avreageTrenchDepth / cnt
    print("Average Trench Depth: " + str(avreageTrenchDepth))
        
    # Calculate height of pipe network
    cnt = 0.0
    averageHeight = 0
    for i in pointsPrim:
        if i[7] != 0:               # only public sewer
            averageHeight += i[3]
            cnt += 1
        
    averageHeight = averageHeight / cnt
    
    print("Average height of sewer network: " + str(averageHeight))
        
        
    # Average wwtp size
    # ----------------------------------
    averageWWTPSize, totFlow = 0, 0
        
    cnt = 0.0
    for i in WWTPs:
        averageWWTPSize += i[1]
        totFlow += i[1]
        cnt += 1
    averageWWTPSize = averageWWTPSize / cnt
        
    # Calculate median wwtpSize
    medianWWTPSize = 0
    
    # Sort list according to size
    sortedWWTPs = sorted(WWTPs, key=lambda WWTPs: WWTPs[1]) #Sort tuple according to size
    anzWWTP = len(sortedWWTPs)
    boolTest = anzWWTP % 2

    # Get median
    '''if boolTest == 0 and anzWWTP > 1: # If odd number
        medianPosition= int(((anzWWTP-1)/2.0) - 0.5)
        medianWWTPSize = sortedWWTPs[medianPosition]
    else:
        medianPosition= int((anzWWTP-1)/2.0)
        medianWWTPSize = (sortedWWTPs[medianPosition][1] + sortedWWTPs[medianPosition + 1][1])/2
    '''
    medianWWTPSize = 0
    # Calculate pipe size distribution
    # ----------------------------------
    p_to10EW, p_10to100EW, p_100to1000EW, p_over1000EW = 0, 0, 0, 0
    to10EW, von10To100EW, von100To1000EW, over1000EW = 0, 0, 0, 0
        
    for i in WWTPs:
        if (i[1]*1000)/EWQuantity <= 10: # Flow in m3 divided by flowPerEW
            to10EW += i[1]
        if (i[1]*1000)/EWQuantity > 10 and (i[1]*1000)/EWQuantity <= 100:
            von10To100EW += i[1]
        if (i[1]*1000)/EWQuantity > 100 and (i[1]*1000)/EWQuantity <= 1000:
            von100To1000EW += i[1]
        if (i[1]*1000)/EWQuantity > 1000:
            over1000EW += i[1]
            
    p_to10EW = (to10EW / totFlow) * 100
    p_10to100EW = (von10To100EW) * 100
    p_100to1000EW = (von100To1000EW / totFlow) * 100
    p_over1000EW = (over1000EW / totFlow) * 100
        
    # Statistics
    statistics = [startnode, sources, sinks, degCen, degCenWeighted, nrOfNeighboursDensity, totalPublicPipeLength, totalPrivatePipeLenth, avreageTrenchDepth, averageHeight, averageWWTPSize, medianWWTPSize, p_to10EW, p_10to100EW, p_100to1000EW, p_over1000EW]
        
    return statistics