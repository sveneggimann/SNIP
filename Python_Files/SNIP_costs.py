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
# Date:      1.1.2015
# Autor:     Eggimann Sven
# ======================================================================================
import math

def calculatePipeCosts(pipeDiameter, distance, averageTrenchDepth, lifeSewers, interestRate, operationCostsPerYear, fc_SewerCost):
    """
    This functions calculates with help of the pipe diameter, the pipe length, average trench Depth and financial parameters,
    the costs of a pipe.
    
    Input Arguments: 
    pipeDiameter              -    Pipe Diameter [cm]
    distance                  -    Pipe Distance [m]
    averageTrenchDepth        -    Average trench depth [m]
    lifeSewers                -    Life span of sewers
    interestRate              -    Real interest rate
    operationCostsPerYear     -    Operation costs per year per meter of pipe [Currency]
    fc_SewerCost              -    Percentage of cost variation
    
    Output Arguments:
    totannuities              -     Annuities of pipe costs, including maintenance
    """
    r = float(interestRate + 1.0)   # calculate r of annuities formula

    # Construction costs
    if pipeDiameter >= 1.2: # If diameter is larger than 1.2m, select parameters of 1200
        pipeDiameter = 1.2

    # CAPEX
    a = 152.51 * pipeDiameter + 173.08                      # Linearly derived function for a & b
    b = 760.31 * pipeDiameter - 78.208                      # Linearly derived function for a & b
    costFactor = 1 + fc_SewerCost                           # Calculate how costs vary
    costPerMeter = a * averageTrenchDepth + b * costFactor  # Calculate cost per meter pipe
    totCost = float(costPerMeter * distance)                # Total costs of whole pipe length
    
    # OPEX
    averageYearlyOperationCosts = operationCostsPerYear * distance  
    totannuities = ((interestRate * r**lifeSewers) / (r**lifeSewers - 1)) * totCost + averageYearlyOperationCosts       # Calculate annuities and add operation costs
    return totannuities

def getPumpCostsDependingOnFlow(Q, heightDifference, pricekWh, nrOfOperatingYears, interestRate):
    """
    This function calculates the costs of a individual pump
    #Source: The sewage pumping handbook, p. 84 ff
    # the pumping investement costs are divided by nr of running years.
    
    Input Arguments: 
    Q                     -    Flow [l / s]
    heightDifference                     --    Slope
    pricekWh              -    Pipe length
    nrOfOperatingYears    -    Strickler coefficient
    interestRate          -    Real Interest rate 
    
    Output Arguments:
    pipeDiameter          -    Needed pipe diameter
    """
    #gravity = 9.81                      # [m / s^2]
    #runninghoursPerYear = 365*24        # [h/year]
    #efficiency = 0.5                    # efficiency of pump plus motor
    
    # Operation costs
    #motorPowerInput = (gravity * Q * heightDifference )/(efficiency*1000)       # [kW]  slower
    #EnergyUsed = motorPowerInput * runninghoursPerYear                          # [kWh] slower
        
    motorPowerInput = (9.81 * Q * heightDifference )/(500)                       # [kW]  faster
    EnergyUsed = motorPowerInput * 8760                                          # [kWh] faster
    
    operationCostsPerYear = EnergyUsed * pricekWh  
    operationCostsOverWholePeriod = operationCostsPerYear * nrOfOperatingYears   # Costs over whole life span
    
    # Error message
    if heightDifference < 0 or operationCostsPerYear < 0:                      
        raise Exception("ERROR: Pumping costs cannot be calculated correctly. " + str(heightDifference) + "" + str(Q))   # Does not make sense if pumped down  

    return operationCostsPerYear, operationCostsOverWholePeriod    #  [running CHF per year + investment costs]

def getPipeDiameter(Q, slope, stricklerC):
    """
    This function calculates the pipe diameter according to Manning-Strickler.
    
    Input Arguments: 
    Q                     -    Flow in pipe
    slope                 -    Slope
    stricklerC            -    Strickler coefficient

    Output Arguments:
    pipeDiameter          -   Needed pipe diameter
    """
    Q = Q/86400.0                                                                           # Convert the flow [m3/day] to [m3/s], 24.0*60.0*60.0  = 86400.0  
    Qmax = 0.8                                                                              # Maximum filling condition
    normDiameterList = (.25, .3, .4, .5, .6, .7, .8, .9, 1, 1.2, 1.5, 2, 2.5, 3, 4, 6, 8)   # [m] Norm pipe diameters 

    #Iterate list with norm diameters until the calculate flow is bigger
    for i in normDiameterList:

        if slope == 0:  # Error if slope is zero
            # If slope is 0 WWTP is pumped and a Diamter of 0.25 is assumed.
            pipeDiameter = normDiameterList[0]
            return pipeDiameter
            #raise Exception("ERROR: Pipe diamater cannot get calculated because slope is zero.") 
        else:
            # Calculate flow in pipe with diameter i [m3/s]
            #Qfull = stricklerC * (i/4.0)**(2.0/3.0) *math.sqrt(abs(slope)) * (math.pi/4.0) * i**2.0                     # slower 
            Qfull = stricklerC * (i/4.0)**(0.6666666666666666) *math.sqrt(abs(slope)) * (0.7853981633974483) * i**2.0    # faster
            
            # If pipe can bear more than flow as input, select this diameter # 80 % condition
            if Qfull * Qmax >= Q:
                pipeDiameter = i
                return pipeDiameter
        
        if i == 8: # Not big enough norm-pipe diameter existing
            pipeDiameter = 10
            return pipeDiameter
    
def costWWTP(flow, EWQuantity, lifeWwtps, interestRate, fc_wwtpOpex, fc_wwtpCapex):
    """
    This function calculates the costs of a wwtp.
    
    Input Arguments: 
    flow                    -    Amount of waste water to be treated [in m3]
    EWQuantity              -    Factor to calculate population equivalent from amount of waste water.
    lifeWwtps               -    Life span of treatment plant
    interestRate            -    Real interest rate [%
    fc_wwtpOpex             -    WWTP operation cost parameter
    fc_wwtpCapex            -    WWTP replacement cost parameter
        
    Output Arguments:
    totalAnnualCosts        -    Total annuities
    """

    r = float(interestRate + 1.0)                                                                       # r of annuities formula
    EW = float(flow) / float(EWQuantity)                                                                # [PE] Calculate flow in population equivalent (Convert liter in EW)
    sensFactor_Operation = 1 + fc_wwtpOpex
    sensFactor_Replacement = 1 + fc_wwtpCapex
    
    # Capex - Annual Operation costs
    replacementCostsPerEW =  13318 * EW**-0.209 * sensFactor_Replacement                                # Source: VSA
    replacementCosts = replacementCostsPerEW * EW 
    annuitiesReplacementCosts = replacementCosts * ((interestRate * r**lifeWwtps)) / (r**lifeWwtps -1)  # Calculate annuities

    # Opex - Annual Operation Costs
    annaulOperationCostsPerEW = 340.82 * EW**-0.171 * sensFactor_Operation                              # Source VSA
    totannaulOperationCosts = annaulOperationCostsPerEW * EW
    totalAnnualCosts = annuitiesReplacementCosts + totannaulOperationCosts                              # Operation Costs & replacement costs
    return totalAnnualCosts

def calculateConnectionCosts(pipeCostI, totPumpCostI, pipeCostII, totPumpCostII, pipeCostIII, totPumpCostIII, WWTPcostsI, interestRate, lifeWwtps, lifeSewers): 
    """
    This function calculates the total replacement costs of the system in case a node is connected to the existing system. 
    This is needed in order to compare the costs of a central connection with the reasonable costs.
    
    Input Arguments: 
    expOrMerge                                --    0: In Expansion module, 1: In Expansion module
    pipeCostI, pumpCostIWholePeriodI          --    Pipe-, pumping costs of option I
    pipeCostII, pumpCostIWholePeriodII        --    Pipe- & pumping costs of option II
    pipeCostIII, pumpCostIWholePeriodVIII     --    Pipe- & pumping costs of option III
    
    WWTPcostsI                                --    WWTP costs option I
    interestRate                              --    Real interest rate [%]
    lifeWwtps                                 --    Life span of treatment plant
    lifeSewers                                --    Life span of sewers
    c_WWTPbefaoreAdding                       --    Costs of WWTP before added flow
    
    Output Arguments:
    costConnection                            --    Costs of lowest central connection
    """
    r = float(interestRate + 1.0)   # calculate r of annuities formula
    var1 = r**lifeSewers
    var2 = var1 -1
    
    # Convert annuities into replacement costs
    totSewerCostI = pipeCostI * var2 / (interestRate * var1)          
    totSewerCostII = pipeCostII * var2 / (interestRate * var1)        
    totSewerCostIII = pipeCostIII * var2 / (interestRate * var1)

    # Calculate total replacement value only of network (including pumps) of the different options
    centralConnectionI = totSewerCostI + totPumpCostI           # total replacement value of pumps and sewer with central connection
    connectionII = totSewerCostII + totPumpCostII               # total replacement value of pumps and sewer without connection
    centralConnectionIII = totSewerCostIII + totPumpCostIII     # total replacement value of pumps and sewer with central connection
    
    # Calculate connection costs of sewers and pumps not including WWTP. The already existing network needs to be substracted (option II).
    costConnectionI = centralConnectionI - connectionII          # Cost central - cost decentral
    costConnectionIII = centralConnectionIII - connectionII      # Cost central - cost decentral

    # Select lowest connection costs
    if costConnectionI < costConnectionIII:
        costConnection = costConnectionI
    else:
        costConnection = costConnectionIII
    return costConnection   

def calculatetotalAnnuities(listWTPs, EW_Q, lifeWwtps, interestRate, pumps, pumpingYears, pricekWh, sewers, flowPoints, edgeList, nodes, stricklerC, lifeSewers, operationCosts, fc_SewerCost, fc_wwtpOpex, fc_wwtpCapex):
    ''' 
    This function calculates the total system costs of a system
    
    Input:
    listWTPs            -    List with wwtps
    EW_Q                -    Waste water per person
    lifeWwtps           -    Lifepsan of wwtp
    interestRate        -    Interest rate
    pumps               -    List with Pumps
    pumpingYears        -    Pump lifespan
    pricekWh            -    Price per kWh   
    sewers              -    Sewers
    flowPoints          -    Nodes with flow
    edgeList            -    List with edges
    nodes               -    Nodes
    stricklerC          -    Strickler Coefficient
    lifeSewers          -    Lifespan of Sewers
    operationCosts      -    Operation costs
    fc_SewerCost        -    Cost factor sewers
    fc_wwtpOpex         -    Cost factor opex WWTP
    
    Output:
    totSystemCosts      -     Total System Costs
    '''
    # calculate WWTPs costs
    completeWWTPCosts = 0
    for i in listWTPs:
        WWTPcostsA1 = costWWTP(i[1], EW_Q, lifeWwtps, interestRate, fc_wwtpOpex, fc_wwtpCapex)   
        completeWWTPCosts += WWTPcostsA1

    # Calculate pump costs
    completePumpCosts = 0
    for pmp in pumps:
        flow, heightDifference = pmp[1], pmp[2] 
        summingPumpCosts, _ = getPumpCostsDependingOnFlow(flow, heightDifference, pricekWh, pumpingYears, interestRate)  # pump is found on path
        completePumpCosts += summingPumpCosts

    # Calculate sewer costs
    completePublicPipeCosts = 0
    for pipe in sewers:
        if sewers[pipe][0] != ():
            oldNode = pipe
            nextNode = sewers[pipe][0]
            
            # Get flow
            for a in nodes:
                if a[0] == oldNode:
                    Q = a[4] + a[8]
                    break
                
            # Get distance, slope
            for edge in edgeList:
                if edge[0][0] == oldNode and edge[1][0] == nextNode:    # Stored inverse, thus slope needs to get inverted
                    distance, slope = edge[2], edge[3] * -1             # distance, # slope needs to be inverted                 
                    break   
    
                if edge[1][0] == oldNode and edge[0][0] == nextNode:
                    distance, slope  = edge[2], edge[3]                 # distance, # slope stays the same      
                    break
            
            # Get Trench Depth
            for punkt in nodes:
                if punkt[0] == oldNode:
                    trenchDepthFrom = punkt[3] - punkt[10]
                    break
            
            for punkt in nodes:
                if punkt[0] == nextNode:
                    trenchDepthTo = punkt[3] - punkt[10]
                    break
            
            averageTrenchDepth = (abs(trenchDepthFrom) + abs(trenchDepthTo)) / 2
            pipeDiameter = getPipeDiameter(Q, slope, stricklerC)
            costsPerYear = calculatePipeCosts(pipeDiameter, distance, averageTrenchDepth, lifeSewers, interestRate, operationCosts, fc_SewerCost) 
            completePublicPipeCosts += costsPerYear
    return completePumpCosts, completeWWTPCosts, completePublicPipeCosts

def costsPrivateSewers(buildings, buildPoints, pipeDiameterPrivateSewer, averageTrenchDepthPrivateSewer, lifeSewers, interestRate, operationCosts, fc_SewerCost):
    '''
    This function calculates the costs of the private sewers. The private sewers are the closest distance to the street network,
    If the street is too far, the whole distance to the building is used.
    
    Input:
    buildings                          -    Buildings
    buildPoints                        -    Coordinates of Buildings
    pipeDiameterPrivateSewer           -    Pipe Diameter
    averageTrenchDepthPrivateSewer     -    Average Trench depth
    lifeSewers                         -    Lifespan of sewers
    interestRate                       -    Interest rate
    operationCosts                     -    Opex
    fc_SewerCost                       -    cost factor sewers
    
    Output:
    totCostPrivateSewer                -    Total replacement value of private sewers
    '''
    costsP_Sewer = 0
    for node in buildings:
        pt_to1_X, pt_to1_Y, gebListe = node[0], node[1], node[2]
    
        for house in gebListe:     
            for geb in buildPoints:  
                if geb[0] == house:
                    _, pt_from1_X, pt_from1_Y, _ = geb[0], geb[1], geb[2], geb[4]
                    break
            
            p0, p1 = (pt_from1_X, pt_from1_Y), (pt_to1_X, pt_to1_Y)
            distance = math.hypot(p0[0] - p1[0], p0[1] - p1[1])     

            privateSewercostsPerYear = calculatePipeCosts(pipeDiameterPrivateSewer, distance, averageTrenchDepthPrivateSewer, lifeSewers, interestRate, operationCosts, fc_SewerCost)
            costsP_Sewer += privateSewercostsPerYear
        
    totCostPrivateSewer = costsP_Sewer * lifeSewers
    return totCostPrivateSewer

def getCostsOfCrossedWWTPs(allNodesToAddToPN, pathBetweenWWTPs, WWTPS_noCon, sewers_NoCon, nodes_noCon, EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement):
    '''
    This function estimates the costs of all crossed wwtps on the path between two wwtps.
    
    Input:
    allNodesToAddToPN    -    All nodes on the path
    pathBetweenWWTPs     -    Path between WWTP
    WWTPS_noCon          -    WWTP
    sewers_NoCon         -    Sewers
    nodes_noCon          -    Nodes
    EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement    -    Cost relevant parameters
    
    Output:
    sumCostcrossedWWTP   -    Costs
    '''
    # Iterate path and get the sum of all flow which flows to WWTPs in the path
    allWWTPsInPath = []             # List to store all crossed wwtp with the flow [[ID, flow]]     
    sumCostcrossedWWTP = 0          # Total costs
    
    # Get all WWTPs in Path
    for i in allNodesToAddToPN:
        for wwtp in WWTPS_noCon:
            if i[0] == wwtp[0]:
                if i[0] in pathBetweenWWTPs:
                    allWWTPsInPath.append([i[0], 0])
                    break

    # Iterate path
    if len(allWWTPsInPath) > 0:
        for i in pathBetweenWWTPs:
            wwtpfound = 0
            
            # get WWTP to which this node flows
            iterate = i
            try:
                while wwtpfound == 0:
                    nextN = sewers_NoCon[iterate]
                    if nextN[0] == ():
                        wwtpfound = 1
                        toWWTP = iterate
                        break
                    iterate = nextN[0]
            except:
                continue # This node was not in network
    
            for wwtp in allWWTPsInPath:
                if wwtp[0] == toWWTP:
                    for n in nodes_noCon:
                        if n[0] == i:
                            fl = n[8]
                            break
                    wwtp[1] += fl
                    break 
        for i in allWWTPsInPath:  
            flowWWTP = i[1]
            costCrossed = costWWTP(flowWWTP, EW_Q, wwtpLifespan, interestRate, fc_wwtpOperation, fc_wwtpReplacement)      
            sumCostcrossedWWTP += costCrossed    
    return sumCostcrossedWWTP
