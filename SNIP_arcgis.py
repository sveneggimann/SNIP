def get_toolboxes(toolboxes, arcpy):
    """
    This functions loads the arcpy module

    Input Arguments: 
    toolboxes          --    Path with ID and distances
    
    Output Arguments:
    arcpy              --    Path with IDs
    """
    sub_folder = "ArcToolbox/Toolboxes/"
    install_dir = arcpy.GetInstallInfo()['InstallDir'].replace("\\", "/")
    tbx_home = os.path.join(install_dir, sub_folder)
    
    for a_tbx in toolboxes:
        try:
            tbx = tbx_home + a_tbx
            arcpy.AddToolbox(tbx)
        except:
            sys.exit()  
    return arcpy

# Importing functions
def rewritePath(inPath):
    ''' Replace Characters in string path to get correct path (avoid bracket problem'''
    
    replaceList = [["\a","REPLACEA"], ["\t", "REPLACET"], ["\r","REPLACER" ], ["\f", "REPLACEF"], ["\b", "REPLACEB"], ["\n", "REPLACEN"], 
                   ["\\", "//"],
                   ["REPLACEA", "/a"], ["REPLACET","/t" ], ["REPLACER", "/r"], ["REPLACEF", "/f"], ["REPLACEB", "/b"], ["REPLACEN", "/n"]]
    
    for c in replaceList:
        inPath = inPath.replace(c[0], c[1])
    return inPath


#import gc; gc.disable() # Disable garbabe collection
import arcpy, os, sys, math, copy, random
from random import choice
from copy import deepcopy
from heapq import heappop, heappush
from datetime import datetime
    
pythonPath = os.getcwd()
sys.path.append(pythonPath)
print("Local Python Path: " + str(pythonPath))

pythonScriptPath = "Q://Abteilungsprojekte/Eng/SWWData/Eggimann_Sven/09-GIS-Python/0-PythonFiles/" 
sys.path.append(pythonScriptPath)                       # Path where python scripts are stored


## Imports
from SNIP_functions import *
from SNIP_astar import *
from SNIP_costs import *
from SNIP_statistics import *

toolboxes = ["Data Management Tools.tbx"]
#arcpy = get_toolboxes(toolboxes, arcpy)
arcpy.env.overwriteOutput = True
desc = arcpy.Describe

# Parameters
inputViaGUI = 0                 # Input via ArcGIS
onlyTxtReadOut = 0              # Criteria if only .txt are read out

if inputViaGUI == 1:
    arcpy.AddMessage("                         ")
    arcpy.AddMessage("Calibrated for ArcGIS GUI")
    arcpy.AddMessage("                         ")
    # ArcGIS System Arguments
    in_FC = sys.argv[1]                                     # Shapefile [line] containing the streets network
    buildings = sys.argv[3]                                 # Shapefile [point] with buildlings
    inDHM = sys.argv[4]                                     # Digital Elevation points. Attention: ID der Raster-Cellen muss greater sein als Anzahl Geb.
    DEMExists = str(sys.argv[5])                            # "YES" (Default) --> DEM Exists, "NO" --> A virtual flat DEM is created. 
    outListFolder = sys.argv[2].replace("\\","/") + "/"     # output folder
    demAlreadyReadOut = 0                                   # DEM is not already read out
    
else:
    print("NO INPUT VIA ARCGIS GUI for Random Case Studies generation")
    print (sys.argv)
    
    # Python System Arguments
    in_FC = sys.argv[5]                                     # Shapefile [line] containing the streets network
    buildings = sys.argv[4]                                 # Shapefile [point] with buildlings
    inDHM = sys.argv[6]                                     # Digital Elevation points. Attention: ID der Raster-Cellen muss greater sein als Anzahl Geb.
    DEMExists = str("YES")                                  # "YES" (Default) --> DEM Exists, "NO" --> A virtual flat DEM is created. 
    outListFolder = sys.argv[7]                             # Path to store txts
    #outListFolder = sys.argv[4]                             # Path to store txts
  
    pathRasterPoints = sys.argv[8]  

    '''
    # SCRAP TEST
    in_FC  = "Q://Abteilungsprojekte/Eng/SWWData/Eggimann_Sven/07-Fallbeispiele/caseStudies/Schwarzenburg/streetNetwork.shp"                                      # Shapefile [line] containing the streets network
    buildings = "Q://Abteilungsprojekte/Eng/SWWData/Eggimann_Sven/07-Fallbeispiele/caseStudies/Schwarzenburg/buildings.shp"   
    inDHM = "Q://Abteilungsprojekte/Eng/SWWData/Eggimann_Sven/07-Fallbeispiele/caseStudies/Schwarzenburg/dem_points.shp"
    outListFolder = "Q://Abteilungsprojekte/Eng/SWWData/Eggimann_Sven/07-Fallbeispiele/caseStudies/Schwarzenburg/SNIP_TESTTER"
    DEMExists = str("YES")                                  # "YES" (Default) --> DEM Exists, "NO" --> A virtual flat DEM is created. 
    pathRasterPoints = ""
    demAlreadyReadOut = 0                                   # Don't read out DEM an read it in from already existing .txt
    '''
    demAlreadyReadOut = 1                                   # IF NCASE DEM IS ALREADY READ OUT
    rasterPoints = readInRasterPoints(pathRasterPoints)     # Read in rasterPoints
    rasterSize = 50
    rasterSizeList = [rasterSize] 
    
outListStep_point = outListFolder + "aggregated_nodes.shp"
outPathPipes = outListFolder + "pipes_original.shp"
outList_nodes = outListFolder + "nodes_original.shp"
outList_pumpes = outListFolder + "pumpes_original.shp"
outList_WWTPs = outListFolder + "WWTPs_original.shp"
aggregatetStreetFile = outListFolder  + "streetGraph.shp"
allNodesPath = outListFolder + "allNodes.shp"
outPath_StartNode = outListFolder + "startnode.shp"


#txtFileInputPath = "Q://Abteilungsprojekte/Eng/SWWData/Eggimann_Sven/07-Fallbeispiele/caseStudies/utzenstorf/"  # Path to folder with .txt files 

#resultPath = txtFileInputPath + "_results" + str(Input_VARIABLE) + "\\"                            # Path to result folder
#resultPath = txtFileInputPath + "_SNIP_VIRTUAL" + "/"     
    
# ========================================================================
# Model parameters
# ========================================================================


runNr = 1                                   # N???

# Sewer related
minTD = 0.25                                # [m] Min trench depth
maxTD = 8                                   # [m] Maximum trench depth
minSlope = 1                                # [%] Criteria of minimum slope without pumps needed
stricklerC = 85                             # [m3/s]

# Algorithm related
f_street = 4                                # [-] How close the sewer follow the road  #I f_street = 1 --> alyays MST (surface MST)
f_merge = 2#3                               # Factor do determine how the OST are reconnected. 
AggregateKritStreet = 1                     # [m] How long the distances on the roads can maximum be before they get aggregated
neighborhood = 1                            # [m] Defines how large the neighbourhood for the a-Star Algorithm
border = 3000                               # How large the virtual dem borders are around topleft and bottom
f_topo = 1.2                                # [-] Factor weighting the dem graph creation for the a* algorithm

# Cost related
resonableCostsPerEW = 7540                  # KT Zurich ihc 7540 pro EG              # [CHF] !Zumutbare AnschlussKosten! If 0: not considered http://www.awel.zh.ch/internet/baudirektion/awel/de/wasserwirtschaft/formular_merkblatt/_jcr_content/contentPar/form_1/formitems/kein_titel_gesetzt_/download.spooler.download.1363330034968.pdf/Richtlinien_Anschlusspflicht.pdf
                                            # The higher, the more are connected                                
pricekWh = 0.20                             # [CHF / kWh] price per kWh of electricity 
pumpingYears = 30                           # [years needed to pump water to compare with sewers
discountYearsSewers = 80                    # [years] How long the pipes do not need to be replaced. Used for discounting calculations
wwtpLifespan = 33                           # [years] How long the wwtps do not need to be replaced. Used for discounting calculations
interestRate = 2                            # [%] Real interest rate

operationCosts = 5                          # [CHF / meter] operation costs per meter pipe per year
pumpInvestmentCosts = 500                   # [CHF] Fix costs of pumps
tileSize = 50                               # [m] for selection of density based starting node

# NEW
fc_SewerCost = 0                             # hOW MUCH more expensivei n percent [ -20% bis + 20%]
fc_wwtpOpex = 0                        # Shift cost curve of sewer constrution up or down!
fc_wwtpCapex = 0                          # 
    
# Cost Assumptions private Sewers
wayPropertyPrivateSewer = 0                 # FIELD construction
pipeDiameterPrivateSewer = 0.25             # Minmium diameter
averageTrenchDepthPrivateSewer = 0.25       # Not deep

# WWTP related
EW_Q = 0.162                                # [m3 / day] 1 EW is equal to 162 liter. This factor must be the same as for the GIS-Files

# WWTP related
EW_Q_Annina = 0.162                         # Scrap
                
# Representation related
drawHouseConnections = 1        # [1 or zero] 1: House Connections are drawn in ArcGIS, 0: House Connections are not drawn in ArcGIS
#intermediateResultsWriteOut =  [10, 20, 30, 40, 50, 60, 70, 80, 90]            # [length of connected nodes] Note all values where an intermediate values is stored

# Representations in ArcGIS
drawFinal, drawIntermediate = True, True           # Draw the 100 % Connected solution # IF intermediate results want to be written out
interestRate = float(interestRate) / 100    # Reformulate real interest rate

'''# Read in txt-Files with informations
pathStreetGraph = txtFileInputPath + "streetGraph.txt"
pathRasterSize = txtFileInputPath + "rastersize.txt"
pathForPrim = txtFileInputPath + "forSNIP.txt"
pathForaggregatedPoints = txtFileInputPath + "aggregatetPoints.txt"
pathForedgesID = txtFileInputPath + "edgesID.txt"
pathRasterPoints = txtFileInputPath + "rasterPoints.txt"
pathstreetVertices = txtFileInputPath + "streetVertices.txt"
pathbuildPoints = txtFileInputPath + "buildPoints.txt"
pathInputParameterDesign = txtFileInputPath + "inputParameters.txt"
pathbuildings = txtFileInputPath + "buildings.txt"

pathList = [pathStreetGraph, pathRasterSize, pathForPrim, pathForaggregatedPoints, pathForedgesID, pathRasterPoints, pathstreetVertices, pathbuildPoints,
            pathInputParameterDesign, pathbuildings]

# Rewrite path in readable form
for path in pathList:
    path = rewritePath(path)
'''
# add parameters into List for SNIP Algorithm
InputParameter = [minTD, maxTD, minSlope, f_merge, 
                  resonableCostsPerEW, neighborhood, f_street, pricekWh, pumpingYears, 
                  discountYearsSewers, interestRate, stricklerC, EW_Q, 
                  wwtpLifespan, operationCosts, pumpInvestmentCosts, f_topo, fc_SewerCost, fc_wwtpOpex, fc_wwtpCapex
                  ]

'''
# Correct flow of Annina which was entered!
for pnt in forSNIP:
    nrOfPerson = pnt[8]/EW_Q_Annina
    pnt[8] = nrOfPerson * EW_Q
    
# Crect flow of Annina
for pnt in aggregatetPoints:
    nrOfPerson = pnt[8]/EW_Q_Annina
    pnt[8] = nrOfPerson * EW_Q

# Crect flow of Annina
for pnt in buildPoints:
    nrOfPerson = pnt[4]/EW_Q_Annina
    pnt[4] = nrOfPerson * EW_Q

'''

# ==========================
# SNIP - Calculations
# ==========================
buildPoints = readBuildingPoints(buildings)                                     # Read out buildings. read ID from Field
arcpy.AddMessage("Building List was read out: Length: " + str(len(buildPoints)))
print("buildings: " + str(buildings))

# If no DEM is given, create flat virtual DEM
if DEMExists == "No":
    rasterSize = 25                                                             # Raster properties for virtual DEM
    rasterPoints = createVirtualDEM(rasterSize, border, buildPoints)            # Create virtual DEM with no heights
    rasterSizeList = [rasterSize]
else:
    anzNumberOfConnections = len(buildPoints)                                   # used for setting new IDs
    if demAlreadyReadOut == 1:
        _ = 0
    else:
        rasterPoints, rasterSize = readRasterPoints(inDHM, anzNumberOfConnections)  # Read out DEM
        rasterSizeList = [rasterSize]
    
nearPoints = readClosestPointsAggregate(buildings)                                  # Read the near_X, near_Y -points of the buildings into a list 

print("TRAH;:S " + str(outListStep_point))
# Aggregate potential streetInlets and create point file. aggregatetPoint become an own ID.
aggregatetPoints, buildings = aggregate(nearPoints, AggregateKritStreet, outListStep_point, minTD)
arcpy.AddMessage("... nodes are aggregated. Total numbers of nodes to be connected: " + str(len(aggregatetPoints)))
print("... nodes are aggregated. Total numbers of nodes to be connected: " + str(len(aggregatetPoints)))

# Update field for the streetInlets
updateFieldStreetInlets(outListStep_point)
writefieldsStreetInlets(outListStep_point, aggregatetPoints)

# Split street line dataset with the streetInlets
arcpy.AddMessage("... the streetdata are splitted")
print("... the streetdata are splitted")

splitStreetwithInlets(in_FC, outListStep_point, aggregatetStreetFile)   # Split with aggregatetPoints

# Update fields in splitted street and add StreetID, the height to each points is assigned from closest DEM-Point
arcpy.AddMessage("... fields are updated")
print("... fields are updated")
updatefieldsPoints(aggregatetStreetFile)

arcpy.AddMessage("... streetVertices are read in and height assigned from aggregatetPoints")
print("... streetVertices are read in and height assigned from aggregatetPoints")
streetVertices = readOutAllStreetVerticesAfterAggregation(aggregatetStreetFile, rasterPoints, rasterSize)

# Build dictionary with vertexes and save which buildlings are connected to which streetInlet
forSNIP, aggregatetPoints = assignStreetVertAggregationMode(aggregatetPoints, streetVertices, minTD)

summe = 0
for i in aggregatetPoints:
    summe += i[8]
arcpy.AddMessage("Summe aggregatetPoints: " + str(summe))     
arcpy.AddMessage(".............----------------------------")    

controlSum = 0
for i in forSNIP:
    if i[8] > 0:
        controlSum += 1
arcpy.AddMessage("ANZALN PRIM MIT BEWOHNER: " + str(controlSum))   

# Write out all relevant nodes
drawAllNodes(streetVertices, allNodesPath)
updateFieldNode(allNodesPath)
writefieldsAllNodes(allNodesPath, streetVertices)

edges = createStreetGraph(aggregatetStreetFile)
arcpy.AddMessage("... edges were read in...")
print("... edges were read in...")

#---Assign id and distance to edges
arcpy.AddMessage("... assign Ids and distance to edges")
print("... assign Ids and distance to edges")
edgesID = addedgesID(edges, streetVertices)


# Add ID to point file and create streetgraph from street lines

streetGraph = appendStreetIDandCreateGraph(edgesID) #alst edges

arcpy.AddMessage("... streetGraph is created")
arcpy.AddMessage("... aggregatetPoints")

aggregatetPoints = assignHighAggregatedNodes(aggregatetPoints, rasterPoints, rasterSize, minTD)  # Assign High to Aggregated Nodes
forSNIP = addBuildingsFarFromRoadTo(aggregatetPoints, forSNIP, rasterPoints, rasterSize, minTD)  # Add all buildings far from the road network to forSNIP
    
# Select random start node
#====================================================
#Create list to store beginngin and end
finalList = {} # Index --> StartNode, number --> length of sewer

while len(finalList) < 1: #DEFINE HOW MANY CALCULATIONS
        arcpy.AddMessage("LENGTH: AGGREGATED NODE: " + str(len(aggregatetPoints)))
        randEntry = choice(aggregatetPoints)
        startnode, startX, startY = randEntry[0], randEntry[1], randEntry[2]
        tileSize = 50 #[m] for selection of density based startnode
        
        _, startnode, startX, startY = densityBasedSelection(aggregatetPoints, tileSize)
        uncertainty = 1
        arcpy.AddMessage("RANDOM STARTNODE: " + str(startnode))
        
        #Alternative if lowest node is starting node
        #startnode = selectLowestNode(aggregatetPoints)  
        
        arcpy.AddMessage("LOWEST STARTNODE: " + str(startnode))

        # Write to .txt files
        #writeTotxt(outListFolder, "intermediateResultsWriteOut", intermediateResultsWriteOut)
        writeTotxt(outListFolder, "inputParameters", InputParameter)
        writeTotxt(outListFolder, "rasterPoints", rasterPoints)
        writeTotxt(outListFolder, "rastersize", rasterSizeList)
        writeTotxt(outListFolder, "buildPoints", buildPoints)
        writeTotxt(outListFolder, "forSNIP", forSNIP)
        writeTotxt(outListFolder, "aggregatetPoints", aggregatetPoints)
        writeTotxt(outListFolder, "forSNIP", forSNIP)
        writeTotxt(outListFolder, "streetVertices", streetVertices)
        writeTotxt(outListFolder, "edgesID", edgesID)
        writeToDoc(outListFolder, "streetGraph", streetGraph)
        writeTotxt(outListFolder, "buildings", buildings)

        # run in arcGIS or not
        if onlyTxtReadOut == 1:
            print("Only read out only the .txt files and abort script ")
            print(sys.argv)
            break
        else:
    
            if startnode in finalList: continue
    
            arcpy.AddMessage("INFO starting node: " + str(startnode))
            
            #---Run Prim and define origin
            #====================================================
            ExpansionTime, MergeTime, VerticGRAPH, pointsPrim, WWTPs, wtpstodraw, pumpList, intermediateResults, edgesID, completePumpCosts, completeWWTPCosts, completePublicPipeCosts, totalSystemCosts = SNIPAlgorithm(runNr, forSNIP, 1, streetGraph, startnode, edgesID, streetVertices, rasterSize, buildPoints, rasterPoints, InputParameter, aggregatetPoints)
            print("Algorithm is successfully calculated.")
            
            arcpy.AddMessage("len intermediateResults:" + str(len(intermediateResults)))
                    
            # Calculate degree of centralization
            sources = len(aggregatetPoints)
            sinks = len(WWTPs)
            degCen = round((float(sources) - float(sinks)) / float(sources), 4)
            
            arcpy.AddMessage("Degree of Centralization")
            arcpy.AddMessage(degCen)

            #TODO: IMplement as well in PYTHON-Script
            sumAggregatedFlow = 0
            for i in aggregatetPoints:
                sumAggregatedFlow += i[8]
                
            summedFlow = 0
            for i in WWTPs:
                summedFlow += i[1]
            
            arcpy.AddMessage("===============================")
            arcpy.AddMessage("SUMMED Aggregated Flow: " + str(sumAggregatedFlow))
            arcpy.AddMessage("===============================")
            
            arcpy.AddMessage("===============================")
            arcpy.AddMessage("SUMMED FLOW: " + str(summedFlow))
            arcpy.AddMessage("===============================")
       
            #---Calculate Length of sub_graph #Stimmt noch nicht!
            sumLength = 0
            for f in VerticGRAPH:
                key = VerticGRAPH[f]            # Dictionary Key des Resultats
                length = key[1]                 # length eines einzelnen Eintrags
                sumLength = sumLength +length
             
            arcpy.AddMessage("sumLent " + str(sumLength))
            finalList[startnode] = sumLength
            
            # ========================
            # Draw final configuration
            # ========================
            if drawFinal == True:
                
                # Write out startnode
                inputStartNod = [[startnode, startX, startY]]
                startNodeToDraw = createPolyPoint(inputStartNod, outPath_StartNode)  # Draw the wtps in ArcGIS
                
                writeStartnode(outPath_StartNode, startNodeToDraw)
                arcpy.AddMessage("startnode is drawn")
                    
                    
                #Draw the graphs in ArcGIS, DrawHouse Connections: First Parameter == 1
                list_GIS = primResultGISList(drawHouseConnections, VerticGRAPH, pointsPrim, streetVertices, buildings, buildPoints, rasterPoints, edgesID)
                
                # TODO: CREATE FILES WITH NETWORK WHICH CAN BE easily read into a GIS
                writeOutPipes(outListFolder, "info_pipes" + str(len(finalList)), list_GIS)
    
                # Distribution
                #==========================================
                arcpy.AddMessage("Final - PipeSitriubiton")
                distributionList = pipeDistribution2d(list_GIS, 10, 10) #liste, interval, klassen
                arcpy.AddMessage(distributionList)
                writeToDoc(outListFolder, "pipeDistribution" + str(len(finalList)), distributionList)
    
                # Create pipes
                createPolyLine(list_GIS, outPathPipes)       # Create pipes
    
                
                # Create wwtps
                #==========================================
                wwtoArDrawn = createPolyPoint(wtpstodraw, outList_WWTPs)  # Draw the wtps in ArcGIS
    
                if wwtoArDrawn == 1:
                    writeARAs(outList_WWTPs, wtpstodraw)
                    arcpy.AddMessage("Final - The wwtps are drawn in ArcGIS")
                else:
                    print("NO WWTP FIELD TO ADD")
                    arcpy.AddMessage("NO WWTP FIELD TO ADD")
    
                
                #Draw the nodes in ArcGIS
                #==========================================
                createPolyPoint(pointsPrim, outList_nodes)
                writefieldsNodes(outList_nodes, pointsPrim)
                arcpy.AddMessage(" Final -The nodes are drawn in ArcGIS")
        
                #Draw the pumps in ArcGIS
                #==========================================
                draw = createPolyPointPump(pumpList, outList_pumpes)
                if draw == True:
                    writeFieldNodesPUMPS(outList_pumpes, pumpList)
                    arcpy.AddMessage("Final - The pumps are drawn in ArcGIS")
                    arcpy.AddMessage("Anzahl Pumpen: " + str(len(pumpList)))
            
            arcpy.AddMessage("============================")
            arcpy.AddMessage(" Final Configuration was drawn")

#Print out list with total length and startnode             
arcpy.AddMessage("finalList")
arcpy.AddMessage(finalList)

#====================================================
arcpy.AddMessage("===========")
arcpy.AddMessage("End Script")
arcpy.AddMessage("===========")
#====================================================
#arcpy.AddMsessage("TTTRAHS")
