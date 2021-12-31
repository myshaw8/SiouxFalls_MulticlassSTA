from link import Link
from node import Node
from path import Path
from od import OD

import sys
import traceback
import utils

FRANK_WOLFE_STEPSIZE_PRECISION = 1e-7

class Network:
    """
    This is the class used for transportation networks.  It uses the following
    dictionaries to store the network; the keys are IDs for the network elements,
    and the values are objects of the relevant type:
        node -- network nodes; see node.py for description of this class
        link -- network links; see link.py for description of this class
        ODpair -- origin-destination pairs; see od.py
        path -- network paths; see path.py.  Paths are NOT automatically generated
                when the network is initialized (you probably wouldn't want this,
                the number of paths is exponential in network size.)
                
        The network topology is expressed both in links (through the tail and head
        nodes) and in nodes (enteringLinks and leavingLinks are Node attributes storing
        the IDs of entering and leaving links in a list).
                
        numNodes, numLinks -- self-explanatory
        numZones -- 
        firstThroughNode -- in the TNTP data format, transiting through nodes with
                            low IDs can be prohibited (typically for centroids; you
                            may not want vehicles to use these as "shortcuts").
                            When implementing shortest path or other routefinding,
                            you should prevent trips from using nodes with lower
                            IDs than firstThroughNode, unless it is the destination.
    """

    def __init__(self, networkFile="", demandFile=""):
        """
        Class initializer; if both a network file and demand file are specified,
        will read these files to fill the network data structure.
        """
        self.numNodes = 0
        self.numLinks = 0
        self.numZones = 0
        self.firstThroughNode = 0

        self.node = dict()
        self.link = dict()
        self.ODpair = dict()
        self.path = dict()

        if len(networkFile) > 0 and len(demandFile) > 0:
            self.readFromFiles(networkFile, demandFile)
            
    def formIncidenceMatrix(self):
        
        self.incidenceMatrix = dict()
        for i in self.node:
            self.incidenceMatrix[i] = dict()
            for ij in self.link:
                self.incidenceMatrix[i][ij] = 0
        
        for ij in self.link:
            self.incidenceMatrix[self.link[ij].tail][ij] = 1
            self.incidenceMatrix[self.link[ij].head][ij] = -1
    
        return self.incidenceMatrix
        
    
    def formAdjacencyMatrix(self):
    
        self.adjacencyMatrix = dict()
        for i in self.node:
            self.adjacencyMatrix[i] = dict()
            for j in self.node:
                self.adjacencyMatrix[i][j] = 0
            
        for ij in self.link:
            self.adjacencyMatrix[self.link[ij].tail][self.link[ij].head] = 1
    
        return self.adjacencyMatrix
    
    
    def calculateDegrees(self):
    
        self.degrees = {"indegree": dict(), "outdegree": dict()}
        for i in self.node:
            self.degrees["indegree"][i] = 0
            self.degrees["outdegree"][i] = 0
        
        for ij in self.link:
            self.degrees["indegree"][self.link[ij].head] += 1
            self.degrees["outdegree"][self.link[ij].tail] += 1
        
        return self.degrees
        
    def beckmannFunction(self):
        """
        This method evaluates the Beckmann function at the current link
        flows.
        """
        beckmann = 0
        for ij in self.link:
            beckmann += self.link[ij].calculateBeckmannComponent() 

        return beckmann
        
    def shiftFlows(self, targetFlows, stepSize):
        for ij in self.link:
            self.link[ij].flow += stepSize * (targetFlows[ij] - self.link[ij].flow)

    def shortestPath(self, origin):
        """
        This method finds the shortest path in a network using the label correcting 
        algorithm.
        
        The implementation discussed in the class uses the immediate predecessor 
        node, p(j), in the shortest path as the one of the two labels of node j. 
        In this assignment, you should use entering link labels instead.  The idea is 
        exactly the same, except you are storing the ID of the last *link* in a 
        shortest path to each node.

        Use the 'cost' attribute of the Links to calculate travel times.  These values
        are given -- do not try to recalculate them based on flows, BPR functions, etc.
                
        The backlink and cost labels are both stored in dict's, whose keys are
        node IDs.

        *** BE SURE YOUR IMPLEMENTATION RESPECTS THE FIRST THROUGH NODE!
        *** Travelers should not be able to use "centroid connectors" as shortcuts
        
        You should use the macro utils.NO_PATH_EXISTS to initialize backlink labels 
        (which corresponds to the \phi we used in class), and utils.INFINITY to 
        initialize cost labels. 
        """
        backlink = dict()
        cost = dict()
        
        for i in self.node:
            backlink[i] = utils.NO_PATH_EXISTS
            cost[i] = utils.INFINITY
        cost[origin] = 0
        
        scanList = [origin]
        
        while len(scanList) > 0:
            i = scanList[0]
            scanList.remove(i)
            for forwardLink in self.node[i].leavingLinks:
                followingNode = self.link[forwardLink].head
                tempCost = cost[i] + self.link[forwardLink].cost
                if tempCost < cost[followingNode]:
                    cost[followingNode] = tempCost
                    backlink[followingNode] = forwardLink
                    if followingNode >= self.firstThroughNode:
                        scanList.append(followingNode)

        return (backlink, cost)

    def allOrNothing(self):
        """
        This method generates an all-or-nothing assignment using the current link
        cost values.  It must do the following:
            1. Find shortest paths from all origins to all destinations
            2. For each OD pairs in the network, load its demand onto the shortest
                path found above. 
        The resulting link flows should be returned in the allOrNothing dict, whose
        keys are the link IDs.

        Be aware that the network files are in the TNTP format, where nodes are numbered
        starting at 1, whereas Python starts numbering at 0.  
        """
        allOrNothing = dict()
        for ij in self.link:
            allOrNothing[ij] = 0
            
        for origin in range(1, self.numZones + 1):
            (backlink, cost) = self.shortestPath(origin)
            for OD in [OD for OD in self.ODpair if self.ODpair[OD].origin == origin]:
                curnode = self.ODpair[OD].destination
                while curnode != self.ODpair[OD].origin:
                    allOrNothing[backlink[curnode]] += self.ODpair[OD].demand
                    curnode = self.link[backlink[curnode]].tail
        
        return allOrNothing

    def FrankWolfeStepSize(self, targetFlows, precision = FRANK_WOLFE_STEPSIZE_PRECISION):
        """
        This method returns the step size alpha used by the Frank-Wolfe algorithm.

        The current link flows are given in the self.link[ij].flow attributes, and the
        target flows are given in the targetFlows dictionary.

        The precision argument dictates how close your method needs to come to finding
        the exact Frank-Wolfe step size: you are fine if the absolute difference
        between the optimal step size, and the value returned by your method, is less than
        precision.
        """
        ### You may need to change network.link[ij].flow and use 
        ### network.link[ij].calculateCost() to calculate the link cost at certain flow 
        ### values. Therefore, here we keep a copy of true values of current flows and 
        ### restore if after the bisection method.
        currentFlows = {ij: self.link[ij].flow for ij in self.link}

        low = 0
        high = 1
        while (high-low) > precision:
            alpha = (low + high)/2
            derivative = 0
            for ij in self.link:
                self.link[ij].flow = alpha*targetFlows[ij] + (1-alpha)*currentFlows[ij]
                derivative += (targetFlows[ij] - currentFlows[ij])*self.link[ij].calculateCost()
                
            if derivative < 0:
                low = alpha
            else:
                high = alpha

        ### restore the current link flows
        for ij in self.link:
            self.link[ij].flow = currentFlows[ij]
                                
        return alpha

    def userEquilibriumFW(self,
                        maxIterations = 100,
                        targetGap = 1e-4):
        """
        This method uses the (link-based) convex combinations algorithm to solve
        for user equilibrium.  Arguments are the following:
            maxIterations -- stop after this many iterations have been performed
            targetGap     -- stop once the gap is below this level
        """
        initialFlows = self.allOrNothing()
        for ij in self.link:
            self.link[ij].flow = initialFlows[ij]
            self.link[ij].updateCost()
            
        iteration = 0
        while iteration < maxIterations:
            iteration += 1

            targetFlows = self.allOrNothing()
            
            stepSize = self.FrankWolfeStepSize(targetFlows)
            self.shiftFlows(targetFlows, stepSize)

            for ij in self.link:
                self.link[ij].updateCost()

            ### convergence criterion 
            SPTT = 0
            TSTT = 0
            for ij in self.link:
                SPTT += self.link[ij].cost * self.allOrNothing()[ij]
                TSTT += self.link[ij].cost * self.link[ij].flow
            currentGap = TSTT/SPTT-1
            print("Iteration %d: gap %f; obj fun %f" % (iteration, currentGap, self.beckmannFunction()))
            if currentGap < targetGap:
                break
        

    def readFromFiles(self, networkFile, demandFile):
        """
        Reads network data from a pair of files (networkFile, containing the topology,
        and demandFile, containing the OD matrix) and do some basic checks on
        the input data (validate)"""
        self.readNetworkFile(networkFile)
        self.readDemandFile(demandFile)
        self.validate()
        self.finalize()


    def readNetworkFile(self, networkFileName):
        """
        Reads network topology data from the TNTP data format.  In keeping with
        this format, the zones/centroids are assumed to have the lowest node
        IDs (1, 2, ..., numZones).
        """
        try:
            with open(networkFileName, "r") as networkFile:
                fileLines = networkFile.read().splitlines()

                # Set default parameters for metadata, then read
                self.numNodes = None
                self.numLinks = None
                self.numZones = None
                self.firstThroughNode = 0
                metadata = utils.readMetadata(fileLines)

                try:
                    self.numNodes = int(metadata['NUMBER OF NODES'])
                    self.numLinks = int(metadata['NUMBER OF LINKS'])
                    if self.numZones != None:
                        if self.numZones != int(metadata['NUMBER OF ZONES']):
                            print("Error: Number of zones does not match in network/demand files.")
                            raise utils.BadFileFormatException
                    else:
                        self.numZones = int(metadata['NUMBER OF ZONES'])
                    self.firstThroughNode = int(metadata['FIRST THRU NODE'])
                except KeyError: # KeyError
                    print("Warning: Not all metadata present, error checking will be limited and code will proceed as though all nodes are through nodes.")
                self.tollFactor = float(metadata.setdefault('TOLL FACTOR', 0))
                self.distanceFactor = float(metadata.setdefault('DISTANCE FACTOR', 0))

                for line in fileLines[metadata['END OF METADATA']:]:
                    # Ignore comments and blank lines
                    line = line.strip()
                    commentPos = line.find("~")
                    if commentPos >= 0: # strip comments
                        line = line[:commentPos]

                    if len(line) == 0:
                        continue

                    data = line.split()
                    if len(data) < 11 or data[10] != ';' :
                        print("Link data line not formatted properly:\n '%s'" % line)
                        raise utils.BadFileFormatException

                    # Create link
                    linkID = '(' + str(data[0]).strip() + "," + str(data[1]).strip() + ')'

                    self.link[linkID] = Link(self,
                            int(data[0]), int(data[1]), # head and tail
                            float(data[2]),   # capacity
                            float(data[3]),   # length
                            float(data[4]),   # free-flow time
                            float(data[5]),   # BPR alpha
                            float(data[6]),   # BPR beta
                            float(data[7]),   # Speed limit
                            float(data[8]),   # Toll
                            data[9])          # Link type

                    # Create nodes if necessary
                    if data[0] not in self.node: # tail
                        self.node[int(data[0])] = Node(True if int(data[0]) <= self.numZones else False)
                    if data[1] not in self.node: # head
                        self.node[int(data[1])] = Node(True if int(data[1]) <= self.numZones else False)

        except IOError:
            print("\nError reading network file %s" % networkFile)
            traceback.print_exc(file=sys.stdout)

    def readDemandFile(self, demandFileName):
        """
        Reads demand (OD matrix) data from a file in the TNTP format.
        """
        try:
            with open(demandFileName, "r") as demandFile:
                fileLines = demandFile.read().splitlines()
                self.totalDemand = 0

                # Set default parameters for metadata, then read
                self.totalDemandCheck = None

                metadata = utils.readMetadata(fileLines)
                try:
                    self.totalDemandCheck = float(metadata['TOTAL OD FLOW'])
                    if self.numZones != None:
                        if self.numZones != int(metadata['NUMBER OF ZONES']):
                            print("Error: Number of zones does not match in network/demand files.")
                            raise utils.BadFileFormatException
                    else:
                        self.numZones = int(metadata['NUMBER OF ZONES'])

                except KeyError: # KeyError
                    print("Warning: Not all metadata present in demand file, error checking will be limited.")

                for line in fileLines[metadata['END OF METADATA']:]:
                    # Ignore comments and blank lines
                    line = line.strip()
                    commentPos = line.find("~")
                    if commentPos >= 0: # strip comments
                        line = line[:commentPos]
                    if len(line) == 0:
                        continue

                    data = line.split()

                    if data[0] == 'Origin':
                        origin = int(data[1])
                        continue

                    # Two possibilities, either semicolons are directly after values or there is an intervening space
                    if len(data) % 3 != 0 and len(data) % 4 != 0:
                        print("Demand data line not formatted properly:\n %s" % line)
                        raise utils.BadFileFormatException

                    for i in range(int(len(data) // 3)):
                        destination = int(data[i * 3])
                        check = data[i * 3 + 1]
                        demand = data[i * 3 + 2]
                        demand = float(demand[:len(demand)-1])
                        if check != ':' :
                            print("Demand data line not formatted properly:\n %s" % line)
                            raise utils.BadFileFormatException
                        ODID = str(origin) + '->' + str(destination)
                        self.ODpair[ODID] = OD(origin, destination, demand)
                        self.totalDemand += demand

        except IOError:
            print("\nError reading network file %s" % networkFile)
            traceback.print_exc(file=sys.stdout)

    def validate(self):
        """
        Perform some basic validation checking of network, link, and node
        data to ensure reasonableness and consistency.
        """
        valid = True

        # Check that link information is valid
        for ij in self.link:
            valid = valid and self.link[ij].head in self.node
            valid = valid and self.link[ij].tail in self.node
            if not valid:
                print("Error: Link tail/head not found: %s %s" % (self.link[ij].tail, self.link[ij].head))
                raise utils.BadFileFormatException
            valid = valid and self.link[ij].capacity >= 0
            valid = valid and self.link[ij].length >= 0
            valid = valid and self.link[ij].freeFlowTime >= 0
            valid = valid and self.link[ij].alpha >= 0
            valid = valid and self.link[ij].beta >= 0
            valid = valid and self.link[ij].speedLimit >= 0
            valid = valid and self.link[ij].toll >= 0
            if not valid:
                print("Link %s has negative parameters." % ij)

        # Then check that all OD pairs are in range
        for ODpair in self.ODpair:
            (origin, destination) = (self.ODpair[ODpair].origin, self.ODpair[ODpair].destination)
            valid = valid and origin in self.node
            valid = valid and destination in self.node
            if not valid:
                print("Error: Origin/destination %s not found" % ODpair)
                raise utils.BadFileFormatException
            valid = valid and self.node[origin].isZone == True
            valid = valid and self.node[destination].isZone == True
            if not valid:
                print("Error: Origin/destination %s does not connect two zones" % str(ODpair))
                raise utils.BadFileFormatException
            valid = valid and self.ODpair[ODpair].demand >= 0
            if not valid:
                print("Error: OD pair %s has negative demand" % ODpair)
                raise utils.BadFileFormatException

        # Now error-check using metadata
        if self.numNodes != None and len(self.node) != self.numNodes:
            print("Warning: Number of nodes implied by network file %d different than metadata value %d" % (len(self.node), self.numNodes))
            self.numNodes = len(self.node)
        if self.numLinks != None and len(self.link) != self.numLinks:
            print("Warning: Number of links given in network file %d different than metadata value %d" % (len(self.link), self.numLinks))
            self.numLinks = len(self.link)
        if self.numZones != None and len([i for i in self.node if self.node[i].isZone == True]) != self.numZones:
            print("Warning: Number of zones given in network file %d different than metadata value %d" % (len([i for i in self.node if self.node[i].isZone == True]), self.numZones))
            self.numLinks = len(self.link)
        if self.totalDemandCheck != None:
            if self.totalDemand != self.totalDemandCheck:
                print("Warning: Total demand is %f compared to metadata value %f" % ( self.totalDemand, self.totalDemandCheck))

    def finalize(self):
        """
        Establish the leavingLinks and enteringLinks lists for nodes, initialize flows and
        costs for links and OD pairs.
        """
        # Establish forward/reverse star lists, set travel times to free-flow
        for i in self.node:
            self.node[i].leavingLinks = list()
            self.node[i].enteringLinks = list()
            
        for ij in self.link:
            self.node[self.link[ij].tail].leavingLinks.append(ij)
            self.node[self.link[ij].head].enteringLinks.append(ij)
            self.link[ij].cost = self.link[ij].freeFlowTime + self.link[ij].length * self.distanceFactor + self.link[ij].toll * self.tollFactor
            self.link[ij].flow = 0
            
        for OD in self.ODpair:
            self.ODpair[OD].leastCost = 0
