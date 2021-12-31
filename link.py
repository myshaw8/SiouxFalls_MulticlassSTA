class Link:
    """
    Class for network links.
    """

    def __init__(self, network, tail, head, capacity = 99999, length = 99999, freeFlowTime = 99999, alpha = 0.15, beta = 4, speedLimit = 99999, toll = 0, linkType = 0):
        """
        Initializer for links; note default values for parameters if not specified.
        """
        # Params you will use in this coding assginment.
        self.network = network
        self.tail = tail
        self.head = head

        # Params you will use in future coding assginments. You don't need to
        # understand their meaning for now.
        self.capacity = capacity
        self.length = length
        self.freeFlowTime = freeFlowTime
        self.alpha = alpha
        self.beta = beta
        self.speedLimit = speedLimit
        self.toll = toll
        self.linkType = linkType
        
    def calculateCost(self):
        """
        Calculates the cost of the link using the BPR relation.
        This cost is returned by the method and NOT stored in the cost attribute.
        """   
        vcRatio = self.flow / self.capacity
        # Protect against negative flows, 0^0 errors.
        if vcRatio <= 0:
            return self.freeFlowTime

        travelTime = self.freeFlowTime * (1 + self.alpha * pow(vcRatio, self.beta))

        return travelTime
        
    def calculateBeckmannComponent(self):
        """
        Calculates the integral of the BPR function for the link, for its
        contribution to the sum in the Beckmann function.
        """
        vcRatio = self.flow / self.capacity
        # Protect against negative flows, 0^0 errors.
        if vcRatio <= 0:
            return 0
        
        beckmannComponent = self.flow * self.freeFlowTime * (1 + self.alpha / (self.beta + 1) * pow(vcRatio, self.beta))
        
        return beckmannComponent

    def updateCost(self):
        self.cost = self.calculateCost()
