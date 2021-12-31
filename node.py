# Node definition
class Node:
    def __init__(self, isZone = False):
        # Some nodes in the network are origins/destinations. In traffic
        # engineering, they are often referred as traffic zones.
        self.isZone = isZone
