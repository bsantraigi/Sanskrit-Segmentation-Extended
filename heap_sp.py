import numpy as np
from word_definite import *
import math

def Parent(i):
    return max(0, math.floor((i - 1)/2))

def Left(i):
    return 2*i + 1

def Right(i):
    return 2*(i + 1)

# nominal node class
class Node:
    def __init__(self, id, dist):
        self.dist = dist
        self.id = id
        self.isConflicted = False
        self.src = -1

class Heap:
    # It's a minHeap
    # Nodes are of type Word_definite
    def __init__(self, nodeList):
        self.nodeList = [n for n in nodeList]
        self.len = len(nodeList)
        self.idLocator = {}
        for i in range(self.len):
            self.idLocator[nodeList[i].id] = i
        self.Build()
            
    def Exchange(self, i, j):
        t = self.nodeList[i]
        self.nodeList[i] = self.nodeList[j]
        self.nodeList[j] = t
        self.idLocator[self.nodeList[i].id] = i
        self.idLocator[self.nodeList[j].id] = j
        
    def Decrease_Key(self, node, newDist, src):
        if node.isConflicted:
            return       
        i = self.idLocator[node.id]
        if newDist > node.dist:
            # relaxation not possible
            return
        else:
            node.dist = newDist
            node.src = src
            parent = Parent(i)
            while ((i > 0) and (self.nodeList[parent].dist > self.nodeList[i].dist)):
                self.Exchange(i, parent)
                i = parent
                parent = Parent(i)
                
    def Pop(self):
        if(self.len == 0):
            return None
        if(self.nodeList[0].isConflicted):
            # print("Pop has seen conflict!!!")
            return None
        
        # Remove the entry from the top of the heap
        nMin = self.nodeList[0]
        self.idLocator[self.nodeList[0].id] = -1
        
        # Put the last node on top of heap and heapify
        self.nodeList[0] = self.nodeList[self.len - 1]
        self.idLocator[self.nodeList[0].id] = 0
        self.len -= 1
        self.Min_Heapify(0)
        return nMin
        
    def Min_Heapify(self, i):
        nMin = self.nodeList[i]
        li = Left(i)
        if(li < self.len):
            if(self.nodeList[li].dist < nMin.dist):
                nMin = self.nodeList[li]
                min_i = li                
        ri = Right(i)
        if(ri < self.len):
            if(self.nodeList[ri].dist < nMin.dist):
                nMin = self.nodeList[ri]
                min_i = ri                
        if(nMin.id != self.nodeList[i].id):
            self.Exchange(i, min_i)
            self.Min_Heapify(min_i)
            
    def Delete(self, node):
        i = self.idLocator[node.id]
        self.nodeList[i].isConflicted = True
        self.nodeList[i].dist = np.inf
        self.Min_Heapify(i)
        
    def Build(self):
        self.len = len(self.nodeList)
        for i in range(int(Parent(self.len - 1)) + 1):
            self.Min_Heapify(i)
    
    def Print(self):
        i = 0
        level = 1
        ilimit = 0
        while(i < self.len):
            print('N(%d, %2.1f)' % (self.nodeList[i].id, self.nodeList[i].dist), end = ' ')
            i += 1
            if(i > ilimit):
                print('\n')
                level *= 2
                ilimit += level

'''
nlist = [
    Node(0, 9),
    Node(1, 13),
    Node(2, 7),
    Node(3, 19),
    Node(4, 23),
    Node(5, 47),
    Node(6, 97),
    Node(7, 37),
    Node(8, np.inf),
    Node(9, np.inf),
    Node(10, np.inf)
]

h = Heap(nlist)
h.Print()
h.idLocator
'''
