import copy
# get_kmers
# Input: A string Genome and an integer k.
# Output: The set of all k-mers in Genome.
from random import randint


def get_kmers(k, Genome):
    kmers = []
    for i in range(len(Genome)-k+1):
        kmers.append(Genome[i:i+k])
    return kmers   

# PathToGenome
# Input: A sequence of k-mers Path = (Pattern1, …, Patternn) 
# such that the last k-1 symbols of Patterni are equal to the 
# first k-1 symbols of Patterni+1 for 1 ≤ i ≤ n-1.  
# Output: The string spelled by Path.
def PathToGenome(Path):
    Genome = Path[0]
    for i in range(1, len(Path)):
        Genome += Path[i][len(Path[i])-1]
    return Genome

# AdjacencyList
# Input: A collection of k-mers Patterns.
# Output: The adjacency list corresponding to Patterns, in the form 
# of a dictionary mapping each k-mer Pattern to a list of k-mers 
# that have an overlap of length k-1 with Pattern.
def AdjacencyList(kmers):
    adjList = {}

    for i in range(len(kmers)):
        adjList[kmers[i]] = []

    for i in range(len(kmers)):
        suffix = kmers[i][1:len(kmers[i])]
        for j in range(len(kmers)):
            prefix = kmers[j][0:len(kmers[j])-1]
            if (prefix == suffix):
                adjList[kmers[i]].append(kmers[j])

    for i in range(len(kmers)):
        if adjList[kmers[i]] == []:
            del adjList[kmers[i]]
    
    order = {'C': 0, 'G': 1, 'A': 2, 'T': 3}
    return dict(sorted(adjList.items(), key=lambda x: [order[c] for c in x[0]]))

# deBruijn
# Input: A collection of k-mers Patterns.
# Output: The de Bruijn graph corresponding to Patterns, in the form of an adjacency list.
def deBruijn(kmers):
    adjList = {}
    
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        
        if prefix not in adjList:
            adjList[prefix] = []
        adjList[prefix].append(suffix)  
    
    return adjList
    #order = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    #return dict(sorted(adjList.items(), key=lambda x: [order[c] for c in x[0]]))

def shift(n, list):
    return list[n:] + list[:n]

# EulerianCycle
# Input: An adjacency list of a directed graph that contains an Eulerian cycle.
# Output: An Eulerian cycle in this graph, in the form of a list of nodes
def eulerianCycle(adjList):
    keys = list(adjList.keys())
    trackList = copy.deepcopy(adjList)
    currentNode = keys[randint(0, len(keys)-1)]
    currentcycle = []
    thiscycle = []

    while any(trackList.values()):
        currentcycle = thiscycle

        while trackList[currentNode] != []:
            if currentcycle == []:
                currentcycle.append(currentNode)
            previousNode = currentNode
            choice = randint(0, len(trackList[currentNode]) - 1)
            currentNode = trackList[currentNode][choice]
            currentcycle.append(currentNode)
            trackList[previousNode].remove(trackList[previousNode][choice])
        
        thiscycle = currentcycle
    
        for key in thiscycle:
            if trackList[key] != []:
                currentNode = key
                n = thiscycle.index(currentNode)
                thiscycle = shift(n, thiscycle[:-1]) + [currentNode]
                break
        
        
        
    
    return thiscycle


# EulerianPath
# Input: An adjacency list of a directed graph that contains an Eulerian path.
# Output: An Eulerian path in this graph, in the form of a list of nodes
def eulerianPath(adjList):
    keys = list(adjList.keys())
    trackList = copy.deepcopy(adjList)
    endNode = -1
    startNode = -1
    for key in keys:
        indegree = 0
        for otherKey in keys:
            if key in trackList[otherKey]:
                indegree += 1
        outdegree = len(trackList[key])
        if outdegree - indegree == 1:
            startNode = key
        if indegree - outdegree == 1:
            endNode = key
    if startNode == -1:
        return False
    if endNode != -1:
        trackList[endNode].append(startNode)

    for key in keys:
        for item in trackList[key]:
            if item not in trackList:
                endNode = item
                trackList[item] = [startNode]

    currentNode = startNode

    currentcycle = []
    thiscycle = []

    while any(trackList.values()):
        currentcycle = thiscycle

        while trackList[currentNode] != []:
            if currentcycle == []:
                currentcycle.append(currentNode)
            previousNode = currentNode
            choice = randint(0, len(trackList[currentNode]) - 1)
            currentNode = trackList[currentNode][choice]
            currentcycle.append(currentNode)
            trackList[previousNode].remove(trackList[previousNode][choice])
        
        thiscycle = currentcycle 
    
        for key in thiscycle:
            if trackList[key] != []:
                currentNode = key
                n = thiscycle.index(currentNode)
                thiscycle = shift(n, thiscycle[:-1]) + [currentNode]
                break

    thiscycle.pop()
    thiscycle = shift(thiscycle.index(startNode), thiscycle)
    
    return thiscycle

# StringReconstruction
# Input: A collection of k-mers Patterns.
# Output: A string Text with k-mer composition equal to Patterns. 
def StringReconstruction(kmers):
    adjList = deBruijn(kmers)  
    path = eulerianPath(adjList)
    return PathToGenome(path) 

def circularStringReconstruction(kmers):
    adjList = deBruijn(kmers)  
    cycle = eulerianCycle(adjList)
    # remove the last node since it is the same as the first node
    for i in range(len(kmers[0]) - 1):
        cycle.pop()
    return PathToGenome(cycle)

def kUniversalCircularString(k):
    kmers = []
    for i in range(2**k):
        kmers.append(to_binary(i).zfill(k))
    return circularStringReconstruction(kmers)
    

def to_binary(k):
    return bin(k)[2:]

# From here on, the genome assembly method will be focused on generating read-pairs
# A k,d-mer is a pair of k-mers separated by a gap of d nucleotides.
def GeneratekdMers(Genome, k, d):
    kdmers = []
    for i in range(len(Genome)-k-d-k+1):
        kdmers.append((Genome[i:i+k] + "|" + Genome[i+k+d:i+k+d+k]))
    return sorted(kdmers)

def deBruijnFromReadPairs(k, kdmers):
    adjList = {}

    for kdmer in kdmers:
        prefix = kdmer[0:k-1] + "|" + kdmer[k+1:len(kdmer)-1]
        suffix = kdmer[1:k] + "|" + kdmer[k+2:len(kdmer)]
    
        if prefix not in adjList:
            adjList[prefix] = []
        adjList[prefix].append(suffix)
    
    return adjList

# Because of the need to have the same base in a column, not all
# eulerian paths will work for genome assembly. 
def ListToString(list):
    string = ""
    for item in list:
        string += item
    return string

# GraphtoGenomeforkdMers
# Input: An integer k, an integer d, and a graph in the form of adjacency list
# Output: A string Text spelled by a path in this graph, where each node in the path is a k,d-mer.
def GraphtoGenomeforkdMers(k, d, graph):
    potentialPaths = []
    i = 0

    # Creating all the potential paths since not all eulerian paths 
    # will work for genome assembly. This is actualized by repeatedly
    # generating eulerian paths and cycles and checking if they are valid for genome assembly.
    # This algorithm could be improved in the future
    while i < 100:
        Path = eulerianPath(graph)
        if Path == False:
            Path = eulerianCycle(graph)
        if Path not in potentialPaths:
            potentialPaths.append(Path)
        i += 1

    # For each potential path, first create a list of the same length as the genome
    # because lists are more mutable. 
    # Before starting to assemble, check if the path is valid for genome assembly
    # by checking if the bases in the same column of the path are the same. 
    # If not, then this path is not valid for genome assembly and I can skip it.
    # If it is valid, then I can start to assemble the genome by extending the genome one base at a time.
    # For the circular case and additional checking, I check if the bases are consistent
    # once the prefix starts to overlap with the suffix. If not, then this path 
    # is not valid for genome assembly and I can skip it.

    for path in potentialPaths:
        Genome = []
        for i in range (len(path)-2+k+d+k):
            Genome.append(" ")
        validPath = True
        pos = 0
        pairCount = 0
        while pairCount < len(path):
            #Check
            if pairCount > 0:
                for prefixpos1 in range(1, k - 1):
                    if path[pairCount][prefixpos1 - 1] != path[pairCount - 1][prefixpos1]:
                        validPath = False
                        break
                    suffixpos1 = prefixpos1 + k - 1 + 1# position of suffix in path[pairCount]（k-1mer | k-1mer）
                    if path[pairCount][suffixpos1 - 1] != path[pairCount - 1][suffixpos1]:
                        validPath = False                
                        break
                if not validPath:
                    break
            if not validPath:
                break
            #Extend
            for prefixpos in range(k - 1):
                currentBase = Genome[pos + prefixpos]
                Genome[pos + prefixpos] = path[pairCount][prefixpos]
                if currentBase != " " and currentBase != Genome[pos + prefixpos]:
                    validPath = False
                    break
                suffixpos = prefixpos + k + d # position of suffix in Genome (kmer _ _ _ _ (d _'s) kmer)
                Genome[pos + suffixpos] = path[pairCount][prefixpos + k]
                
            pos = pos + 1
            pairCount = pairCount + 1

        if not validPath:
            continue
        break

    return ListToString(Genome)
                
        