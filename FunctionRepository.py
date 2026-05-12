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
    
    order = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return dict(sorted(adjList.items(), key=lambda x: [order[c] for c in x[0]]))

def eulerianPath(adjList):
    keys = list(adjList.keys())
    trackList = adjList.copy()
    startNode = keys[randint(0, len(keys)-1)]
    currentNode = startNode
    path = []

    while trackList != {}:
        while trackList[startNode] != []:
            currentNode = trackList[currentNode][0]
            trackList[currentNode].remove(trackList[currentNode][0])
    
        while trackList[startNode] == []:
            for key in trackList:
                if trackList[key] != []:
                    startNode = key
                    break
            break
    
    currentNode = startNode
    while trackList != {}:
        while trackList[currentNode] != []:
            path.append(currentNode)
            currentNode = trackList[currentNode][0]
            trackList[currentNode].remove(trackList[currentNode][0])
    
        while trackList[currentNode] == []:
            for key in trackList:
                if trackList[key] != []:
                    currentNode = key
                    break
            break
    
    return path
        

    
