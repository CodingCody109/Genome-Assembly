# get_kmers
# Input: A string Genome and an integer k.
# Output: The set of all k-mers in Genome.
def get_kmers(k, Genome):
    kmers = []
    for i in range(len(Genome)-k+1):
        kmers.append(Genome[i:i+k])
    return kmers   
   
def PathToGenome(Path):
    Genome = Path[0]
    for i in range(1, len(Path)):
        Genome += Path[i][len(Path[i])-1]
    return Genome

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

def deBruijn(kmers):
    adjList = {}
    
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        
        if prefix not in adjList:
            adjList[prefix] = []
        adjList[prefix].append(suffix)  # append, don't overwrite
    
    order = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return dict(sorted(adjList.items(), key=lambda x: [order[c] for c in x[0]]))