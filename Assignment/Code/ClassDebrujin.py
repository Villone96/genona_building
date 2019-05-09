# class edge of the DeBrujin graph
class Edge:
    def __init__(self, lab):
        self.label = lab

# class vertex of the DeBrujin graph 
class Vertex:
    def __init__(self, lab):
        self.indegree = 0
        self.outdegree = 0
        self.label = lab
        self.occurence = 1

# function to read all the cleaned reads from the correctr file. It is assumed that one line in the input file correspond 
# to one read.
def readReads(fileReads):
    # open the relative file in read mode
    cleanReadsNGS = open(fileReads, "r")
    cleanReads = list()

    # for every row read from the file
    for read in cleanReadsNGS:
        # append it into my structure
        cleanReads.append(read[:len(read)-1])   
     
    return cleanReads

# function that define a De Brujin graph given all reads and k which represent the length of every k-mer to extract from reads
def buildGraph(reads, k):

    # dictionary for graph's vertices
    vertex = dict()

    # dictionary for graph's vertices
    edge = dict()

    # variable to take note if insert a node or not
    insert = True

    # for every read
    for read in reads:

        # index used to create k-mers
        i = 0

        # (i + k) will be the last caracter of the read
        while i + k < len(read):
            
            # extract the [i:i+k] k-mer
            vertex_1 = read[i:i+k]

            # extract the [i+1:i+k+1] k-mer            
            vertex_2 = read[i+1:i+k+1]

            # if the [i:i+k] k-mer is already a node into the graph 
            if vertex_1 in vertex.keys():
                # for every edge of the [i:i+k] k-mer (seen as vertex reachable from [i:i+k] node )
                for e in edge[vertex_1]:

                    # if vertex2 already exists
                    if e.label == vertex_2:
                        # don't insert it again
                        insert = False

                # if vertex2 not already exists
                if insert:
                    # add the relative edge between the [i:i+k] k-mer and the [i+1:i+k+1] k-mer
                    edge[vertex_1] += [Edge(vertex_2)]

                    # increment vertex1 outdegree   
                    vertex[vertex_1].outdegree += 1
                    
                # increment vertex_1 occurrence
                vertex[vertex_1].occurence += 1
            else:
                # otherwise create the vertex for the [i:i+k] k-mer
                vertex[vertex_1] = Vertex(vertex_1)

                # increment its outdegree   
                vertex[vertex_1].outdegree += 1

                # add the relative edge between the [i:i+k] k-mer and the [i+1:i+k+1] k-mer
                edge[vertex_1] = [Edge(vertex_2)]

                # enable the vertex2 indegree increase
                insert = True


            # if the [i+1:i+k+1] k-mer is already a node into the graph 
            if vertex_2 in vertex.keys():

                # if the increase was enabled
                if insert:
                    vertex[vertex_2].indegree += 1

                # increment vertex_2 occurrence
                vertex[vertex_2].occurence += 1

            else:
                # otherwise create the vertex for the [i+1:i+k+1] k-mer
                vertex[vertex_2] = Vertex(vertex_2)

                # increment it's indegree by 1 (because we added an incoming edge coming from vertex1)
                vertex[vertex_2].indegree += 1

                # set outcoming edges from the [i+1:i+k+1] k-mer as an empty list
                edge[vertex_2] = []

            # next k-mers of the same read
            i += 1

            # set default flag value
            insert = True


    return (vertex, edge)   

# function that performs visit in the graph to output genome assembly
def printFinalGenoma(debrujin_g):

    # extract both vertices and edges dicts from the graph
    vertex = debrujin_g[0]
    edge = debrujin_g[1]

    # counter about number of nodes unbalanced in the graph
    unbalanced = 0

    # counter about number of nodes that have 0 out degree (-1 because i accept one aka my real end)
    outDegreeZero = -1

    # counter about number of nodes that have 0 in degree (-1 because i accept one aka my real start)
    inDegreeZero = -1

    
    # if the graph has no vertices
    if len(vertex) == 0:
        print("It's not possible to build De Brujin Graph")
        return None

    # for every vertex
    for v in vertex:
        # count and print information about the number of vertexs that have not balanced degree and not 0 in in/out degree
        if vertex[v].indegree != vertex[v].outdegree and vertex[v].outdegree != 0 and vertex[v].indegree != 0:
            unbalanced += 1
            # print(vertex[v].label)
            # print(vertex[v].indegree)
            # print(vertex[v].outdegree)
            # print("")
    
        # count and print information about the number of vertexs that have 0 in indegree, taken out my real start
        if vertex[v].indegree == 0:
            inDegreeZero += 1
            # if inDegreeZero > 0:
            #    print(vertex[v].label)
            #    print(vertex[v].indegree)
            #    print(vertex[v].outdegree)
            #    print("")

        # count and print information about the number of vertexs that have 0 in out degree, taken out my real end
        if vertex[v].outdegree == 0:
            outDegreeZero += 1
            # if outDegreeZero > 0:
            #    print(vertex[v].label)
            #    print(vertex[v].indegree)
            #    print(vertex[v].outdegree)
            #    print("")

    ### WE ADDED COMMENT TO ABOVE PRINT FOR FOR MORE VISUAL CLARITY DURING THE OUTPUT, HOWEVER IT ARE USEFULL FOR UNDERSTANDING 
    ### ERROR STRUCTURE IN THE GRAPH. DELETE IT FOR SE THE START/END BUBBLES NODES AND TIPS FORK


    # output about all counter 
    print("")
    print("{:>30} {:>2} ".format("unbalanced before subtraction:", unbalanced))
    # for this substraction see documentation to section 5.8
    unbalanced -= (outDegreeZero + inDegreeZero)
    print("{:>30} {:>2} ".format("unbalanced after subtraction:", unbalanced))
    print("{:>30} {:>2} ".format("indegree equals to zero:", inDegreeZero))
    print("{:>30} {:>2} ".format("outdegree equals zero:", outDegreeZero))
    
    


    # if there is a bubbles
    if unbalanced > 0:
        print("There are {:>2} bubbles in the Bruijin graph".format(int((unbalanced)/2)))

    # if there is a tips
    if outDegreeZero + inDegreeZero > 0:
        print("There are {:>2} tips in the Bruijin graph".format((outDegreeZero + inDegreeZero)))

    # Pick first node in the vertex 
    start = list(vertex.keys())[0]
    
    # for every vertex
    for key in vertex.keys():
        # if its outdegree is not zero and its indegree its lower than the actual starting node
        if vertex[key].outdegree != 0 and vertex[key].indegree < vertex[start].indegree:
            # set the actual starting node as the actual "explored" vertex
            start = key

    # the final sequence will start with the starting node label
    contig = start

    # obviously at the first round the starting node will be the current node
    current = start

    # Pick the next node that will become the current node as the first not already visited into current node adjacency list
    # while the adjacency list of the current node is not empty
    while len(edge[current]) > 0:

        # counter for performing max adjacency list node picking 
        maxOccurrence = 0
        maxIndex = 0

        # for every edge outcoming from current node
        for e in edge[current]:
            # find the one with max occurences (it's likely to be correct)
            if vertex[e.label].occurence > maxOccurrence:
                maxOccurrence = vertex[e.label].occurence
                maxIndex = edge[current].index(e) 
        
        # set the next node as the first not already visited into current node adjacency list
        next = edge[current][maxIndex]
        
        # permanently remove that edge from the adjacency list of the current node
        del edge[current][maxIndex]

        # append the last label caracter to extend the sequence
        contig += next.label[-1]

        # Update the current node for the next iteration
        current = next.label
    
    # return the sequence
    return contig