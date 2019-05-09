import numpy as np
#import igraph as ig

# this method create a matrix about score read vs read for all read
def graphMatrix(reads):
    # setting all matrix printable in the screen
    np.set_printoptions(threshold=np.inf)

    # get the number of reads
    nReads = len(reads)

    # instantiate a overlap score matrix  about all reads vs all reads
    graphMat = np.zeros((nReads, nReads))

    # set match/indel/mismatch score
    match = 10
    indel = -4
    mismatch = -1

    # set zero all match/mismatch/indel counter 
    nMatch = 0
    nMismatch = 0
    nIndel = 0

    # initialize a diagonal value for match/mismatch
    diag = 0

    # define two list for creating a align matrix
    alignA = list()
    alignB = list()
    
    # given reads vs reads matrix iterate
    for i in range (0, nReads):
        for j in range (0, nReads):

            # don't compare a read with herself
            if i == j:
                continue
            
            # create a alignment matrix for two reads
            dynamicProgMat = np.zeros((len(reads[i]) + 1, len(reads[j]) + 1))

            # create a direction matrix for rebuild alignment
            # legend: 0 (STOP), 1 (UP), 2 (LEFT), 3 (DIAGONAL)
            directMat = np.zeros((len(reads[i]) + 1, len(reads[j]) + 1))

            # with a couple of read, implement a dynamic programming
            for m in range (0, len(reads[i]) + 1):
                for n in range (0, len(reads[j]) + 1):

                    # base case
                    if m == 0:
                        dynamicProgMat[m][n] = 0
                        continue
                    if n == 0:
                        dynamicProgMat[m][n] = m * indel
                        directMat[m][n] = 1.0
                        continue
                    
                    # general case
                    if reads[i][m - 1] == reads[j][n - 1]:
                        diag = dynamicProgMat[m - 1][n - 1] + match
                    else:
                        diag = dynamicProgMat[m - 1][n - 1] + mismatch  

                    value = max(max(dynamicProgMat[m - 1][n], dynamicProgMat[m][n - 1]) + indel, diag)
                        
                    if value == diag:
                        directMat[m][n] = 3.0
                        dynamicProgMat[m][n] = value 

                    elif value == dynamicProgMat[m][n - 1] + indel:
                        directMat[m][n] = 2.0
                        dynamicProgMat[m][n] = value
                    else:
                        directMat[m][n] = 1.0
                        dynamicProgMat[m][n] = value

            # find max value in the last column (for more details look at algorithm 8)
            m = np.where(dynamicProgMat[:, n] == max(dynamicProgMat[:,n]))[0][0]

            # clear align matrxi
            del alignA[:]
            del alignB[:]

            # set to zero counter for new dynamic program round
            nMatch = 0
            nMismatch = 0
            nIndel = 0
            
            # go back up directMat from max value to first line 
            while m != 0:
                
                if directMat[m][n] == 3.0:
                    m = m - 1
                    n = n - 1
                    if reads[i][m] == reads[j][n]:
                        nMatch += 1
                    else:
                        nMismatch += 1
                    alignA.insert(0, reads[i][m])
                    alignB.insert(0, reads[j][n])

                elif directMat[m][n] == 1.0:
                    m = m - 1
                    alignA.insert(0, reads[i][m])
                    alignB.insert(0, '-')
                    nIndel += 1
                    
                elif directMat[m][n] == 2.0:
                    n = n - 1
                    alignA.insert(0, '-')
                    alignB.insert(0, reads[j][n])
                    nIndel += 1
            
            # delete comment for see the align matrix
            # print(alignA)
            # print(alignB)
            # print("")

            # save with inverse index because in the case of x vs y i have y --> x with overlap
            # (the suffix of y is prefix of x, for more details look at algorithm 8)
            if (nMismatch + nIndel) <= int(0.03*len(reads[j])):
                graphMat[j][i] = nMatch + nMismatch
            else:
                graphMat[j][i] = -1
    return graphMat

# this method create a graph model as (startNode, score, endNode)
def buildGraph(graphMat):

    
    # counter about well knot, we want only one
    cont = 0

    # inizializate graph 
    graph = list()

    # for each row in graphMat get out max overlap score and define and append triple
    for i in range (0, len(graphMat)):
        fromRead = i
        overLap = graphMat.max(axis=1)[i]

        # if is a well knot
        if overLap != 0.0:
            toRead = np.where(graphMat[i, :] == overLap)[0][0]
        else:
            cont += 1
            toRead = -1
        graph.append([fromRead, overLap, toRead])

    # if we have only one well knot return graph otherwise None
    if cont == 1:
        return(graph)
    else:
        return None

# this method check if the graph is available and in the positive case build a genoma and print it and save it in a file
def printFinalGenoma(graph, reads):

    # if we have not a graph to print 
    if graph is None:
        print("With actual quality settings is impossibile to build a genona with overlap graph. There could be more than one end node or no one trim read.")
    else:
        # delete comment to create a visual graph about overlap
        # overlapGraph = ig.Graph(directed = True)
        label = list()
        
        finalGenoma = ""

        # indicies to move in the graph, previus is helpful for overlap graph plot
        toRead = 0
        previous = 0

        # open the file where we will write final genoma
        finalGenomaByOveralp = open("./File/FinalGenomaByOverlap.txt", "w")

        # counter for hamilton path check
        numberOfSteps = 1

        # pick-up firt element for search the real starter
        starter = graph[0][0]

        # search the real starter that is the node that show in "from" but never show in "to"
        i = 0
        while i < len(graph):
            if starter == graph[i][2]:
                starter = graph[i][0]
                i = -1
            i += 1

        # delete comment to add vertex to graph 
        # overlapGraph.add_vertex(starter)
        # overlapGraph.add_vertex(toRead)

        # delete comment to add the read to label for graph plot
        # label.append(reads[starter])

        # add to final genoma the starter read 
        finalGenoma = finalGenoma + reads[starter]

        # indicate the next read to read and the overlap
        toRead = graph[starter][2]
        overlap = graph[starter][1]

        # delete comment to add edge
        # overlapGraph.add_edges([(starter, toRead)])

        # iterate the graph and get the final genoma
        while overlap != 0.0:

            # count number of step for hamilton path
            numberOfSteps += 1
            finalGenoma = finalGenoma + reads[toRead][int(overlap):]

            # label.append(reads[toRead][int(overlap):])

            previous = toRead
            overlap = graph[toRead][1]
            toRead = graph[toRead][2]
            # if toRead != -1:
                #overlapGraph.add_vertex(toRead)
                #overlapGraph.add_edges([(previous, toRead)])

        # check about hamilton path
        if numberOfSteps != len(graph):
            print("Is not possible done hamilton path")
            return None

        

        print(finalGenoma)

        # write final genoma to file
        finalGenomaByOveralp.write(finalGenoma)

        # delete comment to plot the overlap graph
        #ig.plot(overlapGraph, vertex_label = label, bbox=(0, 0, 1800, 1800), margin = (300, 300, 300, 300))

