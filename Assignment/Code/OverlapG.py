import overlapGraphOperation as ogo

def overlapGraph():
    # open file with clean read
    cleanReadsNGS = open("./File/CleanReadsNGS.txt", "r")
    cleanReads = list()

    # save in a list the clean read
    for line in cleanReadsNGS:
        cleanReads.append(line[:len(line)-1])

    # start overlap graph building
    graphMatrix = ogo.graphMatrix(cleanReads)
    graph = ogo.buildGraph(graphMatrix)
    ogo.printFinalGenoma(graph, cleanReads)
    cleanReadsNGS.close()
