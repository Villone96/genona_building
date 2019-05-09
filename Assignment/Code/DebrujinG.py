import readsOperation as ro
import ClassDebrujin as dbgo

def debrujin(): 
    # length of k-mers
    k = 31

    # read clean reads from file
    reads = dbgo.readReads('./File/CleanReadsNGS.txt')

    # build graph given relatives parameters
    g = dbgo.buildGraph(reads, k)

    # print final genoma from graph, return it and save on the file 
    contig = dbgo.printFinalGenoma(g)

    if contig != None:
        print(contig)

        with open("./File/FinalGenomaByDeBrujin.txt", "w") as f:
            f.write(contig)
                                         
