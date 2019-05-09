from DebrujinG import debrujin
from trimming import trimmer
from OverlapG import overlapGraph
from datetime import datetime
from memory_control import measure_memory_usage

def main():

    # begin trimmer with user input
    # trimmer()
    
    # begin overlap block
    start1 = datetime.now()

    # call to overlapGraph function and memory check
    max_mem1 = measure_memory_usage(overlapGraph, "", memory_denominator=1024.0, memory_usage_refresh=0.005)
    print ("{:>32} {:>14}".format("Memory used for overlap step:",str(max_mem1) + " KB"))
    end1 = datetime.now() - start1
    print ("{:>32} {:>14}".format("Execution time for overlap step:", str(end1)))
    
    # begin debrujin block
    start2 = datetime.now()

    # call to debrujin function and memory check
    max_mem2 = measure_memory_usage(debrujin, "", memory_denominator=1024.0, memory_usage_refresh=0.005)
    print ("{:>33} {:>14}".format("Memory used for debrujin step:",str(max_mem2) + " KB"))
    end2 = datetime.now() - start2
    print ("{:>33} {:>14}".format("Execution time for debrujin step:", str(end2)))
    

# this means that if this script is executed, then 
# main() will be executed
if __name__ == '__main__':
    main()
