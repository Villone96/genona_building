#from tkinter import Tk
#from tkinter.filedialog import askopenfilename
import readsOperation as ro

def trimmer():

    # In the final release use this code for choose any file
    # Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    # filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
    # print(filename)
    # readsNGS = open(filename, "r")
    
    # in the case the above code doesn't work for compatibility problems, we open this a ./File/ReadsNGS.txt with a Reads in 
    # FASTQ format
    readsNGS = open("./File/ReadsNGS.txt", "r")

    # open a file for cleanRead
    cleanReadsNGS = open("./File/CleanReadsNGS.txt", "w")

    # delete comments for allow the user input
    stay = True

    while stay:
        try:
            stay = False
            qualityValue = int(input("Enter the required quality value: "))
            if qualityValue > 93:
                qualityValue = 93
            lenghtValue = int(input("Enter the required lenght value: "))
        except:
            stay = True
            print("Bad input format, try again")

    read = ""
    qualityString = ""

    # cont for indicate which read doesn't respect user input (in the case of it)
    i = 1

    # if comments about user input was delete, comment next two line
    # qualityValue = 50
    # lenghtValue = 1

    # for all line in the FASTQ file
    for line in readsNGS:
        # take read
        if line.startswith('@'):
            read = next(readsNGS)
            number = line[5:7]
        # take quality read line
        if line.startswith('+'):
            qualityString = next(readsNGS)
            
            # call function for trimming and retrieve, if possible, the triple
            index = ro.trimOperation(qualityString, qualityValue, lenghtValue)

            # if trimming is not possible, print a message about it and the number of read that caused it 
            # after that stop esecution of trimming
            if index is None:
                print("The minimum lenght is not respected by read number: " +number)
                break
            i += 1
        
            # write in clean read file the clean read
            cleanReadsNGS.write(read[index[0]:index[1]])
            if len(read[index[0]:index[1]]) != len(qualityString):
                cleanReadsNGS.write("\n")

    readsNGS.close()
    cleanReadsNGS.close()
    

