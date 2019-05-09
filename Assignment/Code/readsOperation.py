from operator import itemgetter

# core function, perform trimming operation on qualityString in function of qualityValue e lengthValue 
def trimOperation(qualityString, qualityValue, lengthValue):  

    # initialization of list for all possible substring interval
    subString = []

    # start index
    start = 0

    # end index
    end = 0

    # flag about possibility of all intervall return (from first to end of qualityString) 
    flag = False
    
    for i in range(0, len(qualityString)-1):
        # if the char phared value is < of quality user value
        if convertASCIIToDecimal(qualityString[i]) < qualityValue:
            # if the indexes are overlap 
            if end - start == 0:
                end = i + 1
                start = i + 1
            else:
                # otherwise append the triple
                subString.append([start, end , end - start])
                end = i + 1
                start = i + 1
        else:
            flag = True
            end += 1

    # if the flag is yet false means that no one char has a phred value > of quality
    if not flag:
        return None

    # if the substring is empy but flag is true (all chars have a phred value > of qualityValue)
    if subString == [] and flag:
        # if the qualitySatring is less long of lenghtValue is impossible done trimming
        if len(qualityString) - 1 < lengthValue:
            return None
        else:
            # otherwise return from first index to last 
    	    return [0, len(qualityString), len(qualityString) - 1]
    else:
        # in the case of different intervals in the quality String take longest one and, if it is longer than length Value, return
        # otherwise is impossible done trimming
        LongSubString = sorted(subString, key = itemgetter(2))[-1]
        if LongSubString[1] - LongSubString[0] < lengthValue:
            return None
        else:
            return LongSubString

# this function convert Ascii to Decimal for above check 
def convertASCIIToDecimal(ASCIILetter):
	decimalValue = ord(ASCIILetter) - 33
	return decimalValue