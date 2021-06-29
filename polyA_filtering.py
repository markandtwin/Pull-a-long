import pysam, subprocess
import csv, argparse

parser = argparse.ArgumentParser(description='Filter out reads without polyA tail or with misprimed polyA tail')
parser.add_argument('-i','--input', dest='i', action='store',required=True, help='path to input bam file')
parser.add_argument('-o','--output', dest='o', action='store', default='output.bam', help='path to output bam file')
parser.add_argument('-t','--nthread', dest='t', action='store', default='7', help='number of threads')
args = parser.parse_args()


read_lines = []
totalLines = 0
cigar_list = []
sq_list = []

def fetch_data_from_bam(inputFileName):
    global totalLines, cigar_list, sq_list,read_lines
    
    samfile=pysam.AlignmentFile(inputFileName,'rb')
    for read in samfile.fetch():
        totalLines += 1
        if read.query_sequence !=None:
            read_lines.append(read)
            cigar_list.append(read.cigarstring)
            sq_list.append(read.query_sequence)  
    samfile.close()



def valueFetch(aList):
    for item in aList:
        if item.isdigit()== False:
            aList = aList.replace(item, ',')
    aList = aList.split(',') 
    index1 = int(aList[0])
    index2 = int(aList[-2]) 
    return(index1, index2)

def stringSplit(bString, index1, index2):
    subString1 = bString[:index1]
    subString2 = bString[-index2:]
    return(subString1, subString2)

def second_stringSplit(bString, index1, index2):

    subString1 = bString[index1-10:index1+8]
    subString2 = bString[-index2-8:-index2+10]
    return(subString1, subString2)    

def countT(aString):
      num = len(aString)
      if num < 6:
          return False
      elif num > 6 and num <=8:
          if aString.count('T')>=6:
              return True
      else:
          subString = aString[(num-8):]
          if subString.count('T')>=8:
              return True


def first_round_filtering(aList, bList, cList):
    aList_copy = aList.copy()
    bList_copy = bList.copy()
    cList_copy = cList.copy()
    aList.clear()
    bList.clear()
    cList.clear()
    for i in range(totalLines):
        
        (firstIndex, lastIndex) = valueFetch(aList_copy[i])
        preSplit, postSplit = stringSplit(bList_copy[i], firstIndex, lastIndex)
        if countT(preSplit) or countT(postSplit):
            aList.append(aList_copy[i])
            bList.append(bList_copy[i])
            cList.append(cList_copy[i])
    aList_copy.clear()
    bList_copy.clear()
    cList_copy.clear()


def second_round_filtering(aList, bList, cList):
    aList_copy = aList.copy()
    bList_copy = bList.copy()
    cList_copy = cList.copy()
    aList.clear()
    bList.clear()
    cList.clear()

    for i in range(len(bList_copy)):
    
        firstIndex, lastIndex = valueFetch(aList_copy[i])
        preSplit, postSplit = second_stringSplit(bList_copy[i], firstIndex, lastIndex)
        preCondition =  preSplit[:10].count('T')>=8 and preSplit[-8:].count('T')<8
        postCondition =  postSplit[-10:].count('A')>=8 and postSplit[:8].count('A')<8
        if  (preCondition or postCondition):
            if ((preSplit[:10].count('T')>=8) != (postSplit[-10:].count('A')>=8)):
                aList.append(aList_copy[i])
                bList.append(bList_copy[i])
                cList.append(cList_copy[i])

    aList_copy.clear()
    bList_copy.clear()
    cList_copy.clear()



def write_to_text(cList, outfile):
    with open(outfile, 'w') as csvfile:

        for item in cList:
            csvfile.write(item)
            csvfile.write('\n')


def write_to_bamfile(inputFileName, outputFileName, cList):
    samfile = pysam.AlignmentFile(inputFileName, "rb")
    with pysam.AlignmentFile(outputFileName, "wb", template=samfile) as ofile:
        for item in cList:
            ofile.write(item)
    samfile.close()


if __name__=='__main__':
    inputFileName = args.i
    outputFileName = args.o
    fetch_data_from_bam(inputFileName)
  
    write_to_text(sq_list, 'sq.text')
    write_to_text(cigar_list, 'cigar.text')
    #first_round_filtering(cigar_list, sq_list, read_lines)
    second_round_filtering(cigar_list, sq_list, read_lines)
    
    write_to_bamfile(inputFileName, outputFileName, read_lines)
    print('Total lines: %d, valid lines: %d' %(totalLines, len(cigar_list)))
    subprocess.call(['samtools', 'index', '-@', args.t, args.o])
    subprocess.call(['rm', 'sq.text', 'cigar.text'])

