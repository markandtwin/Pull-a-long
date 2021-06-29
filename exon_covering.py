import pysam, subprocess, argparse
from csv import writer, reader


pre_index = []
post_index = []
counter = 0
sum = []
pre_delta = []
post_delta = []
deal_all_valid = []
deal_all_invalid = []
parser = argparse.ArgumentParser(description='Filter out reads covering less than 80% of given region.')
parser.add_argument('-i','--input', dest='i', action='store',required=True, help='path to input bam file')
parser.add_argument('-o','--output', dest='o', action='store', default='output.bam', help='path to output bam file')
parser.add_argument('-r','--rest', dest='r', action='store', default='rest.bam', help='path to failed reads bam file')
parser.add_argument('-b','--bed', dest='b', action='store', required=True, help='bed4 file of filtering regions')
parser.add_argument('-t','--nthread', dest='t', action='store', default='7', help='number of threads')
parser.add_argument('-p','--percentage', dest='p', action='store', default= 0.8, help='coverage ratio')
args = parser.parse_args()


def read_file_to_csv(cname,firstindex,lastindex,inputFileName):
    global pre_index, post_index, counter, post_delta, post_delta
    samfile=pysam.AlignmentFile(inputFileName,'rb')
    items = samfile.fetch(cname,firstindex, lastindex)

    
    with open('data.csv', 'w') as infile:
        write_to_csv = writer(infile)
        for read in items:
            counter += 1
            write_to_csv.writerow(read.get_blocks())
            for i in range(len(read.get_blocks())):
                if read.get_blocks()[len(read.get_blocks())-1][0] < firstindex:
                    pre_delta.append(0)
                    pre_index.append(len(read.get_blocks())-1)
                    break 
                else:
                    if read.get_blocks()[i][0] > firstindex or read.get_blocks()[i][0] == firstindex:
                        if read.get_blocks()[i-1][1] > firstindex and i > 0:
                            delta_front = read.get_blocks()[i-1][1]-firstindex
                        else:
                            delta_front = 0
                        pre_delta.append(delta_front)
                        pre_index.append(i)
                        break
            for j in range(len(read.get_blocks())):
                if read.get_blocks()[len(read.get_blocks())-1][1] < lastindex:
                    post_delta.append(0)
                    post_index.append(len(read.get_blocks()))
                    break
                else:
                    if read.get_blocks()[j][1] > lastindex or read.get_blocks()[j][1] == lastindex:
                        if read.get_blocks()[j][0] < lastindex:
                            delta_back = lastindex - read.get_blocks()[j][0]
                        else:
                            delta_back = 0
                        post_delta.append(delta_back)
                        post_index.append(j)
                        break
    samfile.close()



def deal_file(cname,firstindex,lastindex,inputFileName):
    reads = []
    deal =[]
 
    samfile=pysam.AlignmentFile(inputFileName,'rb')
    items = samfile.fetch(cname,firstindex, lastindex)
    for read in items:
        reads.append(read.get_blocks())

    with open('data2.csv', 'w') as valuefile:
        write_to_csv = writer(valuefile)
        for read in reads:
            for i in read:
                deal.append(i[1]-i[0])
            write_to_csv.writerow(deal)
            deal.clear()

def calculate():
    global pre_index, post_index, sum, counter, pre_delta, post_delta
    sum_list =[]
    for i in range(counter):
        sum.append(0)
    t =0
        
    with open("data2.csv") as csvfile:
        reader_csv = reader(csvfile)
        k = 0
        for row in reader_csv:
            sum_list.append((row[pre_index[k]: post_index[k]]))
            k = k+1
        for ros in sum_list:
            for j in ros:
                sum[t] += int(j)
            t = t + 1
           
def read_bam_file(cname,firstindex,lastindex,inputFileName, outputFileName_valid, outputFileName_invalid, gname):
    global sum, counter, pre_delta, post_delta, deal_all_valid, deal_all_invalid
    deal =[]
    invalid_deal = []

    invalidNums = 0
    invalidLine = []
    for i in range(counter):
        rate = float(pre_delta[i]+post_delta[i]+sum[i])/(lastindex - firstindex)
        if rate < float(args.p):
            invalidNums += 1 
            invalidLine.append(i)
    print(gname, "valid lines: ", (counter - invalidNums))
    print(gname, "total lines: ", counter)
 
    samfile=pysam.AlignmentFile(inputFileName,'rb')
    items = samfile.fetch(cname,firstindex, lastindex)
    for read in items:
        deal.append(read)

    
    for i in range(invalidNums):
        invalid_deal.append(deal[invalidLine[i]])

    for line_invalid in invalid_deal:
        deal_all_invalid.append(line_invalid)
    #write each section into a file
    # with pysam.AlignmentFile(outputFileName_invalid, "wb", template=samfile) as outfile:
    #     for ips in invalid_deal:
    #         outfile.write(ips)

    ind = 0
    for i in range(invalidNums):
        deal.remove(deal[invalidLine[i] - ind])
        ind =  ind +1
    
    #append valid lines into the list of deal_all_valid
    for line_valid in deal:
        deal_all_valid.append(line_valid)


    # with pysam.AlignmentFile(outputFileName_valid, "wb", template=samfile) as outf:
    #     for ps in deal:
    #         outf.write(ps)
    samfile.close()


def composed_func(cname,firstindex,lastindex, gname,inputFileName, outputFileName_valid, outputFileName_invalid):
    global sum, counter, pre_delta, post_delta
    deal_file(cname,firstindex,lastindex,inputFileName)
    read_file_to_csv(cname,firstindex,lastindex,inputFileName)
    calculate()
    read_bam_file(cname,firstindex,lastindex,inputFileName, outputFileName_valid, outputFileName_invalid, gname)
    pre_index.clear()
    post_index.clear()
    counter =0
    sum.clear()
    pre_delta.clear()
    post_delta.clear()



def read_bed_file(inputFileName, outputFileName_valid, outputFileName_invalid, bedFileName):
    fileNum = 0
    sp1 = outputFileName_valid.split(".")
    sp2 = outputFileName_invalid.split(".")
    
    with open(bedFileName, 'r') as bedfile:
        lines = bedfile.readlines()

    for line in lines:
        chname, first, final, name = line.split()
        first = int(first)
        final = int(final)
        outputFileName_valid_num = sp1[0] + "_" + str(fileNum)+ "." + sp1[1]
        outputFileName_invalid_num = sp2[0] + "_" + str(fileNum)+ "." + sp2[1]
        composed_func(chname, first, final, name, inputFileName, outputFileName_valid_num, outputFileName_invalid_num)
        fileNum += 1

        

def output_composed_bam(inputFileName,outputFileName_valid, outputFileName_invalid):
    global deal_all_invalid, deal_all_valid
    samfile=pysam.AlignmentFile(inputFileName,'rb')
    with pysam.AlignmentFile("pass.bam", "wb", template=samfile) as oct1file:
        for ps in deal_all_valid:
            oct1file.write(ps)
    subprocess.call(['samtools','sort','-@', args.t,'pass.bam', '-o', outputFileName_valid])
    subprocess.call(['samtools', 'index', '-@', args.t, outputFileName_valid])

    with pysam.AlignmentFile("fail.bam", "wb", template=samfile) as oct2file:
        for ips in deal_all_invalid:
            oct2file.write(ips)
    subprocess.call(['samtools','sort','-@', args.t,'fail.bam', '-o', outputFileName_invalid])
    subprocess.call(['samtools', 'index', '-@', args.t, outputFileName_invalid])
    subprocess.call(['rm', '-f', 'data.csv', 'data2.csv', 'pass.bam', 'fail.bam'])
    samfile.close()
        


if __name__=="__main__":

    iFileName = args.i
    oFileName_valid = args.o
    oFileName_invalid = args.r
    bFileName = args.b
    read_bed_file(iFileName, oFileName_valid, oFileName_invalid, bFileName)
    output_composed_bam(iFileName,oFileName_valid, oFileName_invalid)
 
