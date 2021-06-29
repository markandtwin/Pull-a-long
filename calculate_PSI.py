import pysam,  argparse
from csv import writer, reader

parser = argparse.ArgumentParser(description='Filter out reads without polyA tail or misprimed')
parser.add_argument('-i','--input', dest='i', action='store',required=True, help='path to input bam file')
parser.add_argument('-o','--output', dest='o', action='store', default='PSI.csv', help='path to output csv file')
parser.add_argument('-b','--bed', dest='b', action='store', required=True, help='bed4 file of cassette exon')
parser.add_argument('-t','--nthread', dest='t', action='store', default='7', help='number of threads')
parser.add_argument('-p','--percentage', dest='p', action='store', default= 0.8, help='coverage ratio')
args = parser.parse_args()

total_blocks_list = []

def read_data_from_bedFile(bedFileName):
    ch_index_lists = []
    with open(bedFileName, 'r') as bedfile:
        lines = bedfile.readlines()
        for line in lines:
            ch_index_lists.append(line.split())
    return ch_index_lists


def fetch_blocks_from_bamfile(chname,firstindex,lastindex,inputFileName):
    global total_blocks_list
    firstindex = firstindex - 50
    lastindex = lastindex + 50
    samfile=pysam.AlignmentFile(inputFileName,'rb')
    items = samfile.fetch(chname,firstindex, lastindex)
    for read in items:
        total_blocks_list.append(read.get_blocks())
    samfile.close()



def fetch_overlap_value(blocks_lists, firstindex, lastindex):
    pre_index = []
    post_index = [] 
    pre_delta = [] 
    post_delta = []

    for blocks_row in blocks_lists:
        for i in range(len(blocks_row)):
            if blocks_row[len(blocks_row)-1][0] < firstindex:
                pre_delta.append(0)
                pre_index.append(len(blocks_row)-1)
                break 
            else:
                if blocks_row[i][0] > firstindex or blocks_row[i][0] == firstindex:
                    if blocks_row[i-1][1] > firstindex and i > 0:
                        delta_front = blocks_row[i-1][1]-firstindex
                    else:
                        delta_front = 0
                    pre_delta.append(delta_front)
                    pre_index.append(i)
                    break
        for j in range(len(blocks_row)):
            if blocks_row[len(blocks_row)-1][1] < lastindex:
                post_delta.append(0)
                post_index.append(len(blocks_row))
                break
            else:
                if blocks_row[j][1] > lastindex or blocks_row[j][1] == lastindex:
                    if blocks_row[j][0] < lastindex:
                        delta_back = lastindex - blocks_row[j][0]
                    else:
                        delta_back = 0
                    post_delta.append(delta_back)
                    post_index.append(j)
                    break
    return pre_index, post_index, pre_delta, post_delta


def calculate_difference_of_block(blocks_lists):

    diff_row = []
    diff_lists = []

    for blocks_row in blocks_lists:
        for block in blocks_row:
            diff_row.append(block[1]-block[0])
        diff_lists.append(diff_row)
    return diff_lists


def calculate_sum_of_nonoverlap(individual_value_list, pre_index_list, post_index_list):

    sum_value_list = []

    for individual_value_row, i in zip(individual_value_list, range(len(pre_index_list))):
        sub_value_row = individual_value_row[pre_index_list[i]: post_index_list[i]]
        sub_value_sum = sum(sub_value_row)
        sum_value_list.append(sub_value_sum)
    return sum_value_list


def first_round_filtering(first, last):
    global total_blocks_list
    total_blocks_list_copy = total_blocks_list.copy()
    total_blocks_list.clear()
    pre_index_l, post_index_l, pre_delta_l, post_delta_l = fetch_overlap_value(total_blocks_list_copy, first-50, first-1)
    pre_index_r, post_index_r, pre_delta_r, post_delta_r = fetch_overlap_value(total_blocks_list_copy, last+1, last+50)
    diff_list = calculate_difference_of_block(total_blocks_list_copy)
    sum_list_l = calculate_sum_of_nonoverlap(diff_list, pre_index_l, post_index_l)
    sum_list_r = calculate_sum_of_nonoverlap(diff_list, pre_index_r, post_index_r)

    for i in range(len(total_blocks_list_copy)):
        pre_rate = pre_delta_l[i]+post_delta_l[i]+sum_list_l[i]/50.0
        post_rate = pre_delta_r[i]+post_delta_r[i]+sum_list_r[i]/50.0
        if pre_rate < 0.1 and post_rate < 0.1:
            total_blocks_list.append(total_blocks_list_copy[i])


def second_round_filtering(first, last):
    global total_blocks_list
    total_blocks_list_copy = total_blocks_list.copy()
    total_blocks_list.clear()
    pre_index, post_index, pre_delta, post_delta = fetch_overlap_value(total_blocks_list_copy, first, last)
    diff_list = calculate_difference_of_block(total_blocks_list_copy)
    sum_list= calculate_sum_of_nonoverlap(diff_list, pre_index, post_index)

    a = 0
    b = 0

    for i in range(len(total_blocks_list_copy)):
        rate = pre_delta[i]+post_delta[i]+sum_list[i]/float(last-first)
        if rate > float(args.p):
            a += 1
        elif rate < 0.1:
            b  += 1
    return a ,b


def calculate_psi(a, b):
    if (a+b) !=0:
        psi = float(a)/(a+b)
        return psi


def create_csv_for_result(outputResultFile):
    with open (outputResultFile, 'w') as firstfile:
        pass

def append_result_to_csv(totalLines, validLines, psi, resultFile, gene):
    with open(resultFile, 'a') as outfile:
        write_csv = writer(outfile, delimiter=',')
        write_csv.writerow([gene, validLines, totalLines, psi])





if __name__=="__main__":

    bedFileName = args.b
    inputFileName = args.i
    outputResultFile = args.o

    create_csv_for_result(outputResultFile)

    ch_index_lists = read_data_from_bedFile(bedFileName)

    for ch_index_list in ch_index_lists:
        cname = ch_index_list[0]
        firstindex = int(ch_index_list[1])
        lastindex = int(ch_index_list[2])
        gene = ch_index_list[3]

        fetch_blocks_from_bamfile(cname,firstindex,lastindex,inputFileName)
        first_round_filtering(firstindex, lastindex)
        a, b = second_round_filtering(firstindex, lastindex)
        psi = calculate_psi(a, b)
        total_blocks_list.clear()
        append_result_to_csv((a+b), a, psi, outputResultFile, gene)


   
 
