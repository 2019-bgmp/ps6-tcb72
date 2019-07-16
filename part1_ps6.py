###***UNIT TEST MANUAL RESULTS***###

#SEE: Unit_test.fa

#N50
#sorted sequence (descending) below
#TTCGACGTTCGACGTTCGACGTTCGACGTTCGACGTTCGACGTTCGACGTTCGACGTTCGACGTTCGACGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGCTGAGCTGAGCTGAGCTGAGCTGAGCTGAGCTGAGCTGAGCTGAGCTGAGCT
# total length of 182
# N50 is at 91st nucleotide
# 91st nucleotide is in the second sequence
# length of second sequence is 57
# therefore, N50 is 57.

#maximum contig length
#max of 70, 55, 57, so max is 70.

#mean contig length
#mean is (70+55+57) / 3, which is approximately 60.67

#total length
#total length is 70 + 55 + 57, which is 182.

#number of contigs
#number of contigs is 3

#mean depth of coverage
#coverage for single contig is Ck * L / (L - Lk + 1) where L is physical length
#first: L is 70, Ck is 19.4929384, Lk is 49. Plug that in, C = 62.02
#second: L is 55, Ck is 52.4303022, Lk is 49. Plug that in, C = 411.95
#third: L is 57, Ck is 43.1929483, Lk is 49. Plug that in, C = 273.55
#average: (62.02+411.95+273.55)/3 = 249.173
###************************************###

import re

def get_physical_lengths(filename, KMER_LENGTH):
    physical_lengths = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                kmer_length_header = int(re.findall('length_[0-9]+', line)[0][7:])
                kmer_coverage = float(re.findall('cov_[0-9]+\.[0-9]+', line)[0][4:])
                physical_lengths.append(kmer_length_header + KMER_LENGTH - 1)
    return(physical_lengths)

def get_kmer_coverages(filename):
    kmer_coverages = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                kmer_coverages.append(float(re.findall('cov_[0-9]+\.[0-9]+', line)[0][4:]))
    return(kmer_coverages)

def get_num_of_contigs(filename):
    count = 0
    with open(filename,'r') as f:
        for line in f:
            if line.startswith('>'):
                count+=1
    return(count)

def get_max_contig_length(filename):
    seq_dict = {}
    max_length = 0
    with open(filename,'r') as f:
        for index, line in enumerate(f):
            line = line.strip()
            if line.startswith('>'):
                current_key = line
                seq_dict[current_key] = ''
            else:
                seq_dict[current_key] += line
    for key in seq_dict:
        current_value = len(seq_dict[key])
        if current_value > max_length:
            max_length = current_value
    return(max_length)

def get_total_length(filename):
    nuc_sum = 0
    with open(filename,'r') as f:
        for index, line in enumerate(f):
            line = line.strip()
            if not line.startswith('>'):
                for char in line:
                    nuc_sum += 1
    return(nuc_sum)


def get_mean_depth_coverage(filename, kmer_coverages, physical_lengths, KMER_LENGTH):
    cov_sum = 0
    for index, length in enumerate(physical_lengths):
        current_cov = (kmer_coverages[index] * length) / (length - KMER_LENGTH + 1)
        cov_sum += current_cov
    return(cov_sum/len(physical_lengths))

def get_n50(filename, total_length):
    max_length = 0
    nuc_count = 0
    seq_dict = {}
    with open(filename,'r') as f:
        for index, line in enumerate(f):
            line = line.strip()
            if line.startswith('>'):
                current_key = line
                seq_dict[current_key] = ''
            else:
                seq_dict[current_key] += line
        all_contigs = list(seq_dict.values())
        all_contigs.sort(key=len,reverse=True)
        perc_50 = total_length/2
        for contig in all_contigs:
            for char in contig:
                nuc_count += 1
                if nuc_count >= perc_50:
                    return(len(contig))

def get_contig_dist(filename):
    dist = {}
    seq_dict = {}
    with open(filename, 'r') as f:
        for index,line in enumerate(f):
            line = line.strip()
            if line.startswith('>'):
                current_key = line
                seq_dict[current_key] = ''
            else:
                seq_dict[current_key] += line
        for key in seq_dict:
            current_sequence = seq_dict[key]
            current_length = len(current_sequence)
            lowest_100 = int(current_length/100) * 100
            if lowest_100 in dist:
                dist[lowest_100] += 1
            else:
                dist[lowest_100] = 1
    print('Contig Length' + '\t' + 'Number of contigs in this category')
    for key in sorted(dist.keys()):
        print(str(key) + '\t' + str(dist[key]))
    return(dist)



def main():
    filename = '/gpfs/projects/bgmp/tbiondi/k49/contigs.fa'
    KMER_LENGTH = 49
    num_contigs = get_num_of_contigs(filename)
    print('\n********************')
    print('Processed ' + filename + '\n')
    print('There are ' + str(num_contigs) + ' contigs in this file.' )
    max_contig_length = get_max_contig_length(filename  )
    print('The maximum contig length is ' + str(max_contig_length) + ' nt.')
    total_length = get_total_length(filename)
    print('The total length of this file is ' + str(total_length) + ' nt.')
    mean_contig_length = total_length / num_contigs
    print('The mean contig length is approximately ' + str(round(mean_contig_length,2)) + ' nt.')
    kmer_coverages = get_kmer_coverages(filename)
    physical_lengths = get_physical_lengths(filename, KMER_LENGTH)
    mean_depth_coverage = get_mean_depth_coverage(filename, kmer_coverages, physical_lengths, KMER_LENGTH)
    print('The mean depth of coverage is approximately ' + str(round(mean_depth_coverage,2)) + 'x.')
    n50 = get_n50(filename, total_length)
    print('The N50 for this dataset is ' + str(n50) + ' nt.')
    print('********************')
    print('\n***CONTIG DISTRIBUTION***')
    contig_dist = get_contig_dist(filename)

if __name__ == "__main__":
    main()
