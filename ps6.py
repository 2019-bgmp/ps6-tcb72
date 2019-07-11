import re

def get_physical_lengths(filename, KMER_LENGTH):
    physical_lengths = []
    with open(filename, 'r') as f:
        for index,line in enumerate(f):
            if line[0] == ">":
                kmer_length_header = int(re.findall('length_[0-9]+', line)[0][7:])
                kmer_coverage = float(re.findall('cov_[0-9]+\.[0-9]+', line)[0][4:])
                physical_lengths.append(kmer_length_header + KMER_LENGTH - 1)
    return(physical_lengths)

def get_kmer_coverages(filename):
    kmer_coverages = []
    with open(filename, 'r') as f:
        for index,line in enumerate(f):
            if line[0] == ">":
                kmer_coverages.append(float(re.findall('cov_[0-9]+\.[0-9]+', line)[0][4:]))
    return(kmer_coverages)

def get_num_of_contigs(filename):
    with open(filename,'r') as f:
        count = 0
        for line in f:
            if '>' in line:
                count+=1
    return(count)

def get_max_contig_length(filename):
    max_length = 0
    with open(filename,'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                current_sequence = []
            elif '>' not in line:
                current_sequence.append(line.strip())
            else:
                current_sequence = ''.join(current_sequence)
                current_length = len(current_sequence)
                if current_length > max_length:
                    max_length = current_length
                current_sequence = []
    return(max_length)

def get_total_length(filename):
    max_length = 0
    sum_length = 0
    lengths = []
    with open(filename,'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                total_sequence = []
            elif '>' not in line:
                total_sequence.append(line.strip())
    total_sequence = ''.join(total_sequence)
    len_sequence = len(total_sequence)
    return(len_sequence)

def get_mean_depth_coverage(filename, kmer_coverages, physical_lengths, KMER_LENGTH):
    cov_sum = 0
    for index, length in enumerate(physical_lengths):
        current_cov = (kmer_coverages[index] * length) / (length - KMER_LENGTH + 1)
        cov_sum += current_cov
    return(cov_sum/len(physical_lengths))

def get_n50(filename, total_length):
    max_length = 0
    all_contigs = []
    nuc_count = 0
    with open(filename,'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                current_sequence = []
            elif '>' not in line:
                current_sequence.append(line.strip())
            else:
                current_sequence = ''.join(current_sequence)
                all_contigs.append(current_sequence)
                current_sequence = []
        all_contigs.sort(key=len,reverse=True)
        perc_50 = total_length//2
        for contig in all_contigs:
            for char in contig:
                nuc_count += 1
                if nuc_count >= perc_50:
                    return(len(contig))

def get_contig_dist(filename):
    dist = {}
    with open(filename, 'r') as f:
        for index,line in enumerate(f):
            if index == 0:
                current_sequence = []
            elif '>' not in line:
                current_sequence.append(line.strip())
            else:
                current_sequence = ''.join(current_sequence)
                current_length = len(current_sequence)
                lowest_100 = int(current_length/100) * 100
                if lowest_100 in dist:
                    dist[lowest_100] += 1
                else:
                    dist[lowest_100] = 1
                current_sequence = []
    print('Contig Length' + '\t' + 'Number of contigs in this category')
    for key in sorted(dist.keys()):
        print(str(key) + '\t' + str(dist[key]))
    return(dist)


def main():
    filename = '/gpfs/projects/bgmp/shared/Bi621/contigs.fa'
    KMER_LENGTH = 49
    num_contigs = get_num_of_contigs(filename)
    print('\n')
    print('There are ' + str(num_contigs) + ' contigs in this file.' )
    max_contig_length = get_max_contig_length(filename)
    print('The maximum contig length is ' + str(max_contig_length) + ' nucleotides.')
    total_length = get_total_length(filename)
    print('The total length of this file is ' + str(total_length) + ' nucleotides.')
    mean_contig_length = total_length / num_contigs
    print('The mean contig length is approximately ' + str(round(mean_contig_length,2)) + ' nucleotides.')
    kmer_coverages = get_kmer_coverages(filename)
    physical_lengths = get_physical_lengths(filename, KMER_LENGTH)
    mean_depth_coverage = get_mean_depth_coverage(filename, kmer_coverages, physical_lengths, KMER_LENGTH)
    print('The mean depth of coverage is approximately ' + str(round(mean_depth_coverage,2)) + 'x.')
    n50 = get_n50(filename, total_length)
    print('The N50 for this dataset is ' + str(n50) + '.')
    print('\n***CONTIG DISTRIBUTION***')
    contig_dist = get_contig_dist(filename)

if __name__ == "__main__":
    main()
