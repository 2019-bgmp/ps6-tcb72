import re
import argparse
def get_args():
    # set parser
    parser = argparse.ArgumentParser(description='A program to complete PS6.')
    # -k will be for kmer_length, must be int
    parser.add_argument('-k', '--kmer_length', type=int, help = 'How long you want the kmer size to be.')
    # -f will be for filename, must be string
    parser.add_argument('-f','--filename', type=str, help = 'Filename of your contigs.fa file.')
    # parse the user input
    args = parser.parse_args()
    # get filename and kmer length out
    filename = args.filename
    kmer_length = args.kmer_length
    #return args in array
    return(filename, kmer_length)

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
    filename, KMER_LENGTH = get_args()
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
