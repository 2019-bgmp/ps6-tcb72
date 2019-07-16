
def calculate_coverage(filenames):
    genome_length = 40000*50
    nuc_sum = 0
    for file in filenames:
        with open(file,'r') as f:
            for index, line in enumerate(f):
                if index%4==1:
                    line = line.strip()
                    nuc_sum+=len(line)
    return(nuc_sum/genome_length)

def calculate_mean_read_length(filenames):
    nuc_sum = 0
    record_counter = 0
    for file in filenames:
        with open(file,'r') as f:
            for index, line in enumerate(f):
                if index%4==1:
                    record_counter += 1
                    line = line.strip()
                    nuc_sum+=len(line)

    return(nuc_sum/record_counter)

def calculate_kmer_coverage(coverage, Lmean, kmer_lengths=[31,41,49]):
    for kmer_length in kmer_lengths:
        kmer_coverage = (coverage * (Lmean - kmer_length + 1)) / Lmean
        print('Kmer coverage for k = ' + str(kmer_length) +' is: ' + str(kmer_coverage))

if __name__ == "__main__":
    filenames = ['/gpfs/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1','/gpfs/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2','/gpfs/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched']
    coverage = calculate_coverage(filenames)
    print('Coverage is: ' + str(coverage))
    Lmean = calculate_mean_read_length(filenames)
    calculate_kmer_coverage(coverage, Lmean)
