import matplotlib.pyplot as plt
from matplotlib import use
import os


def make_graphs(out_file):
    data_dicts = []
    temp_dict = {}
    not_data = 0
    with open(out_file, 'r') as f:
        lines = f.readlines()
    beginning_indexes = [i+1 for i, x in enumerate(lines) if 'Contig Length' in x]
    end_indexes = [i for i, x in enumerate(lines) if (('[0.000000]' in x) | ('[0.000001]' in x)) & (('Reading FastQ file' in x) | ('Reading graph file' in x))]
    end_indexes.append(len(lines)-1)
    print(beginning_indexes)
    print(end_indexes)





if __name__ == "__main__":
    make_graphs('PS6.out')
