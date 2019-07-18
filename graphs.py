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
    title_indexes = [i for i, x in enumerate(lines) if 'Printing' in x]
    end_indexes.append(len(lines)-1)
    print(title_indexes)
    use('TkAgg')
    for i in range(len(beginning_indexes)):
        current_lines = [i.strip().split('\t') for i in lines[beginning_indexes[i]:end_indexes[i]]]
        x_values = [int(i[0]) for i in current_lines]
        y_values = [int(i[1]) for i in current_lines]
        plt.figure(figsize=(8,6))
        plt.bar(x_values,y_values,width=50)
        plt.xlabel('Contig length (bp)')
        plt.ylabel('Count')
        plt.title(lines[title_indexes[i]][21:])
        plt.savefig('plot'+str(i)+'.png')
        plt.show()





if __name__ == "__main__":
    make_graphs('PS6.out')
