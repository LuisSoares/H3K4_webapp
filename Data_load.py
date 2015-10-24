__author__ = 'Luis'

import numpy as np
from collections import OrderedDict
import h5py


# Genomic stuff
# TODO: Maybe create a file containing the names of datasets to make it trully costumizable upon deployement
# The lists of chromosome names and of size is for S.cerevisiae

list_of_chromosomes = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV',
                          'XVI', 'MT']
size_of_chromosomes = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816,
                           1078177, 924431, 784333, 1091291, 948066, 85779]

# Loads chromosomes with track from BDG file
# the BDG file is a normal BED FILE, the SEQ_NAME collumn should contain the names of the chromosomes without "chr"

def load_gene_list(gene_file,name=True):
    gene_list = []
    genes = open(gene_file)
    genes.readline()
    for line in genes:
        temp = line.split('\t')
        if name==True:
            gene_list.append([temp[4].strip().upper(), temp[0][3:], int(temp[1]), int(temp[2]), temp[3]])

        else:
            gene_list.append(['CUT', temp[0][3:], int(temp[1]), int(temp[2]), temp[3].strip()])
    return gene_list


genes = load_gene_list('genes.txt')
cuts=load_gene_list('CUTS.txt',name=False)


data_sets_names = ['wt_me3', 'wt_me2','spp1_me3','spp1_me2']
location_sets=['genes','CUTS']
data_sets_legend = ['H3K4me$^3$', 'H3K4me$^2$','H3K4me$^3$ ($\Delta$spp1)','H3K4me$^2$ ($\Delta$spp1)']

hdf5_file=h5py.File("datasets.hdf5","w")
for item in data_sets_names:
    temp=hdf5_file.create_group(item)
    for n,chr in enumerate(list_of_chromosomes):
        temp.create_dataset(chr,(size_of_chromosomes[n],),dtype='f')
for item in location_sets:
    temp=hdf5_file.create_group(item)
    for n,chr in enumerate(list_of_chromosomes):
        temp.create_dataset(chr,(2,size_of_chromosomes[n],),dtype='b')

def load_track_from_bdg(track_file,name):
    chromosome_dict = OrderedDict()
    for n in range(len(list_of_chromosomes)):
        chromosome_dict[list_of_chromosomes[n]] = np.zeros(size_of_chromosomes[n])
    with open(track_file) as data:
        for line in data:
            # if 'MT' in line:continue
            line_list = line.rstrip().split('\t')
            chromosome_dict[line_list[0]][int(line_list[1]):int(line_list[2])] = float(line_list[3])
    for key,value in chromosome_dict.items():
            print(key,value)
            hdf5_file[name][key][...] = value

def create_genes_tracks(gene_list,name):
    chromosome_dict = OrderedDict()
    for n in range(len(list_of_chromosomes)):
        chromosome_dict[list_of_chromosomes[n]] = np.zeros((2,size_of_chromosomes[n]))
    for line in gene_list:
        if line[4] == '+':
            chromosome_dict[line[1]][0,(line[2]):int(line[3])] = 1
        else:
            chromosome_dict[line[1]][1,int(line[2]):int(line[3])] = -1
    for key,value in chromosome_dict.items():
            print(key,value)
            hdf5_file[name][key][...] = value

load_track_from_bdg('bar40_norm.bdg','wt_me3')
load_track_from_bdg('bar49_norm.bdg','wt_me2')
load_track_from_bdg('bar41_norm.bdg','spp1_me3')
load_track_from_bdg('bar50_norm.bdg','spp1_me2')
create_genes_tracks(genes,'genes')
create_genes_tracks(cuts,'CUTS')