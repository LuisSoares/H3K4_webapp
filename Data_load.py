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
suts=load_gene_list('suts.txt',name=False)


data_sets_names = ['wt_me3', 'wt_me2','spp1_me3','spp1_me2','swd2_me3','swd2_me2'
    ,'set2_me3','set2_me2','tbp','h4ac','h3','ncb2']
location_sets=['genes','CUTS','SUTS']


hdf5_file=h5py.File("datasets.hdf5","w")

def create_hdf5_track(dataset_name,hdf5_file_handler):
    temp=hdf5_file_handler.create_group(dataset_name)
    for n,chr in enumerate(list_of_chromosomes):
        temp.create_dataset(chr,(size_of_chromosomes[n],),dtype='f')


def create_location_datatsets():
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

for item in data_sets_names:
    create_hdf5_track(item,hdf5_file)
create_location_datatsets()
load_track_from_bdg('bar40_norm.bdg','wt_me3')
load_track_from_bdg('bar49_norm.bdg','wt_me2')
load_track_from_bdg('bar41_norm.bdg','spp1_me3')
load_track_from_bdg('bar50_norm.bdg','spp1_me2')
load_track_from_bdg('bar42_norm.bdg','swd2_me3')
load_track_from_bdg('bar51_norm.bdg','swd2_me2')
load_track_from_bdg('BAR28_treat_pileup.bdg','set2_me3')
load_track_from_bdg('BAR29_treat_pileup.bdg','set2_me2')
load_track_from_bdg('TBP_treat_pileup.bdg','tbp')
load_track_from_bdg('H4ac_treat_pileup.bdg','h4ac')
load_track_from_bdg('H3_treat_pileup.bdg','h3')
load_track_from_bdg('Ncb2_treat_pileup.bdg','ncb2')
create_genes_tracks(genes,'genes')
create_genes_tracks(cuts,'CUTS')
create_genes_tracks(suts,'SUTS')
