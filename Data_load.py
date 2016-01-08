__author__ = 'Luis'

import numpy as np
from collections import OrderedDict
import h5py
import os
import sys

PATH=os.path.dirname(os.path.realpath(__file__))
DEFAULT_TRACK_PATH=os.path.join(PATH,"Features")

if not sys.argv[1]:
    print('Usage: Data_load <dataset initialization file')
    sys.exit()

if not os.path.isfile(os.path.join(PATH,sys.argv[1])):
    print(os.path.isfile(os.path.join(PATH,sys.argv[1])))
    print ('Please indicate valid dataset file')
    sys.exit()

if not os.path.exists(os.path.join(PATH,'Datasets')):
    os.mkdir(os.path.join(PATH,'Datasets'))

HDF5_PATH=os.path.join(PATH,'Datasets')

# Genomic stuff
# TODO: Maybe create a file containing the names of datasets to make it trully costumizable upon deployement
# The lists of chromosome names and of size is for S.cerevisiae

list_of_chromosomes = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV',
                          'XVI', 'MT']
size_of_chromosomes = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816,
                           1078177, 924431, 784333, 1091291, 948066, 85779]

# Loads chromosomes with track from BDG file
# the BDG file is a normal BED FILE, the SEQ_NAME collumn should contain the names of the chromosomes without "chr"

def initialization(file_name):
    datasets=[]
    with open(file_name) as params_file:
        header=params_file.readline()
        #assert header=='Dset\tRep\tFile'
        for line in params_file:

            temp=line.split(' ')
            print(temp)
            datasets.append((temp[0],temp[1],temp[2]))

        return datasets


def load_gene_list(gene_file,name=True):
    gene_list = []
    genes = open(gene_file)
    genes.readline()
    for line in genes:
        temp = line.split('\t')
        if name==True:
            gene_list.append([temp[4].strip().upper(), temp[0][3:], int(temp[1]), int(temp[2]), temp[3]])

        else:
            gene_list.append(['ND', temp[0][3:], int(temp[1]), int(temp[2]), temp[3].strip()])
    return gene_list




'''data_sets_names = ['wt_me3', 'wt_me2','spp1_me3','spp1_me2','swd2_me3','swd2_me2'
    ,'set2_me3','set2_me2','tbp','h4ac','h3','ncb2']'''
location_sets=['genes','CUTS','SUTS','ORFS']




def create_hdf5_track(dataset_name,dataset_representation,hdf5_file_handler):
    temp=hdf5_file_handler.create_group(dataset_name)
    temp.attrs['dataset_representation']=dataset_representation
    for n,chr in enumerate(list_of_chromosomes):
        temp.create_dataset(chr,(size_of_chromosomes[n],),dtype='f',compression='gzip', compression_opts=9)


def create_location_datatsets(dataset_file):
    for item in location_sets:
        temp=dataset_file.create_group(item)
        for n,chr in enumerate(list_of_chromosomes):
            temp.create_dataset(chr,(2,size_of_chromosomes[n],),dtype='b',compression='gzip', compression_opts=9)

def load_track_from_bdg(track_file,name,fast=True):
    """Fast should be used for dense data,
    slow should be used for large genomes with sparse data
    TODO: maybe slow is too slow to be usable"""
    with open(track_file) as data:
        chromosome_dict = OrderedDict()
        if fast:
            for n in range(len(list_of_chromosomes)):
                chromosome_dict[list_of_chromosomes[n]] = np.zeros(size_of_chromosomes[n])
            for line in data:
                line_list = line.rstrip().split('\t')
                chromosome_dict[line_list[0]][int(line_list[1]):int(line_list[2])]=float(line_list[3])
            for key,value in chromosome_dict.items():
                print(key,value)
                hdf5_file[name][key][...] = value
        else:
            for line in data:
                # if 'MT' in line:continue
                line_list = line.rstrip().split('\t')
                hdf5_file[name][line_list[0]][int(line_list[1]):int(line_list[2])] = float(line_list[3])


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

hdf5_file=h5py.File(os.path.join(HDF5_PATH,"datasets_test.hdf5"),"w")

genes = load_gene_list(os.path.join(DEFAULT_TRACK_PATH,'genes.txt'))
cuts=load_gene_list(os.path.join(DEFAULT_TRACK_PATH,'CUTS.txt'),name=False)
suts=load_gene_list(os.path.join(DEFAULT_TRACK_PATH,'suts.txt'),name=False)
orfs=load_gene_list(os.path.join(DEFAULT_TRACK_PATH,'ORF.txt'),name=False)

if __name__ == '__main__':
    create_location_datatsets(hdf5_file)
    create_genes_tracks(genes,'genes')
    create_genes_tracks(cuts,'CUTS')
    create_genes_tracks(suts,'SUTS')
    create_genes_tracks(orfs,'ORFS')
    print(initialization(sys.argv[1]))
    for item in initialization(sys.argv[1]):
        print(item)
        create_hdf5_track(item[0],item[1],hdf5_file_handler=hdf5_file)
        load_track_from_bdg(item[2],item[0])


