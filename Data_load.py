__author__ = 'Luis'
from intermine.webservice import Service
from spyre import server
import codecs
#import requests
# from spyre.server import Site, App
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from collections import OrderedDict
import os
import pickle


#Setting plotting parameters
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['lines.color'] = 'r'
mpl.rcParams['axes.facecolor'] = "#E8EBEA"
mpl.rcParams['axes.edgecolor'] = "white"
mpl.rcParams['axes.grid'] = True
mpl.rcParams['grid.color'] = 'white'
mpl.rcParams['grid.linestyle'] = '-'
mpl.rcParams['grid.alpha'] = 0.5
mpl.rcParams['grid.linewidth'] = 1.5
mpl.rcParams['axes.axisbelow'] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'HELVETICA'
mpl.rcParams['text.color'] = '#3B3B3B'
mpl.rcParams['axes.color_cycle'] = ['#E55934', '#5BC0EB']
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.fontsize'] = "medium"
mpl.rcParams['legend.shadow'] = False
mpl.rcParams['xtick.major.size'] = 0
mpl.rcParams['ytick.major.size'] = 0

# Setting encoding for HTML File
ENCODING = 'utf-8'

# This block allows to override the Default Spyre HTML Template
# The template in now called view.html but is present in the home directory
class CustomView(server.View.View):
    def getHTML(self):
        file_path = 'view.html'
        f = codecs.open(file_path, 'r', ENCODING)
        html = f.read()
        f.close()
        return html

server.View.View = CustomView

# Genomic stuff
# TODO: Maybe create a file containing the names of datasets to make it trully costumizable upon deployement
# The lists of chromosome names and of size is for S.cerevisiae

list_of_chromosomes = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV',
                          'XVI', 'MT']
size_of_chromosomes = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816,
                           1078177, 924431, 784333, 1091291, 948066, 85779]

# Loads chromosomes with track from BDG file
# the BDG file is a normal BED FILE, the SEQ_NAME collumn should contain the names of the chromosomes without "chr"
def load_track_from_bdg(track_file):
    chromosome_dict = OrderedDict()
    for n in range(len(list_of_chromosomes)):
        chromosome_dict[list_of_chromosomes[n]] = np.zeros(size_of_chromosomes[n])
    with open(track_file) as data:
        for line in data:
            # if 'MT' in line:continue
            line_list = line.rstrip().split('\t')
            chromosome_dict[line_list[0]][int(line_list[1]):int(line_list[2])] = float(line_list[3])
    return chromosome_dict


# Only runs the function the first time, otherwise relies on pickled list
if not os.path.exists('datasets.pkl'):
    wt_me3 = load_track_from_bdg('bar40_norm.bdg')
    wt_me2 = load_track_from_bdg('bar49_norm.bdg')
    data_sets = [wt_me3, wt_me2]
    pickle.dump(data_sets, open('datasets.pkl', 'wb'))


data_sets = pickle.load(open('datasets.pkl', 'rb'))

data_sets_names = ['H3K4me3', 'H3K4me2']
data_sets_legend = ['H3K4me$^3$', 'H3K4me$^2$']

# Loads gene list from BED file
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
# Loads chromosomes with Genes from BDG file, each chromosome is a 2D numpy array with first row
# being the + strand genes, and the second row the genes in the - strand
def create_genes_tracks(gene_list):
    chromosome_dict = OrderedDict()
    for n in range(len(list_of_chromosomes)):
        chromosome_dict[list_of_chromosomes[n]] = np.zeros((2,size_of_chromosomes[n]), dtype='b')
    for line in gene_list:
        if line[4] == '+':
            chromosome_dict[line[1]][0,(line[2]):int(line[3])] = 1
        else:
            chromosome_dict[line[1]][1,int(line[2]):int(line[3])] = -1
    return chromosome_dict

# TODO:maybe change to native numpy array save
# Only runs the function the first time, otherwise relies on pickled array
if not os.path.exists('genes.pkl'):
    genes_track = create_genes_tracks(genes)
    pickle.dump(genes_track, open('genes.pkl', 'wb'))
if not os.path.exists('cuts.pkl'):
    genes_track = create_genes_tracks(cuts)
    pickle.dump(genes_track, open('cuts.pkl', 'wb'))

genes_track = pickle.load(open('genes.pkl', 'rb'))
cuts_track = pickle.load(open('cuts.pkl', 'rb'))

