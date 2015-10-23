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

list_of_chromosomes = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV',
                          'XVI', 'MT']
size_of_chromosomes = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816,
                           1078177, 924431, 784333, 1091291, 948066, 85779]

# Loads chromosomes with track from BDG file
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


class SimpleSineApp(server.App):
    title = "Gene plotter"
    tabs = ["Plot", "Summary", "Contact"]
    inputs = [{"input_type": "text",
               "label":'Gene Code/Name',
               "variable_name": "freq",
               "value": 'YLR249W',
               # "action_id":"sine_wave_plot"
              },
              {"input_type": 'checkboxgroup',
               "label":"Track",
               "options": [
                   {"label": "H3K4me3", "value": 0, "checked": True},
                   {"label": "H3K4me2", "value": 1}
               ],
               "variable_name": 'check_boxes'},
              {"input_type": 'slider',
               "label": 'Extend Range',
               "min": 0,
               "max": 2000,
               "value": 0,
               "variable_name": 'range',
               "action_id": 'plot'}]

    controls = [{"control_type": "button",
                 "control_id": "button1",
                 "label": "Plot",
                }]

    outputs = [{"output_type": "plot",
                "control_id": "button1",
                "output_id": "sine_wave_plot",
                "on_page_load": False, 'tab': 'Plot'}
        , {"output_type": "html",
           "control_id": "button1",
           "output_id": "custom_html",
           "tab": "Contact"},
               {"output_type": "table",
                "output_id": "table_id",
                "control_id": "button1",
                "tab": "Summary"}]

    def search_SGD(self,gene_code=None):
        service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")
        query = service.new_query("Gene")
        query.add_view(
        "chromosome.primaryIdentifier", "chromosomeLocation.end",
        "chromosomeLocation.start", "chromosomeLocation.strand",
        "secondaryIdentifier"
        )
        query.add_constraint("symbol", "=",gene_code, code = "A")
        for row in query.rows():
            print([row["secondaryIdentifier"],row["chromosome.primaryIdentifier"], row["chromosomeLocation.start"],row["chromosomeLocation.end"],
                    '+' if row["chromosomeLocation.strand"] else '-' ])
            return [row["secondaryIdentifier"],row["chromosome.primaryIdentifier"][3:], row["chromosomeLocation.start"],row["chromosomeLocation.end"],
                    '+' if row["chromosomeLocation.strand"] else '-' ]


    def getPlot(self, params):

        print(params)
        return self.plot(params)


    def plot(self, params):
        global current_gene
        fig = plt.figure(figsize=(9, 6))
        gs = gridspec.GridSpec(3, 2, height_ratios=[5, 1,1])
        gs.update(hspace=0.05)
        splt1 = plt.subplot(gs[0, :])
        gene = params['freq'].upper()
        range_ext = int(params['range'])
        searched='False'
        if gene.startswith('chr'):
            temp=gene.split(':')
            current_gene=['Region',temp[0][3:],int(temp[1]),int(temp[2]),'+']
        else:
            for item in genes:
                if item[0] == gene:
                    current_gene = item
                    break

            else:
                temp_gene=self.search_SGD(gene)
                print(temp_gene)
                searched='True'
                for item in genes:
                    if item[0] == temp_gene[0]:
                        current_gene = item
                        break
                else:
                    current_gene=temp_gene
                    print (gene)
#Todo add possibility to search SGD defaulting to SGD coordinates
        print('*'*40+'\nQuery by user:'+str(current_gene)+';\nrange extension:'
            +str(range_ext)+' ;Checked Boxes:'+str(params['check_boxes'])+';SGD used:'+searched+'\n'+'*'*40)
        #print(params)
        splt1.set_title('{}\nChromosome {}:{}:{}, strand:{}'.format(current_gene[0], current_gene[1], current_gene[2],
                                                                    current_gene[3], current_gene[4]))
        # splt1.set_xlabel('Distance from TSS (nt)')
        splt1.axes.xaxis.set_ticklabels([])
        splt1.set_ylabel('Normalized counts per million reads')
        n=0
        fill_colors=['#E55934', '#5BC0EB']
        for item in params['check_boxes']:
            if current_gene[4] == '+':
                splt1.plot(
                    data_sets[int(item)][current_gene[1]][current_gene[2] - range_ext:current_gene[3] + range_ext],
                    label=data_sets_legend[int(item)])
                splt1.fill_between(range(0,len(data_sets[int(item)][current_gene[1]][current_gene[2] - range_ext:current_gene[3] + range_ext])),
                0, data_sets[int(item)][current_gene[1]][current_gene[2] - range_ext:current_gene[3] + range_ext], alpha=0.2,color=fill_colors[n])
            else:
                splt1.plot(
                    data_sets[int(item)][current_gene[1]][current_gene[3] + range_ext:current_gene[2] - range_ext:-1],
                    label=data_sets_legend[int(item)])
                splt1.fill_between(range(0,len(data_sets[int(item)][current_gene[1]][current_gene[3] + range_ext:current_gene[2] - range_ext:-1])),
                0, data_sets[int(item)][current_gene[1]][current_gene[3] + range_ext:current_gene[2] - range_ext:-1], alpha=0.2,color=fill_colors[n])
            n+=1
        splt1.legend(loc='upper left', bbox_to_anchor=(1, 1))
        size = (current_gene[3] + range_ext) - (current_gene[2] - range_ext)
        splt1.set_xlim(0, size)
        #print(str(plt.ylim()[1] - 5))
        splt1.text(0.02, 0.02, 'Luis Soares', fontsize=6, transform=splt1.transAxes, verticalalignment='top')

        splt2 = plt.subplot(gs[1, :])
        if current_gene[4] == '+':
            y_values1 = genes_track[current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext]
            y_values2 = genes_track[current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext]
        else:
            y_values1 = genes_track[current_gene[1]][0,current_gene[3] + range_ext:current_gene[2] - range_ext:-1]
            y_values2 = genes_track[current_gene[1]][1,current_gene[3] + range_ext:current_gene[2] - range_ext:-1]
        x_values = range(0, len(y_values1))
        #splt2.plot(y_values, 'grey', linewidth=2)

        splt2.plot([0, len(y_values1)], [0, 0], 'black', linewidth=1.5)
        splt2.fill_between(x_values, y_values1/4, y_values1, color='#FFB30F')
        splt2.fill_between(x_values, y_values2/4, y_values2, color='#FD151B')
        splt2.set_ylim(-1.5, 1.5)
        splt2.text(1.01, 0.5, 'GENES', fontsize=12, transform=splt2.transAxes, verticalalignment='center')
        splt2.yaxis.set_tick_params(color='w')
        splt2.set_yticks([-0.5, 0.5])
        splt2.spines['left'].set_position(('outward', 10))
        splt2.set_yticklabels([])
        splt2.set_xlim(0, size)
        splt2.axes.xaxis.set_ticklabels([])

        splt3 = plt.subplot(gs[2, :])
        if current_gene[4] == '+':
            y_values1 = cuts_track[current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext]
            y_values2 = cuts_track[current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext]
            splt3.text(0.01, 0.9, '>'*46, fontsize=13, transform=splt3.transAxes, verticalalignment='top',alpha=0.2)
            splt3.text(0.01, 0.46, '<'*46, fontsize=13, transform=splt3.transAxes, verticalalignment='top',alpha=0.2)
            splt2.text(0.01, 0.9, '>'*46, fontsize=13, transform=splt2.transAxes, verticalalignment='top',alpha=0.2)
            splt2.text(0.01, 0.46, '<'*46, fontsize=13, transform=splt2.transAxes, verticalalignment='top',alpha=0.2)
        else:
            y_values1 = cuts_track[current_gene[1]][0,current_gene[3] + range_ext:current_gene[2] - range_ext:-1]
            y_values2 = cuts_track[current_gene[1]][1,current_gene[3] + range_ext:current_gene[2] - range_ext:-1]
            splt3.text(0.01, 0.9, '<'*46, fontsize=13, transform=splt3.transAxes, verticalalignment='top',alpha=0.2)
            splt3.text(0.01, 0.46, '>'*46, fontsize=13, transform=splt3.transAxes, verticalalignment='top',alpha=0.2)
            splt2.text(0.01, 0.9, '<'*46, fontsize=13, transform=splt2.transAxes, verticalalignment='top',alpha=0.2)
            splt2.text(0.01, 0.46, '>'*46, fontsize=13, transform=splt2.transAxes, verticalalignment='top',alpha=0.2)
        x_values = range(0, len(y_values1))
        #splt2.plot(y_values, 'grey', linewidth=2)
        splt3.plot([0, len(y_values1)], [0, 0], 'black', linewidth=1.5)
        splt3.fill_between(x_values, y_values1/4, y_values1, color='#028090',alpha=0.75)
        splt3.fill_between(x_values, y_values2/4, y_values2, color='#6EEB83',alpha=0.75)
        splt3.set_ylim(-1.5, 1.5)
        splt3.text(1.01, 0.5, 'CUTS', fontsize=12, transform=splt3.transAxes, verticalalignment='center')
        splt3.spines['bottom'].set_position(('outward', 10))
        splt3.spines['bottom'].set_color('black')
        splt3.spines['bottom'].set_linewidth(2)
        splt3.spines["bottom"].axis.axes.tick_params(direction="outward", length=5, color='black', width=2)
        splt3.xaxis.tick_bottom()
        splt3.yaxis.set_tick_params(color='w')
        splt3.set_yticks([])
        splt3.set_xlim(0, size)
        for tick in splt3.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        return fig

    def getHTML(self, params):

        return '<p>Send me an email:<a href="mailto:luis.miguel.mendes.soares@gmail.com?Subject=Gene%20Plotter" target="_top">Luis Soares</a></p>'
        # return '<iframe src="https://www.google.com" width="600" height="600"></iframe>'


    def getTable(self, params):
        df = pd.DataFrame(columns=['Track', 'Code', 'Chromosome', 'Start', 'End', 'Strand', 'Max', 'Max Pos'])
        gene = params['freq'].upper()
        if gene.startswith('chr'):
            temp=gene.split(':')
            current_gene=['Region',temp[0][3:],int(temp[1]),int(temp[2]),'+']
        else:
            for item in genes:
                if item[0] == gene:
                    current_gene = item
                    break

            else:
                temp_gene=self.search_SGD(gene)
                print(temp_gene)
                searched='True'
                for item in genes:
                    if item[0] == temp_gene[0]:
                        current_gene = item
                        break
                else:
                    current_gene=temp_gene
                    print (gene)
        for n, item in enumerate(params['check_boxes']):
            wt_data = data_sets[int(item)][current_gene[1]][current_gene[2]:current_gene[3]]
            max_peak = np.max(wt_data)
            if current_gene[4] == '+':
                max_pos = np.argmax(wt_data)
            else:
                wt_data = wt_data[::-1]
                max_pos = np.argmax(wt_data)
            df.loc[n] = [data_sets_names[int(item)], current_gene[0], current_gene[1], current_gene[2], current_gene[3],
                         current_gene[4], max_peak, max_pos]
        return df

    #Overrides the CSS with custom_style.css located on home directory
    def getCustomCSS(self):
        """Override this function
		returns:
		string of css to insert on page load
		"""
        with open('custom_style.css') as style:
            return style.read()

    # Java script for google analytics
    def getCustomJS(self):
        return '''
(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');
  ga('create', 'UA-55629344-3', 'auto');
  ga('send', 'pageview');
  '''

# class Index(App):
# def getHTML(self, params):
# return "Title Page Here"

if __name__ == '__main__':
    # site = Site(SimpleSineApp)
    app = SimpleSineApp()
    app.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))
    # site.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))
    # app.launch(host='0.0.0.0')