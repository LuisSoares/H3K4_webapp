__author__ = 'Luis'
from intermine.webservice import Service
from spyre import server
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import h5py
import time


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
mpl.rcParams['axes.color_cycle'] = ['#921B16', '#D6701C','#251F47']
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.fontsize'] = "medium"
mpl.rcParams['legend.shadow'] = False
mpl.rcParams['xtick.major.size'] = 0
mpl.rcParams['ytick.major.size'] = 0

# Setting encoding for HTML File
ENCODING = 'utf-8'




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


genes = load_gene_list('/home/lint78/web_site/genes.txt')



data_sets=h5py.File("/home/lint78/web_site/datasets.hdf5","r")
data_sets_names = ['wt_me3', 'wt_me2','spp1_me3','spp1_me2',
                   'swd2_me3','swd2_me2','set2_me3','set2_me2','tbp','h4ac','h3','ncb2']
location_sets=['genes','CUTS','SUTS','ORFS']
data_sets_legend = ['H3K4me$^3$', 'H3K4me$^2$','H3K4me$^3$ ($\Delta$spp1)','H3K4me$^2$ ($\Delta$spp1)'
    ,'H3K4me$^3$ ($\Delta$swd2)','H3K4me$^2$ ($\Delta$swd2)','H3K4me$^3$ ($\Delta$set2)'
    ,'H3K4me$^2$ ($\Delta$set2)','TBP','H4Ac','H3','NCB2']

class SimpleApp(server.App):
    title = "Gene plotter"
    tabs = ["Plot", "Summary", "SGD"]
    inputs = [{"input_type": "text",
               "label":'Gene Code/Name',
               "variable_name": "freq",
               "value": 'YLR249W',

              },
              {"input_type": 'checkboxgroup',
               "label":"Track",
               "options": [
                   {"label": "H3K4me3", "value": 0, "checked": True},
                   {"label": "H3K4me2", "value": 1,},
                   {"label": "H3K4me3 (&#916;spp1)","value":2},
                   {"label": "H3K4me2 (&#916;spp1)","value":3},
                   {"label": "H3K4me3 (&#916;swd2)","value":4},
                   {"label": "H3K4me2 (&#916;swd2)","value":5},
                   {"label": "H3K4me3 (&#916;set2)","value":6},
                   {"label": "H3K4me2 (&#916;set2)","value":7},
                   {"label": "TBP(Yoo Jin)","value":8},
                   {"label": "H4Ac(Yoo Jin)","value":9},
                   {"label": "H3(Yoo Jin)","value":10},
                   {"label": "NCB2(Yoo Jin)","value":11}
               ],
               "variable_name": 'check_boxes'},
              {"input_type": 'slider',
               "label": 'Extend Range',
               "min": 0,
               "max": 4000,
               "value": 0,
               "variable_name": 'range',
               "action_id": 'plot'}]

    controls = [{"control_type": "button",
                 "control_id": "button1",
                 "label": "Plot",
                }

                ]

    outputs = [{"output_type": "plot",
                "control_id": "button1",
                "output_id": "sine_wave_plot",
                "on_page_load": False,
                'tab': 'Plot'}
        ,       {"output_type": "html",
                "control_id": "button1",
                "output_id": "custom_html",
                "tab": "SGD"},
               {"output_type": "table",
                "output_id": "table_id",
                "control_id": "button1",
                "tab": "Summary"}

               ]

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
        start_time=time.time()
        global current_gene
        fig = plt.figure(figsize=(9, 7))
        gs = gridspec.GridSpec(4, 2, height_ratios=[6, 1,1,1])
        gs.update(hspace=0.05)
        splt1 = plt.subplot(gs[0, :])
        gene = params['freq'].upper()
        range_ext = int(params['range'])
        searched='False'
        if gene.startswith('CHR'):
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
        open('/home/lint78/web_site/log.txt','a').write(str(current_gene)+'\t'+str(params['check_boxes'])+'\t'+searched+'\n')
        splt1.set_title('{}\nChromosome {}:{}:{}, strand:{}'.format(current_gene[0], current_gene[1], current_gene[2],
                                                                    current_gene[3], current_gene[4]))
        # splt1.set_xlabel('Distance from TSS (nt)')
        splt1.axes.xaxis.set_ticklabels([])
        splt1.set_ylabel('Normalized counts per million reads')
        n=0
        fill_colors=['#921B16', '#D6701C','#251F47','#68A691','#68A691',
                     '#447604','#232621','#0BBE30','#BC3908','#191923','#391923',
                     '#591923']
        for item in params['check_boxes']:
            if current_gene[4] == '+':
                splt1.plot(
                    data_sets[data_sets_names[int(item)]][current_gene[1]][current_gene[2] - range_ext:current_gene[3] + range_ext],
                    label=data_sets_legend[int(item)],color=fill_colors[n])
                splt1.fill_between(range(0,len(data_sets[data_sets_names[int(item)]][current_gene[1]][current_gene[2] - range_ext:current_gene[3] + range_ext])),
                0, data_sets[data_sets_names[int(item)]][current_gene[1]][current_gene[2] - range_ext:current_gene[3] + range_ext],
                                   alpha=0.2,color=fill_colors[n])
            else:
                splt1.plot(
                    data_sets[data_sets_names[int(item)]][current_gene[1]][current_gene[2] -range_ext:current_gene[3] +range_ext][::-1],
                    label=data_sets_legend[int(item)],color=fill_colors[n])
                splt1.fill_between(range(0,len(data_sets[data_sets_names[int(item)]][current_gene[1]][current_gene[2] - range_ext:current_gene[3] + range_ext][::-1])),
                0, data_sets[data_sets_names[int(item)]][current_gene[1]][current_gene[2] - range_ext:current_gene[3] + range_ext][::-1],
                                   alpha=0.2,color=fill_colors[n])
            n+=1
        leg=splt1.legend(loc='upper left', bbox_to_anchor=(1, 1))
        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)
        size = (current_gene[3] + range_ext) - (current_gene[2] - range_ext)
        splt1.set_xlim(0, size)
        splt1.text(0.02, 0.02, 'Luis Soares', fontsize=6, transform=splt1.transAxes, verticalalignment='top')

        splt2 = plt.subplot(gs[1, :])
        if current_gene[4] == '+':
            y_values1 = data_sets['genes'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext]
            y_values2 = data_sets['genes'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext]
            y_values3 = data_sets['ORFS'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext]
            y_values4 = data_sets['ORFS'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext]
        else:
            y_values1 =  data_sets['genes'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
            y_values2 =  data_sets['genes'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
            y_values3 =  data_sets['ORFS'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
            y_values4 =  data_sets['ORFS'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
        x_values = range(0, len(y_values1))
        #splt2.plot(y_values, 'grey', linewidth=2)

        splt2.plot([0, len(y_values1)], [0, 0], 'black', linewidth=1.5)
        if current_gene[4] == '+':
            text_pos=0.71
        else: text_pos=0.25
        splt2.text(0.5, text_pos,current_gene[0],horizontalalignment='center',verticalalignment='center',transform = splt2.transAxes)
        splt2.fill_between(x_values, y_values1/4.0, y_values1, color='#FFB30F',alpha=0.5)
        splt2.fill_between(x_values, y_values2/4.0, y_values2, color='#FD151B',alpha=0.5)
        splt2.fill_between(x_values, y_values3/4.0, y_values3, color='#FFB30F',alpha=0.8)
        splt2.fill_between(x_values, y_values4/4.0, y_values4, color='#FD151B',alpha=0.8)
        splt2.set_ylim(-1.5, 1.5)
        splt2.text(1.01, 0.5, 'GENES', fontsize=12, transform=splt2.transAxes, verticalalignment='center')
        splt2.yaxis.set_tick_params(color='w')
        splt2.set_yticks([-0.5, 0.5])
        splt2.spines['left'].set_position(('outward', 10))
        splt2.set_yticklabels([])
        splt2.set_xlim(0, size)
        splt2.axes.xaxis.set_ticklabels([])
        #plot 3
        splt3 = plt.subplot(gs[2, :])
        if current_gene[4] == '+':
            y_values1 =  data_sets['CUTS'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext]
            y_values2 =  data_sets['CUTS'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext]
            splt3.text(0.01, 0.9, '>'*46, fontsize=13, transform=splt3.transAxes, verticalalignment='top',alpha=0.2)
            splt3.text(0.01, 0.46, '<'*46, fontsize=13, transform=splt3.transAxes, verticalalignment='top',alpha=0.2)
            splt2.text(0.01, 0.9, '>'*46, fontsize=13, transform=splt2.transAxes, verticalalignment='top',alpha=0.2)
            splt2.text(0.01, 0.46, '<'*46, fontsize=13, transform=splt2.transAxes, verticalalignment='top',alpha=0.2)
        else:
            y_values1 = data_sets['CUTS'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
            y_values2 = data_sets['CUTS'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
            splt3.text(0.01, 0.9, '<'*46, fontsize=13, transform=splt3.transAxes, verticalalignment='top',alpha=0.2)
            splt3.text(0.01, 0.46, '>'*46, fontsize=13, transform=splt3.transAxes, verticalalignment='top',alpha=0.2)
            splt2.text(0.01, 0.9, '<'*46, fontsize=13, transform=splt2.transAxes, verticalalignment='top',alpha=0.2)
            splt2.text(0.01, 0.46, '>'*46, fontsize=13, transform=splt2.transAxes, verticalalignment='top',alpha=0.2)
        x_values = range(0, len(y_values1))
        splt3.plot([0, len(y_values1)], [0, 0], 'black', linewidth=1.5)
        splt3.fill_between(x_values, y_values1/4.0, y_values1, color='#028090',alpha=0.75)
        splt3.fill_between(x_values, y_values2/4.0, y_values2, color='#6EEB83',alpha=0.75)
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
        splt3.axes.xaxis.set_ticklabels([])
        #plot 4
        splt4 = plt.subplot(gs[3, :])
        if current_gene[4] == '+':
            y_values1 =  data_sets['SUTS'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext]
            y_values2 =  data_sets['SUTS'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext]
            splt4.text(0.01, 0.9, '>'*46, fontsize=13, transform=splt4.transAxes, verticalalignment='top',alpha=0.2)
            splt4.text(0.01, 0.46, '<'*46, fontsize=13, transform=splt4.transAxes, verticalalignment='top',alpha=0.2)

        else:
            y_values1 = data_sets['SUTS'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
            y_values2 = data_sets['SUTS'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
            splt4.text(0.01, 0.9, '<'*46, fontsize=13, transform=splt4.transAxes, verticalalignment='top',alpha=0.2)
            splt4.text(0.01, 0.46, '>'*46, fontsize=13, transform=splt4.transAxes, verticalalignment='top',alpha=0.2)
        x_values = range(0, len(y_values1))
        splt4.plot([0, len(y_values1)], [0, 0], 'black', linewidth=1.5)
        splt4.fill_between(x_values, y_values1/4.0, y_values1, color='#028090',alpha=0.75)
        splt4.fill_between(x_values, y_values2/4.0, y_values2, color='#6EEB83',alpha=0.75)
        splt4.set_ylim(-1.5, 1.5)
        splt4.text(1.01, 0.5, 'SUTS', fontsize=12, transform=splt4.transAxes, verticalalignment='center')
        splt4.spines['bottom'].set_position(('outward', 10))
        splt4.spines['bottom'].set_color('black')
        splt4.spines['bottom'].set_linewidth(2)
        splt4.spines["bottom"].axis.axes.tick_params(direction="outward", length=5, color='black', width=2)
        splt4.xaxis.tick_bottom()
        splt4.yaxis.set_tick_params(color='w')
        splt4.set_yticks([])
        splt4.set_xlim(0, size)
        for tick in splt4.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        end_time=time.time()
        print('Elapsed time:'+str(end_time-start_time))
        return fig

    def getHTML(self, params):

        #return '<p>Send me an email:<a href="mailto:luis.miguel.mendes.soares@gmail.com?Subject=Gene%20Plotter" target="_top">Luis Soares</a></p>'
        return '<iframe src="http://www.yeastgenome.org/search?query={}" width=100% height="600"></iframe>'.format(current_gene[0])


    def getTable(self, params):
        df = pd.DataFrame(columns=['Track', 'Code', 'Chromosome', 'Start', 'End', 'Strand', 'Max', 'Max Pos'])

        gene = params['freq'].upper()
        if gene.startswith('CHR'):
            temp=gene.split(':')
            current_gene=['Region',temp[0][3:],int(temp[1]),int(temp[2]),'+']
        else:
            for item in genes:
                if item[0] == gene:
                    current_gene = item
                    break

            else:
                temp_gene=self.search_SGD(gene)
                searched='True'
                for item in genes:
                    if item[0] == temp_gene[0]:
                        current_gene = item
                        break
                else:
                    current_gene=temp_gene
        for n, item in enumerate(params['check_boxes']):
            wt_data = data_sets[data_sets_names[int(item)]][current_gene[1]][current_gene[2]:current_gene[3]]
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
        with open('/home/lint78/web_site/custom_style.css') as style:
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
    app = SimpleApp()
    app.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))
    # site.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))
    # app.launch(host='0.0.0.0')
