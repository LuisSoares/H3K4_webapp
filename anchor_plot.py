__author__ = 'Luis'
from intermine.webservice import Service
from spyre import server
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import h5py

PATH=os.path.dirname(os.path.realpath(__file__))
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


genes2 = load_gene_list('{}/genes2.txt'.format(PATH))


data_sets=h5py.File("{}/datasets.hdf5".format(PATH),'r')

data_sets_names = ['wt_me3', 'wt_me2','spp1_me3','spp1_me2',
                   'swd2_me3','swd2_me2','set2_me3','set2_me2',1,2,3,4,'h3k36_me3','h3k4_me1']
location_sets=['genes','CUTS','SUTS','ORFS']
data_sets_legend = ['H3K4me$^3$', 'H3K4me$^2$','H3K4me$^3$ ($\Delta$spp1)','H3K4me$^2$ ($\Delta$spp1)'
    ,'H3K4me$^3$ ($\Delta$swd2)','H3K4me$^2$ ($\Delta$swd2)','H3K4me$^3$ ($\Delta$set2)'
    ,'H3K4me$^2$ ($\Delta$set2)',1,2,3,4,'H3K36me$^3$','H3K4me$^1$']

class AnchorPlot(server.App):
    title = "Anchor Plot"
    inputs = [
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
                   {"label": "H3K36me3","value":12},
                   {"label": "H3K4me1","value":13}
               ],
               "variable_name": 'check_boxes'},
              {"input_type": 'slider',
               "label": 'Upstream Range (Relative to TSS)',
               "min": -2000,
               "max": 0,
               "value": 1000,
               "variable_name": 'upstream_range',
               "action_id": 'plot'},
              {"input_type": 'slider',
               "label": 'Downstream Range (Relative to TSS)',
               "min": 0,
               "max": 4000,
               "value": 1000,
               "variable_name": 'downstream_range',
               "action_id": 'plot'},
              {"input_type": "text",
               "label":'Restriction (search word present in YGD description of gene, e.g.:ribosomal)',
               "variable_name": "input_gene",
                             },]

    controls = [{"control_type": "button",
                 "control_id": "button1",
                 "label": "Plot",
                }

                ]

    outputs = [{"output_type": "plot",
                "control_id": "button1",
                "output_id": "anchor_plot",
                "on_page_load": False,
                }
        ]

    def getPlot(self, params):
        fill_colors=['#921B16', '#D6701C','#251F47','#68A691','#68A691',
                     '#447604','#232621','#0BBE30','#BC3908','#191923','#391923',
                     '#591923','#921B1C']
        fig,ax = plt.subplots(1,1,figsize=(6, 6))
        up = int(params['upstream_range'])
        down=int(params['downstream_range'])
        gene_restr=params['input_gene']
        print('*'*100)
        print('gene_list',gene_restr)
        gene_list=[]
        if gene_restr != u'':
            service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")
            query = service.new_query("Gene")
            query.add_view("secondaryIdentifier")
            query.add_constraint("description", "CONTAINS", str(gene_restr), code = "A")
            for row in query.rows():
                gene_list.append(row["secondaryIdentifier"])
        print(len(gene_list),gene_list)
        datasets_to_plot=params['check_boxes']
        for n,item in enumerate(datasets_to_plot):
            current_data=self.create_anchor_data(item,up,down,[])
            ax.plot(current_data[0],label=data_sets_legend[int(item)],color=fill_colors[n],alpha=0.7,lw=2)
            if gene_list!=[]:
                current_data2=self.create_anchor_data(item,up,down,gene_list)
                ax.plot(current_data2[0],'--',label=data_sets_legend[int(item)]+'\n constrained to {} genes'.format(current_data2[1])
                                        ,color=fill_colors[n],alpha=0.7,lw=2)

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        ax.set_xlim(0,abs(up)+abs(down))
        ax.set_title("Average of {} genes plotted".format(current_data[1]))
        ax.set_ylabel('Average normalized counts per million reads')
        ax.set_xlabel('Distance in nucleotides ({} is TSS)'.format(-up))
        return ax

    def create_anchor_data(self,dataset,up,down,restrictions=None):
        data=np.zeros(abs(up)+abs(down))
        number_of_genes=0
        number_of_genes_plotted=0
        if restrictions==[]:
            for gene in genes2:
                number_of_genes+=1
                if gene[4]=='+':
                    try:
                        to_add=data_sets[data_sets_names[int(dataset)]][gene[1]][gene[2]+up:gene[2]+down]
                        number_of_genes_plotted+=1
                        data=data+to_add
                    except Exception as e:
                        pass
                else:
                    try:
                        to_add=data_sets[data_sets_names[int(dataset)]][gene[1]][gene[3]-down:gene[3]-up][::-1]
                        number_of_genes_plotted+=1
                        data=data+to_add
                    except Exception as e:
                        pass
        else:
            for gene in genes2:
                if gene[0] in restrictions:
                    number_of_genes+=1
                    if gene[4]=='+':
                        try:
                            to_add=data_sets[data_sets_names[int(dataset)]][gene[1]][gene[2]+up:gene[2]+down]
                            number_of_genes_plotted+=1
                            data=data+to_add
                        except Exception as e:
                            pass
                    else:
                        try:
                            to_add=data_sets[data_sets_names[int(dataset)]][gene[1]][gene[3]-down:gene[3]-up][::-1]
                            number_of_genes_plotted+=1
                            data=data+to_add
                        except Exception as e:
                            pass
        print(number_of_genes_plotted)
        return (data/number_of_genes_plotted,number_of_genes_plotted)





    #Overrides the CSS with custom_style.css located on home directory
    def getCustomCSS(self):
        """Override this function
		returns:
		string of css to insert on page load
		"""
        with open('{}/custom_style.css'.format(PATH)) as style:
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




if __name__ == '__main__':
    app = AnchorPlot()
    app.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))

