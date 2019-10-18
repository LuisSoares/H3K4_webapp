__author__ = 'Luis'

from spyre import server
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
from matplotlib.colors import PowerNorm
from scipy.cluster.vq import kmeans,vq
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde
from bokeh import plotting
#from bokeh.models import OpenURL,Callback,Rect,ColumnDataSource,CustomJS,HoverTool,Range1d,PanTool,WheelZoomTool,TapTool,BoxSelectTool,BoxZoomTool
from bokeh.embed import components
from bokeh.resources import INLINE
from bokeh.resources import CDN
from matplotlib.cm import Reds as Reds

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % (rgb[0]*255,rgb[1]*255,rgb[2]*255)

hex=[]
for i in range(255):
    hex.append(rgb_to_hex(Reds(i)))

PATH=os.path.dirname(os.path.realpath(__file__))
#Setting plotting parameters
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['lines.color'] = 'r'
mpl.rcParams['axes.facecolor'] = "#E8EBEA"
mpl.rcParams['axes.edgecolor'] = "black"
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
with open('{}/Steinmetz_YPD_levels.txt'.format(PATH)) as data:
    steinmetz=dict()
    try:
        while True:
            temp=data.readline().split('\t')
            #print(temp)
            steinmetz[temp[0].upper()]=float(temp[1].rstrip())
    except IndexError:pass
#print(genes[:10])
data_sets=h5py.File("{}/datasets.hdf5".format(PATH))

data_sets_names = ['wt_me3', 'wt_me2','spp1_me3','spp1_me2',
                   'swd2_me3','swd2_me2','set2_me3','set2_me2',1,2,3,4,'h3k36_me3','h3k4_me1']
location_sets=['genes','CUTS','SUTS','ORFS']
data_sets_legend = ['H3K4me$^3$', 'H3K4me$^2$','H3K4me$^3$ ($\Delta$spp1)','H3K4me$^2$ ($\Delta$spp1)'
    ,'H3K4me$^3$ ($\Delta$swd2)','H3K4me$^2$ ($\Delta$swd2)','H3K4me$^3$ ($\Delta$set2)'
    ,'H3K4me$^2$ ($\Delta$set2)',1,2,3,4,'H3K36me$^3$','H3K4me$^1$']

class HeatMap(server.App):
    def __init__(self):
        self.data=None
    tabs = ["Static", "Interactive"]
    title = "Heat Map"
    inputs = [
        {"input_type":'dropdown',
				"label": 'track',
				"options" : [
					{"label": "Choose Dataset", "value":"empty"},
					{"label": "H3K4me3", "value":0},
					{"label": "H3K4me2", "value":1},
					{"label": "H3K4me3 spp1", "value":2},
                    {"label": "H3K4me2 spp1", "value":3},
                    {"label": "H3K4me3 swd2", "value":4},
                    {"label": "H3K4me2 swd2", "value":5},
                    {"label": "H3K4me3 set2", "value":6},
                    {"label": "H3K4me2 set2", "value":7},
                    {"label": "H3K36me3", "value":12},
                {"label": "H3K4me1", "value":13}],
				"key": 'ticker',
				"linked_key":'custom_ticker',
				"linked_type":'text'},
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
            {	"input_type":'dropdown',
				"label": 'Sorting',
				"options" : [
					{"label": "None", "value":None},
					{"label": "Max value", "value":"max"},
					{"label": "K means clustering", "value":"clustering"},
                    {"label": "Expression Levels","value":"expression"},
                ],
                 "key": 'cluster_ticker',
				"linked_key":'cluster_ticker',
				"linked_type":'text'},
        {"input_type": "text",
               "label":'Number of clusters',
               "variable_name": "number_of_clusters",
               "value": 1,
         }
              ]

    controls = [{"control_type": "button",
                 "control_id": "button1",
                 "label": "Plot",
                }

                ]

    outputs = [{"output_type": "plot",
                "control_id": "button1",
                "output_id": "anchor_plot",
                "on_page_load": False,
                "tab":"Static"
                },
               {"output_type": "html",
                "control_id": "button1",
                "output_id": "anchor_plot2",
                "on_page_load": False,
                "tab":"Interactive"
                }
        ]
    #self.data=[]
    def getPlot(self, params):
        mpl.rcParams['axes.edgecolor'] = "black"
        fig = plt.figure(figsize=(3, 6))
        gs = gridspec.GridSpec(2, 2, height_ratios=[6,2])
        gs.update(hspace=0.05)
        ax1 = plt.subplot(gs[0, :])
        ax2 = plt.subplot(gs[1, :])
        up = int(params['upstream_range'])
        down=int(params['downstream_range'])
        sorting_mode=params['cluster_ticker']
        number_of_clusters=params['number_of_clusters']
        #print('*'*100)
        #print(number_of_clusters)
        #print(genes[:10])
        dataset=params['ticker']
        data=np.zeros((7000,(abs(up)+abs(down))))
        #print(data)
        number_of_genes=0
        list_of_genes=[]
        for gene in genes2:
                if gene[4]=='+':
                    try:
                        to_add=data_sets[data_sets_names[int(dataset)]][gene[1]][gene[2]+up:gene[2]+down]
                        data[number_of_genes]=to_add
                        number_of_genes+=1
                        list_of_genes.append(gene[0])
                    except Exception as e:
                        #print(e)
                        pass
                else:
                    try:
                        to_add=data_sets[data_sets_names[int(dataset)]][gene[1]][gene[3]-down:gene[3]-up][::-1]
                        data[number_of_genes]=to_add
                        number_of_genes+=1
                        list_of_genes.append(gene[0])
                    except Exception as e:
                        pass
        #print(data)
        #print(number_of_genes)
        data=data[:number_of_genes]
        #print("Enter rolling")
        #print(data.shape)
        data=data[:,:data.shape[1]-data.shape[1]%10]
        #print(data.shape)
        def movingaverage (values):
            return np.mean(values.reshape(-1, 10), axis=1)
        data = np.apply_along_axis(movingaverage,1,data)
        #print("exit rolling")
        if sorting_mode=='max':
            data=data[np.argsort(np.array([np.mean(item) for item in data]))][::-1]
        elif sorting_mode=="clustering":
            #print("Enter cluster")
            n=int(number_of_clusters)
            centroids,_ = kmeans(data,n)
            idx,_ = vq(data,centroids)
            clusters=np.vstack((data[idx==index,:] for index in range(n)))
            data=clusters
            #print("Exit cluster")
        elif sorting_mode=='expression':
            expression=[]
            for item in list_of_genes:
                expression.append(steinmetz.get(item,0))
            data=data[np.argsort(np.array(expression))][::-1]
            #print(expression)
        ax1.imshow(data,aspect='auto',cmap='Reds')
        ax1.grid(False)
        ax1.axes.get_xaxis().set_visible(False)
        self.data=data
        kde_min=0
        kde_max=np.max(data)
        kde_range=np.linspace(kde_min,kde_max,300)
        def kde(item):
            kernel=gaussian_kde(item)
            return kernel(kde_range)
        #print(data.shape)
        #print("enter kde")
        data3=np.apply_along_axis(kde,1,data.T)
        #print("exit kde")

        #data3=data3.T

        ax2.imshow(data3.T[::-1],cmap="hot_r",aspect='auto',norm=PowerNorm(gamma=1./3.),interpolation="none")
        #norm=PowerNorm(gamma=1./2.)
        ax2.grid(False)
        start, end = ax2.get_xlim()
        starty, endy = ax2.get_ylim()
        ax2.xaxis.set_ticks(np.linspace(start, end, 4))
        ax2.xaxis.set_ticklabels([int(x) for x in np.linspace(up, down, 4)])
        ax2.set_xlabel('Distance from TSS (nt)')
        ax2.yaxis.set_ticks(np.linspace(starty, endy, 4))
        ax2.yaxis.set_ticklabels([int(x) for x in np.linspace(np.min(data),np.max(data), 4)])
        ax1.set_title('{}'.format(data_sets_legend[int(dataset)]))
        return fig





    #Overrides the CSS with custom_style.css located on home directory
    def getCustomCSS(self):
        """Override this function
		returns:
		string of css to insert on page load
		"""
        with open('{}/custom_style.css'.format(PATH)) as style:
            return style.read()+INLINE.css_raw[0]

    # Java script for google analytics
    def getCustomJS(self):
        return INLINE.js_raw[0]+'''
(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');
  ga('create', 'UA-55629344-3', 'auto');
  ga('send', 'pageview');
  '''

    def getHTML(self, params):
        #Reds9=["#140000","#290000","#3d0000","#520000","#660000","#7a0000","#8f0000","#a80000","#b80000","#cc0000","#e00000","#f50000",
        #    "#ff0a0a","#ff1f1f","#ff3333","#ff4747","#ff5c5c","#ff7070","#ff8585","#ff9999","#ffadad","#ffc2c2","#ffd6d6","#ffebeb","#ffffff"]
        up = int(params['upstream_range'])
        down=int(params['downstream_range'])
        sorting_mode=params['cluster_ticker']
        number_of_clusters=params['number_of_clusters']
        #print('*'*100)
        #print(number_of_clusters)
        #print(genes[:10])
        dataset=params['ticker']
        while self.data !=None:
            p = plotting.figure(width=250, height=500,x_range=(0,self.data.shape[1]*10), y_range=(0,self.data.shape[0]))

# must give a vector of image data for image parameter
            p.image(image=[self.data[::-1]], x=0, y=0, dw=self.data.shape[1]*10, dh=self.data.shape[0], palette=hex)
            print(self.data.shape)

            script, div = components(p,CDN)
            html = "%s\n%s"%(script,div)
            #print(script)
            return html



if __name__ == '__main__':
    app = HeatMap()
    app.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))

