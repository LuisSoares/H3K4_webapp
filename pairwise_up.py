__author__ = 'Luis'
import os
import h5py
from spyre import server
import numpy as np
from bokeh import plotting
from bokeh.models import OpenURL,Callback,Rect,ColumnDataSource,CustomJS,HoverTool,Range1d,PanTool,WheelZoomTool,TapTool,BoxSelectTool,BoxZoomTool
from bokeh.embed import components
from bokeh.resources import INLINE
from bokeh.resources import CDN
from collections import Counter

PATH=os.path.dirname(os.path.realpath(__file__))

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



genes = load_gene_list('{}/genes.txt'.format(PATH))


data_sets=h5py.File("{}/datasets.hdf5".format(PATH))
data_sets_names = ['wt_me3', 'wt_me2','spp1_me3','spp1_me2',
                   'swd2_me3','swd2_me2','set2_me3','set2_me2','tbp','h4ac','h3','ncb2']
data_sets_display_names=['WT me3','WT me2','SPP1 me3','SPP1 me2',
                         'SWD2 me3','SWD2 me2','SET2 me3','SET2 me2','TBP','H4Ac','H3','NCB2']
location_sets=['genes']

def get_max(dataset1,dataset2=None,range1=None,range2=None,restrict=False):
    x,y,z,w,wx,wy=[],[],[],[],[],[]
    w=['#0029DE' for _ in range(len(genes))]
    for n,gene in enumerate(genes):
        if gene[4]=='+':
            limit1=gene[2]+range1
            limit2=gene[2]+range2
            if restrict:

                if gene[3]-gene[2]<range1:
                    limit1=gene[3]
                if gene[3]-gene[2]<range2:
                    limit2=gene[3]
            try:
                temp=np.max(data_sets[dataset1][gene[1]][gene[2]:limit1])
                temp2=np.argmax(data_sets[dataset1][gene[1]][gene[2]:limit1])
                x.append(temp)
                wx.append(temp2)
                if not restrict and temp>np.max(data_sets[dataset1][gene[1]][gene[2]:gene[3]]):
                    w[n]='#E8053A'

            except:
                temp=np.max(data_sets[dataset1][gene[1]][gene[2]:])
                temp2=np.argmax(data_sets[dataset1][gene[1]][gene[2]:])
                x.append(temp)
                wx.append(temp2)
                if not restrict and temp>np.max(data_sets[dataset1][gene[1]][gene[2]:gene[3]]):
                    w[n]='#E8053A'
            try:
                temp=np.max(data_sets[dataset2][gene[1]][gene[2]:limit2])
                temp2=np.argmax(data_sets[dataset2][gene[1]][gene[2]:limit2])
                y.append(temp)
                wy.append(temp2)
                if not restrict and temp>np.max(data_sets[dataset2][gene[1]][gene[2]:gene[3]]):
                    w[n]='#E8053A'
            except:
                temp=np.max(data_sets[dataset2][gene[1]][gene[2]:])
                temp2=np.argmax(data_sets[dataset2][gene[1]][gene[2]:])
                y.append(temp)
                wy.append(temp2)
                if not restrict and temp>np.max(data_sets[dataset2][gene[1]][gene[2]:gene[3]]):
                    w[n]='#E8053A'
            z.append(gene[0])
        else:
            limit1=gene[3]-range1
            limit2=gene[3]-range2
            if restrict:
                if gene[3]-gene[2]<range1:
                    limit1=gene[2]
                if gene[3]-gene[2]<range2:
                    limit2=gene[2]
            try:
                temp=np.max(data_sets[dataset1][gene[1]][limit1:gene[3]])
                temp2=range1-np.argmax(data_sets[dataset1][gene[1]][limit1:gene[3]])
                x.append(temp)
                wx.append(temp2)
                if not restrict and temp>np.max(data_sets[dataset1][gene[1]][gene[2]:gene[3]]):
                    w[n]='#E8053A'
            except:
                temp=np.max(data_sets[dataset1][gene[1]][:gene[3]])
                temp2=range1-np.argmax(data_sets[dataset1][gene[1]][:gene[3]])
                x.append(temp)
                wx.append(temp2)
                if not restrict and temp>np.max(data_sets[dataset1][gene[1]][gene[2]:gene[3]]):
                    w[n]='#E8053A'
            try:
                temp=np.max(data_sets[dataset2][gene[1]][limit2:gene[3]])
                temp2=range2-np.argmax(data_sets[dataset2][gene[1]][limit2:gene[3]])
                y.append(temp)
                wy.append(temp2)
                if not restrict and temp>np.max(data_sets[dataset2][gene[1]][gene[2]:gene[3]]):
                    w[n]='#E8053A'
            except:
                temp=np.max(data_sets[dataset2][gene[1]][:gene[3]])
                temp2=range2-np.argmax(data_sets[dataset2][gene[1]][:gene[3]])
                y.append(temp)
                wy.append(temp2)
                if not restrict and temp>np.max(data_sets[dataset2][gene[1]][gene[2]:gene[3]]):
                    w[n]='#E8053A'

            z.append(gene[0])
    return x,y,z,w,wx,wy

def test(range):
    import matplotlib.pyplot as plt
    x,y,z,w,wx,wy=get_max('wt_me3','wt_me2',range,range)
    #plt.scatter(x,y)
    #plt.savefig('temp.png')
    return x,y,z,w,wx,wy

class pairwise(server.App):
    title="Pairwise Comparisons"
    tabs=['Comparison']

    inputs = [{	"type":'dropdown',
				"label": 'X-axis',
				"options" : [
					{"label": "Choose Dataset", "value":"empty"},
					{"label": "H3K4me3", "value":"wt_me3"},
					{"label": "H3K4me2", "value":"wt_me2"},
					{"label": "H3K4me3 spp1", "value":"spp1_me3"},
                    {"label": "H3K4me2 spp1", "value":"spp1_me2"},
                    {"label": "H3K4me3 swd2", "value":"swd2_me3"},
                    {"label": "H3K4me2 swd2", "value":"swd2_me2"},
                    {"label": "H3K4me3 set2", "value":"set2_me3"},
                    {"label": "H3K4me2 set2", "value":"set2_me2"},
                    {"label": "TBP", "value":"tbp"},
                    {"label": "H4Ac", "value":"h4ac"},
                    {"label": "H3", "value":"h3"},
                     {"label": "NCB2", "value":"ncb2"}],
				"key": 'ticker',
				"linked_key":'custom_ticker',
				"linked_type":'text'},
              {"input_type": 'slider',
               "label": 'Extend Range',
               "min": 0,
               "max": 2000,
               "value": 500,
               "variable_name": 'range1',
               },
              {	"type":'dropdown',
				"label": 'Y-axis',
				"options" : [
					{"label": "Choose Dataset", "value":"empty"},
					{"label": "H3K4me3", "value":"wt_me3"},
					{"label": "H3K4me2", "value":"wt_me2"},
					{"label": "H3K4me3 spp1", "value":"spp1_me3"},
                    {"label": "H3K4me2 spp1", "value":"spp1_me2"},
                    {"label": "H3K4me3 swd2", "value":"swd2_me3"},
                    {"label": "H3K4me2 swd2", "value":"swd2_me2"},
                    {"label": "H3K4me3 set2", "value":"set2_me3"},
                    {"label": "H3K4me2 set2", "value":"set2_me2"},
                    {"label": "TBP", "value":"tbp"},
                    {"label": "H4Ac", "value":"h4ac"},
                    {"label": "H3", "value":"h3"},
                    {"label": "NCB2", "value":"ncb2"}],
				"key": 'ticker2',
				"action_id": "update_data",
				"linked_key":'custom_ticker2',
				"linked_type":'text'},
              {"input_type": 'slider',
               "label": 'Extend Range',
               "min": 0,
               "max": 2000,
               "value": 500,
               "variable_name": 'range2',
               },
              {"input_type": 'checkboxgroup',
               "label":"<strong>Option</strong><br>(This option controls if the maximum value is restricted to only the"+
                        " size of the gene if the range is longer than that size"+
                       ".If not selected, genes in which the maximum is mapped outside of transcriptional"+
                  " unit will be colored red)",
               "options": [
                   {"label": "Restrict to gene lenght", "value": 1, "checked": False}],
               "variable_name": 'restrict'}
              ]

    controls = [{"control_type": "button",
                 "control_id": "button1",
                 "label": "Plot",
                }

                ]
    outputs = [{"type" : "html",
					"id" : "html_id",
					"control_id" : "button1",
					"tab" : "Comparison",
                "on_page_load" : False }]


    def getHTML(self,params):
        restriction=bool(params['restrict'])
        print (restriction)

        x,y,z,w,wx,wy=get_max(params['ticker'],params['ticker2'],int(params['range1']),int(params['range2']),restrict=restriction)
        data_sets_names = ['wt_me3', 'wt_me2','spp1_me3','spp1_me2',
                   'swd2_me3','swd2_me2','set2_me3','set2_me2','tbp','h4ac','h3','ncb2']
        #print(Counter(w))
        labels={'wt_me3':'WT H3K4me3',
                'wt_me2':'WT H3K4me2',
                'spp1_me3':'SPP1 H3K4me3',
                'spp1_me2':'SPP1 H3K4me2',
                'swd2_me3':'SWD2 H3K4me3',
                'swd2_me2':'SWD2 H3K4me2',
                'set2_me3':'SET2 H3K4me3',
                'set2_me2':'SET2 H3K4me2',
                'tbp':'TBP',
                'h4ac':'H4Ac',
                'h3':'H3',
                'ncb2':'NCB2'}

        source2=ColumnDataSource({'x': x, 'y': y,'z':z,'wx':wx,'wy':wy})
        p1=plotting.figure(tools='hover,box_zoom,reset,tap,save',
                           plot_width=400,plot_height=400,
                           toolbar_location='right',
                           x_axis_label = labels[params['ticker']],
                           y_axis_label = labels[params['ticker2']])


        url="http://www.yeastgenome.org/search?query=@z"
        taptool = p1.select(type=TapTool)
        taptool.callback = OpenURL(url=url)
        hover = p1.select(dict(type=HoverTool))
        hover.tooltips = [("Gene", "@z"),('Dataset1','@wx'),('Dataset2','@wy')]
        hover.mode = 'mouse'
        #colors = ["#%02x%02x%02x" % (r, g, 125) for r, g in zip(np.floor(150+x), np.floor(150+y))]
        p1.scatter('x','y',source=source2,color=w,line_alpha=0.1,
        size=8, alpha=0.3,name="mystuff")
        script, div = components(p1,CDN)
        html = "%s\n%s"%(script,div)
        return html

    def getCustomJS(self):
       return INLINE.js_raw[0]

    def getCustomCSS(self):
        with open('{}/custom_style.css'.format(PATH)) as style:
            previous=style.read()
        return previous+INLINE.css_raw[0]





if __name__ == '__main__':
    app=pairwise()
    app.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))