__author__ = 'Luis'
from intermine.webservice import Service
from spyre import server
import codecs
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import h5py
import jinja2
import json
import mpld3
from mpld3 import plugins, utils
import cStringIO
from scipy import ndimage
number_of_samples=0
#Setting plotting parameters
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)
def create_java(n):
        java_string='''obj[%s].elements().
            on("mouseover", function(d, i){
            d3.select(this).transition().duration(50).style("stroke-opacity", alpha_fg);
            obj[%s].elements().transition().duration(200).style("stroke-opacity", alpha_fg);})
        .on("mouseout", function(d, i){
            d3.select(this).transition().duration(200).style("stroke-opacity", alpha_bg);
             obj[%s].elements().transition().duration(200).style("stroke-opacity", alpha_bg);});'''
        full_string=''
        count=0
        while count<n:
            full_string=full_string+(java_string)%(count,count+1,count+1)
            count=count+2
        return full_string





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

#server.View.View = CustomView


data_sets=h5py.File("/home/lint78/web_site/datasets.hdf5","r")
data_sets_names = ['wt_me3', 'wt_me2','spp1_me3','spp1_me2',
                   'swd2_me3','swd2_me2','set2_me3','set2_me2','tbp','h4ac','h3','ncb2']
location_sets=['genes','CUTS','SUTS']
data_sets_legend = ['H3K4me3', 'H3K4me2','H3K4me3 (spp1)','H3K4me2 (spp1)'
    ,'H3K4me3 (swd2)','H3K4me2 (swd2)','H3K4me3 (set2)'
    ,'H3K4me2 (set2)','TBP','H4Ac','H3','NCB2']

class GenePlotter2(server.App):
    title = "Gene plotter2"
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

    outputs = [{"output_type": "html",
                "control_id": "button1",
                "output_id": "sine_wave_plot",
                "on_page_load": False,
                'tab': 'Plot'}]

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


    def getHTML(self, params):
        #print(params)
        global current_gene
        data_max=0
        fig= plt.figure(figsize=(9, 7))
        gs = gridspec.GridSpec(4, 1, height_ratios=[6, 1,1,1])
        gene = params['freq'].upper()
        range_ext = int(params['range'])
        searched='False'
        lines=[]
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
        print('*'*40+'\nQuery by user:'+str(current_gene)+';\nrange extension:'
            +str(range_ext)+' ;Checked Boxes:'+str(params['check_boxes'])+';SGD used:'+searched+'\n'+'*'*40)
        #print(params)
        splt1 = plt.subplot(gs[0, :])
        splt1.set_title('{}-Chromosome {}:{}:{}, strand:{}'.format(current_gene[0], current_gene[1], current_gene[2],
                                                                    current_gene[3], current_gene[4])
                        ,size=20)
        # splt1.set_xlabel('Distance from TSS (nt)')
        splt1.axes.xaxis.set_ticklabels([])
        splt1.set_ylabel('Normalized counts per million reads\n ',size=16,labelpad=15)
        n=0
        fill_colors=['#921B16', '#D6701C','#251F47','#68A691','#68A691',
                     '#447604','#232621','#0BBE30','#BC3908','#191923','#391923',
                     '#591923']
        for item in params['check_boxes']:
            if current_gene[4] == '+':
                data=data_sets[data_sets_names[int(item)]][current_gene[1]][current_gene[2] - range_ext:current_gene[3] + range_ext]
                factor=100./len(data)

                data=ndimage.interpolation.zoom(data,factor)

                _=splt1.plot(data,
                    label=data_sets_legend[int(item)],color=fill_colors[n],alpha=0.3,lw=3)
                lines.extend(_)
                splt1.fill_between(range(0,len(data)),0, data,alpha=0.2,color=fill_colors[n])
            if current_gene[4] == '-':
                data=data_sets[data_sets_names[int(item)]][current_gene[1]][current_gene[2] - range_ext:current_gene[3] + range_ext]
                factor=100./len(data)

                data=ndimage.interpolation.zoom(data[::-1],factor)
                print(data)
                _=splt1.plot(data,
                    label=data_sets_legend[int(item)],color=fill_colors[n],alpha=0.3,lw=3)
                lines.extend(_)
                splt1.fill_between(range(0,len(data)),0, data,alpha=0.2,color=fill_colors[n])
            if np.max(data)>data_max:data_max=np.max(data)
            n+=1
        #plot2
        if current_gene[4] == '+':
            y_values1 = data_sets['genes'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext]
            y_values2 = data_sets['genes'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext]
        else:
            y_values1 =  data_sets['genes'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
            y_values2 =  data_sets['genes'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]

        #splt2.plot(y_values, 'grey', linewidth=2)
        factor=100./len(y_values1)
        splt2 = plt.subplot(gs[1, :])
        y_values1=ndimage.interpolation.zoom(y_values1,factor)
        y_values2=ndimage.interpolation.zoom(y_values2,factor)
        x_values = range(0, len(y_values1))
        splt2.plot([0, len(y_values1)], [0, 0], 'black', linewidth=1.5)
        splt2.fill_between(x_values, y_values1/4.0, y_values1, color='#FFB30F')
        splt2.fill_between(x_values, y_values2/4.0, y_values2, color='#FD151B')
        splt2.set_ylim(-1.5, 1.5)
        splt2.text(1.01, 0.5, 'GENES', fontsize=12, transform=splt2.transAxes, verticalalignment='center')
        splt2.yaxis.set_tick_params(color='w')
        splt2.set_yticks([-0.5, 0.5])
        splt2.spines['left'].set_position(('outward', 10))
        splt2.set_yticklabels([])
        splt2.axes.xaxis.set_ticklabels([])
        #plot3

        if current_gene[4] == '+':
            y_values1 = data_sets['CUTS'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext]
            y_values2 = data_sets['CUTS'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext]
        else:
            y_values1 =  data_sets['CUTS'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
            y_values2 =  data_sets['CUTS'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]

        #splt2.plot(y_values, 'grey', linewidth=2)
        factor=100./len(y_values1)
        splt3 = plt.subplot(gs[2, :])
        y_values1=ndimage.interpolation.zoom(y_values1,factor)
        y_values2=ndimage.interpolation.zoom(y_values2,factor)
        x_values = range(0, len(y_values1))
        splt3.plot([0, len(y_values1)], [0, 0], 'black', linewidth=1.5)
        splt3.fill_between(x_values, y_values1/4.0, y_values1, color='#028090')
        splt3.fill_between(x_values, y_values2/4.0, y_values2, color='#6EEB83')
        splt3.set_ylim(-1.5, 1.5)
        splt3.text(1.01, 0.5, 'CUTS', fontsize=12, transform=splt3.transAxes, verticalalignment='center')
        splt3.yaxis.set_tick_params(color='w')
        splt3.set_yticks([-0.5, 0.5])
        splt3.spines['left'].set_position(('outward', 10))
        splt3.set_yticklabels([])
        splt3.axes.xaxis.set_ticklabels([])
        #plot4

        if current_gene[4] == '+':
            y_values1 = data_sets['SUTS'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext]
            y_values2 = data_sets['SUTS'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext]
        else:
            y_values1 =  data_sets['SUTS'][current_gene[1]][0,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]
            y_values2 =  data_sets['SUTS'][current_gene[1]][1,current_gene[2] - range_ext:current_gene[3] + range_ext][::-1]

        #splt2.plot(y_values, 'grey', linewidth=2)
        factor=100./len(y_values1)
        splt4 = plt.subplot(gs[3, :])
        y_values1=ndimage.interpolation.zoom(y_values1,factor)
        y_values2=ndimage.interpolation.zoom(y_values2,factor)
        x_values = range(0, len(y_values1))
        splt4.plot([0, len(y_values1)], [0, 0], 'black', linewidth=1.5)
        splt4.fill_between(x_values, y_values1/4.0, y_values1, color='#028090')
        splt4.fill_between(x_values, y_values2/4.0, y_values2, color='#6EEB83')
        splt4.set_ylim(-1.5, 1.5)
        splt4.text(1.01, 0.5, 'SUTS', fontsize=12, transform=splt4.transAxes, verticalalignment='center')
        splt4.yaxis.set_tick_params(color='w')
        splt4.set_yticks([-0.5, 0.5])
        splt4.spines['left'].set_position(('outward', 10))
        splt4.set_yticklabels([])
        splt4.axes.xaxis.set_ticklabels([])
        leg=splt1.legend(ncol=3)
        leg.set_title("")
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)
        size = (current_gene[3] + range_ext) - (current_gene[2] - range_ext)
        splt1.set_ylim(0, np.max(data_max)*1.3)
        #print(str(plt.ylim()[1] - 5))
        splt1.text(0.02, 0.02, 'Luis Soares', fontsize=6, transform=splt1.transAxes, verticalalignment='top')

        leg_line=leg.get_lines()
        number_of_samples=len(leg_line)
        class HighlightLines(plugins.PluginBase):
            """A plugin to highlight lines on hover"""

            JAVASCRIPT = """
    mpld3.register_plugin("linehighlight", LineHighlightPlugin);
    LineHighlightPlugin.prototype = Object.create(mpld3.Plugin.prototype);
    LineHighlightPlugin.prototype.constructor = LineHighlightPlugin;
    LineHighlightPlugin.prototype.requiredProps = ["legend_ids","line_ids"];
    LineHighlightPlugin.prototype.defaultProps = {alpha_bg:0.3, alpha_fg:1.0}
    function LineHighlightPlugin(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    LineHighlightPlugin.prototype.draw = function(){
      var obj=[]
      for(var i=0; i<this.props.legend_ids.length; i++){
         obj=obj.concat([mpld3.get_element(this.props.legend_ids[i], this.fig),
         mpld3.get_element(this.props.line_ids[i], this.fig)]);

        alpha_fg = this.props.alpha_fg;
        alpha_bg = this.props.alpha_bg;
}
%s

        }
    """%(create_java(number_of_samples*2))


            def __init__(self, legend,lines,samples):
                self.lines = lines
                self.legend=legend
                self.dict_ = {"type": "linehighlight",
                      "legend_ids": [utils.get_id(line) for line in legend],
                      "line_ids": [utils.get_id(line) for line in lines],
                      "alpha_bg": lines[0].get_alpha(),
                      "alpha_fg": 1.0}
        class TopToolbar(plugins.PluginBase):
            """Plugin for moving toolbar to top of figure"""

            JAVASCRIPT = """
    mpld3.register_plugin("toptoolbar", TopToolbar);
    TopToolbar.prototype = Object.create(mpld3.Plugin.prototype);
    TopToolbar.prototype.constructor = TopToolbar;
    function TopToolbar(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    TopToolbar.prototype.draw = function(){
      // the toolbar svg doesn't exist
      // yet, so first draw it
      this.fig.toolbar.draw();

      // then change the y position to be
      // at the top of the figure
      this.fig.toolbar.toolbar.attr("y", 2);

      // then remove the draw function,
      // so that it is not called again
      this.fig.toolbar.draw = function() {}
    }
    """
            def __init__(self):
                self.dict_ = {"type": "toptoolbar"}

        plugins.clear(fig)
        plugins.connect(fig,
                        TopToolbar(),
                        HighlightLines(leg_line,lines,number_of_samples)
                        ,plugins.Reset()
                        ,plugins.BoxZoom())
        file = cStringIO.StringIO()
        mpld3.save_html(fig,file)
        return file.getvalue()




    def getTable(self, params):
        df = pd.DataFrame(columns=['Track', 'Code', 'Chromosome', 'Start', 'End', 'Strand', 'Max', 'Max Pos'])

        gene = params['freq'].upper()
        print(gene)
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
    app = GenePlotter2()
    app.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))
    # site.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))
    # app.launch(host='0.0.0.0')
