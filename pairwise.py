__author__ = 'Luis'

from spyre import server
import codecs
import os
import pandas as pd
import numpy as np
from bokeh import plotting
from bokeh.models import OpenURL,Callback,Rect,ColumnDataSource,CustomJS,HoverTool,Range1d,PanTool,WheelZoomTool,TapTool,BoxSelectTool,BoxZoomTool
from bokeh.embed import components
from bokeh.resources import INLINE
from bokeh.resources import CDN


data=pd.read_csv('/home/lint78/web_site/temp.csv')

x=data['Max_Peak_1000_WT_H3K4me3']
y=data['Max_Peak_1000_WT_H3K4me2']
z=data['ID']
source2=ColumnDataSource({'x': x, 'y': y,'z':z})




class pairwise(server.App):
    title="Pairwise Comparisons"
    tabs=['Comparison']
    inputs = [{	"type":'dropdown',
				"label": 'X-axis',
				"options" : [
					{"label": "Choose Dataset", "value":"empty"},
					{"label": "H3K4me3", "value":"Max_Peak_1000_WT_H3K4me3"},
					{"label": "H3K4me2", "value":"Max_Peak_1000_WT_H3K4me2"},
					{"label": "H3K4me3 spp1", "value":"Max_Peak_1000_SPP1_H3K4me3"},
                    {"label": "H3K4me2 spp1", "value":"Max_Peak_1000_SPP1_H3K4me2"},
                    {"label": "H3K4me3 bre2", "value":"Max_Peak_1000_BRE2_H3K4me3"},
                    {"label": "H3K4me2 bre2", "value":"Max_Peak_1000_BRE2_H3K4me2"},
                    {"label": "H3K4me3 sdc1", "value":"Max_Peak_1000_SDC1_H3K4me3"},
                     {"label": "H3K4me2 sdc1", "value":"Max_Peak_1000_SDC1_H3K4me2"}],
				"key": 'ticker',
				"linked_key":'custom_ticker',
				"linked_type":'text'},
              {	"type":'dropdown',
				"label": 'Y-axis',
				"options" : [
					{"label": "Choose Dataset", "value":"empty"},
					{"label": "H3K4me3", "value":"Max_Peak_1000_WT_H3K4me3"},
					{"label": "H3K4me2", "value":"Max_Peak_1000_WT_H3K4me2"},
					{"label": "H3K4me3 spp1", "value":"Max_Peak_1000_SPP1_H3K4me3"},
                    {"label": "H3K4me2 spp1", "value":"Max_Peak_1000_SPP1_H3K4me2"},
                    {"label": "H3K4me3 bre2", "value":"Max_Peak_1000_BRE2_H3K4me3"},
                    {"label": "H3K4me2 bre2", "value":"Max_Peak_1000_BRE2_H3K4me2"},
                    {"label": "H3K4me3 sdc1", "value":"Max_Peak_1000_SDC1_H3K4me3"},
                     {"label": "H3K4me2 sdc1", "value":"Max_Peak_1000_SDC1_H3K4me2"}],
				"key": 'ticker2',
				"action_id": "update_data",
				"linked_key":'custom_ticker2',
				"linked_type":'text'}]

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
        x=data[params['ticker']]
        y=data[params['ticker2']]
        labels={'Max_Peak_1000_WT_H3K4me3':'WT H3K4me3',
                'Max_Peak_1000_WT_H3K4me2':'WT H3K4me2',
                'Max_Peak_1000_SPP1_H3K4me3':'SPP1 H3K4me3',
                'Max_Peak_1000_SPP1_H3K4me2':'SPP1 H3K4me2',
                'Max_Peak_1000_BRE2_H3K4me3':'BRE2 H3K4me3',
                'Max_Peak_1000_BRE2_H3K4me2':'BRE2 H3K4me2',
                'Max_Peak_1000_SDC1_H3K4me3':'SDC1 H3K4me3',
                'Max_Peak_1000_SDC1_H3K4me2':'SDC1 H3K4me2'}
        z=data['ID']
        source2=ColumnDataSource({'x': x, 'y': y,'z':z})
        p1=plotting.figure(tools='hover,box_zoom,reset,tap,save',
                           plot_width=400,plot_height=400,
                           toolbar_location='right',
                           x_axis_label = labels[params['ticker']],
                           y_axis_label = labels[params['ticker2']])
        #url = "http://h3k4.sciencelint.org/?freq=@z&check_boxes=__list__,0,1,2&range=__float__0&"
        url="http://www.yeastgenome.org/search?query=@z"
        taptool = p1.select(type=TapTool)
        taptool.callback = OpenURL(url=url)
        hover = p1.select(dict(type=HoverTool))
        hover.tooltips = [("Gene", "@z")]
        hover.mode = 'mouse'
        colors = ["#%02x%02x%02x" % (r, g, 125) for r, g in zip(np.floor(150+x), np.floor(150+y))]
        p1.scatter('x','y',source=source2,
        radius=1, alpha=0.2,color=colors,name="mystuff")
        script, div = components(p1,CDN)
        html = "%s\n%s"%(script,div)
        return html

    def getCustomJS(self):
       return INLINE.js_raw[0]

    def getCustomCSS(self):
        with open('/home/lint78/web_site/custom_style.css') as style:
            previous=style.read()
        return previous+INLINE.css_raw[0]





if __name__ == '__main__':
    app=pairwise()
    app.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))
