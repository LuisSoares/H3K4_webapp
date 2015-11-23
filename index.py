__author__ = 'Luis'
from spyre import server


class Index(server.App):
    title='Index'
    outputs = [{"output_type": "html",
                "output_id": "Index",
                "on_page_load": True}]
    def getHTML(self, params):

        return '''<head>
        <meta name="description" content="Histone methylation, Chip-Seq Plotter, Buratowski Lab">
    <META NAME="ROBOTS" CONTENT="INDEX, FOLLOW">
    </head>
    <style>
        a:link {color:black;
    text-decoration: none;
}

a:visited {
    text-decoration: none;
}

a:hover {
    text-decoration: underline;
}

a:active {
    text-decoration: underline;
}</style>
<body><h1>Description</h1><p style="font-size:16px;font-weigth:900"><span><img src="/static/images/857953-Compass.png" alt="COMPASS COMPLEX" style="width:304px;height:228px;float:right"></span>This website contains web apps developed
        and used in the Buratowski Lab for the visualization of ChIP-Seq results and Next Generation Sequencing.

        The datasets used were obtained by Luis Soares and YooJin Joo using standard
        Buratowski lab protocols. Sequencing was performed at the Harvard Univeristy Bauer Center
        Sequencing Facility.</p>
        <p style="font-size:16px;font-weigth:900">This website is entirely coded in python and all data analysis is also done with python.There are currently 3 tools available at the website:GenePlotter,GenePlotter2 and Comparisons. The difference between GenePlotter and GenePlotter2 relates to
        the way the plots are displayed, while the former produces a static .png image embedded (which you can save using the mouse right button), in the latter the plot is
        javascript embeded making it interactive but with no option to save it (there are workarounds but are too complicated to be stated here). The data in GenePlotter2
        is also downsampled to make it able to be loaded in a decent time line.</p>
        <h2>Contact and Suggestions:</h2>

        <p style="font-size:16px;font-weigth:900"><a href="http://www.sciencelint.org">Luis Soares</a></p></body>'''

    def getCustomCSS(self):
        with open('/home/lint78/web_site/custom_style.css') as style:
            return style.read()+'''\n .left-panel{display: none;}
        .right-panel{width:80%;margin: 5px}'''

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
	app = Index()
	app.launch()