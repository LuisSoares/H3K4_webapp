__author__ = 'Luis'
from spyre import server


class References(server.App):
    title='References'
    tabs=['Experimental','Informatics']
    outputs = [{"output_type": "html",
                "id": "Experimental",
                "on_page_load": True},
               {"output_type": "html",
                "id": "Informatics",
                "on_page_load": True}]

    def Experimental(self, params):

        return ''''''

    def Informatics(self, params):
        with open('informatics.html') as page:
            return page.read()

    def getCustomCSS(self):
        with open('custom_style.css') as style:
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
	app = References()
	app.launch()