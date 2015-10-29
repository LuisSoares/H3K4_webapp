__author__ = 'Luis'
from spyre import server
import os
from index import Index
from pairwise import pairwise
from simple_spyre import SimpleApp
from references import References


site = server.Site(Index)


site.root.templateVars['custom_head'] ='''<div class="banner">
		<h1 style="font-family: Helvetica, sans-serif;font-size: 30px
		;margin-left:2%;margin-top:2%">Buratowski Lab Yeast H3K4 methylation studies</h1></div>'''

site.addApp(SimpleApp, '/app2')
site.addApp(pairwise, '/app3')
site.addApp(References, '/app4')
for fullRoute, _ in site.site_app_bar[1:]:
			parent, route = site.get_route(fullRoute)
			parent.__dict__[route].templateVars['custom_head'] ='''<div class="banner">
		<h1 style="font-family: Helvetica, sans-serif;font-size: 30px
		;margin-left:2%;margin-top:2%">Buratowski Lab Yeast H3K4 methylation studies</h1></div>'''


if __name__ == '__main__':
	site.launch(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')))
