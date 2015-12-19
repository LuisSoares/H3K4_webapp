#!/usr/bin/python
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib as mpl

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

def parser():
    genes1=[]
    datasets1=[]
    data=open('/home/lint78/web_site/log.txt')
    for line in data:
        temp=line.split('\t')
        genes1.append(eval(temp[0])[0])
        datasets1.extend(eval(temp[1]))
    data.close()
    return Counter(genes1).most_common(10),Counter(datasets1).most_common(12)

def plot_stats():
    colors_pie = ['yellowgreen', 'gold',
              'lightskyblue', 'lightcoral',
              'lightgreen','lightsalmon',
              'palevioletred','lightslategray',
              'khaki','wheat','crimson','cornsilk']
    names_of_sets=['H3K4me$^3$','H3K4me$^2$',
                   'H3K4me$^3$($\Delta$spp1)',
                   'H3K4me$^2$($\Delta$spp1)',
                   'H3K4me$^3$($\Delta$swd2)',
                   'H3K4me$^2$($\Delta$swd2)',
                   'H3K4me$^3$($\Delta$set2)',
                   'H3K4me$^2$($\Delta$set2)',
                   'TBP','H4Ac','H3','NCB2']
    genes1,datasets1=parser()
    names1,counts1=zip(*genes1)
    sets1,set_counts1=zip(*datasets1)
    fig1,(ax11,ax21)=plt.subplots(1,2,figsize=(9,5))
    fig1.subplots_adjust(bottom=0.35,right=0.75)
    ax11.bar(range(len(names1)),counts1,width=0.95,linewidth=0,color='orange')
    ax11.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5])
    ax11.set_xticklabels(names1,rotation='vertical')
    ax11.set_ylim(0,max(counts1)*1.2)
    ax11.set_xlabel('Genes')
    ax11.set_ylabel('# of searches')
    ax11.set_title('Top 10 searched genes')
    ax21.set_title('Searched Datasets')
    ax21.pie(set_counts1,explode=[0.05 for item in range(len(set_counts1))],colors=colors_pie,shadow=True)
    ax21.legend([names_of_sets[int(x)] for x in sets1],bbox_to_anchor=(1.7, 1),fontsize=10)
    plt.savefig('/home/lint78/web_site/log.png')


plot_stats()