# Adapted from bioinfokit
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
from random import sample
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import seaborn as sns
from adjustText import adjust_text

plt.rcParams.update({'font.family':'sans-serif'})
plt.rcParams.update({'font.sans-serif':'Arial'})
    

class general:
    def __init__(self):
        pass

    rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                   '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                   '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                   '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')

    def get_figure(show, r, figtype, fig_name):
        if show:
            plt.show()
        else:
            plt.savefig(fig_name+'.'+figtype, format=figtype, bbox_inches='tight', dpi=r)
        #plt.close()

    def axis_labels(x, y, axlabelfontsize=None, axlabelfontname=None):
        plt.xlabel(x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.ylabel(y, fontsize=axlabelfontsize, fontname=axlabelfontname)

    def axis_ticks(xlm=None, ylm=None, axtickfontsize=None, axtickfontname=None, ar=None):
        if xlm:
            plt.xlim(left=xlm[0], right=xlm[1])
            plt.xticks(np.arange(xlm[0], xlm[1], xlm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.xticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.yticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

    def depr_mes(func_name):
        print("This function is deprecated. Please use", func_name )
        print("Read docs at https://reneshbedre.github.io/blog/howtoinstall.html")

    def check_for_nonnumeric(pd_series=None):
        if pd.to_numeric(pd_series, errors='coerce').isna().sum() == 0:
            return 0
        else:
            return 1


    
class gene_exp:
    def __init__(self):
        pass
    def geneplot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle, plotlabelrotation):
        texts = list()
        if genenames is not None and genenames == "deg":
            for i in d[geneid].unique():
                if (d.loc[d[geneid] == i, lfc].iloc[0] >= lfc_thr and d.loc[d[geneid] == i, pv].iloc[0] < pv_thr) or \
                        (d.loc[d[geneid] == i, lfc].iloc[0] <= -lfc_thr and d.loc[d[geneid] == i, pv].iloc[0] < pv_thr):
                    if gstyle==1:
                        texts.append(plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0], i, fontsize=gfont, rotation = plotlabelrotation))
                    elif gstyle==2:
                        plt.annotate(i, xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        elif genenames is not None and type(genenames) is tuple:
            for i in d[geneid].unique():
                if i in genenames:
                    if gstyle==1:
                        texts.append(plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0], i, fontsize=gfont, rotation = plotlabelrotation))
                    elif gstyle==2:
                        plt.annotate(i, xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        elif genenames is not None and type(genenames) is dict:
            for i in d[geneid].unique():
                if i in genenames:
                    if gstyle==1:
                        texts.append(plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0],
                                      genenames[i], fontsize=gfont, rotation = plotlabelrotation))
                    elif gstyle == 2:
                        plt.annotate(genenames[i], xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        adjust_text(texts, arrowprops=dict(arrowstyle='-|>, head_width=.1', color="#A4C2EA", shrinkA = 3, shrinkB = 3))


    def volcano(df="dataframe", lfc=None, pv=None, edgecolorcategory = 0, lfc_thr=1, pv_thr=0.05, 
                color=("green", "grey", "red"), valpha=1, figtitle = 'Volcano Plot',
                geneid=None, genenames=None, gfont=8, dim=(5, 5), r=300, ar=90, dotsize=8, markerdot="o",
                sign_line=False, gstyle=1, show=False, figtype='png', axtickfontsize=9,
                axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
                axylabel=None, xlm=None, ylm=None, plotlegend=False, legendpos='best',
                figname='volcano', legendanchor=None, legendlabels=['significant up', 'not significant', 'significant down'],
                plotlabelrotation = 0):
        _x = r'$ log_{2}(Fold Change)$'
        _y = r'$ -log_{10}(P-value)$'
        color = color
        # check if dataframe contains any non-numeric character
        assert general.check_for_nonnumeric(df[lfc]) == 0, 'dataframe contains non-numeric values in lfc column'
        assert general.check_for_nonnumeric(df[pv]) == 0, 'dataframe contains non-numeric values in pv column'
        # this is important to check if color or logpv exists and drop them as if you run multiple times same command
        # it may update old instance of df
        df = df.drop(['color_add_axy', 'logpv_add_axy'], axis=1, errors='ignore')
        #assert len(set(color)) == 3, 'unique color must be size of 3'
        df.loc[(df[lfc] >= lfc_thr) & (df[pv] < pv_thr), 'color_add_axy'] = color[0]  # upregulated
        
        
        df.loc[(df[lfc] <= -lfc_thr) & (df[pv] < pv_thr), 'color_add_axy'] = color[2]  # downregulated
        df['color_add_axy'].fillna(color[1], inplace=True)  # intermediate
        df['logpv_add_axy'] = -(np.log10(df[pv]))
        # print(df[df['color']==color[0]].count(), 'zzzz')
            
        # plot
        assign_values = {col: i for i, col in enumerate(color)}
        color_result_num = [assign_values[i] for i in df['color_add_axy']]
        fig, ax = plt.subplots(figsize=dim)
        if plotlegend:
            s = plt.scatter(df[lfc], df['logpv_add_axy'], c=color_result_num, cmap=ListedColormap(color), alpha=valpha, s=dotsize,
                    marker=markerdot)
            if edgecolorcategory !=0:
                temp = df[df.color_add_axy == color[0]] 
                temp = temp.append(df[df.color_add_axy == color[2]], ignore_index = True)
                
                colorcycle = pl.cm.hsv(np.linspace(0,1,temp[edgecolorcategory].value_counts().shape[0]))
                colorcycle[:,:3]*=.9
                colorcycle = list(colorcycle)
                ax.set_prop_cycle('color', colorcycle)
                groups = temp.groupby(edgecolorcategory)
                for name, group in groups:
                    ax.plot(group[lfc], group.logpv_add_axy, marker=markerdot, linestyle='', ms = dotsize, label=name)
                    
                ax.legend()
            else: 
                plt.legend(handles=s.legend_elements()[0], labels=legendlabels, loc=legendpos, bbox_to_anchor=legendanchor)
        else:
            plt.scatter(df[lfc], df['logpv_add_axy'], c=color_result_num, cmap=ListedColormap(color), alpha=valpha, s=dotsize,
                        marker=markerdot)
            
        plt.title(figtitle);
        if sign_line:
            plt.axhline(y=-np.log10(pv_thr), linestyle='--', color='#7d7d7d', linewidth=1)
            plt.axvline(x=lfc_thr, linestyle='--', color='#7d7d7d', linewidth=1)
            plt.axvline(x=-lfc_thr, linestyle='--', color='#7d7d7d', linewidth=1)
        gene_exp.geneplot(df, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle, plotlabelrotation)

        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        return fig, ax
        #general.get_figure(show, r, figtype, figname)
        
class marker:
    def __init__(self):
        pass

    def geneplot_mhat(df, markeridcol, chr, pv, gwasp, markernames, gfont, gstyle, ax, plotlabelrotation):
        if markeridcol is not None:
            if markernames is not None and markernames is True:
                for i in df[markeridcol].unique():
                    if df.loc[df[markeridcol] == i, pv].iloc[0] <= gwasp:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0], str(i), fontsize=gfont, rotation = plotlabelrotation)
                        elif gstyle == 2:
                            plt.annotate(i, xy=(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
            elif markernames is not None and isinstance(markernames, (tuple, list)):
                for i in df[markeridcol].unique():
                    if i in markernames:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0], str(i), fontsize=gfont, rotation = plotlabelrotation)
                        elif gstyle == 2:
                            plt.annotate(i, xy=(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
            elif markernames is not None and isinstance(markernames, dict):
                for i in df[markeridcol].unique():
                    if i in markernames:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0], markernames[i], fontsize=gfont, rotation = plotlabelrotation)
                        elif gstyle == 2:
                            plt.annotate(markernames[i], xy=(
                            df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
        else:
            raise Exception("provide 'markeridcol' parameter")
            
    def geneplot_mhat_logp(df, markeridcol, chr, tpval, gwasp, markernames, gfont, gstyle, ax, plotlabelrotation, vertalign = 'bottom'):
        loggwasp = np.log10(gwasp)
        if markeridcol is not None:
            if markernames is not None and markernames is True:
                for i in df[markeridcol].unique():
                    if abs(df.loc[df[markeridcol] == i, tpval].iloc[0]) >= abs(loggwasp):
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0],
                                    str(i), fontsize=gfont, rotation = plotlabelrotation, va = vertalign)
                        elif gstyle == 2:
                            plt.annotate(i, xy=(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
            elif markernames is not None and isinstance(markernames, (tuple, list)):
                for i in df[markeridcol].unique():
                    if i in markernames:
                        if gstyle == 1:
                           plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0],
                                str(i), fontsize=gfont, rotation = plotlabelrotation)
                        elif gstyle == 2:
                            plt.annotate(i, xy=(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
            elif markernames is not None and isinstance(markernames, dict):
                for i in df[markeridcol].unique():
                    if i in markernames:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0],
                                 markernames[i], fontsize=gfont, rotation = plotlabelrotation)
                        elif gstyle == 2:
                            plt.annotate(markernames[i], xy=(
                            df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
        else:
            raise Exception("provide 'markeridcol' parameter")
    
    def mhat(df="dataframe", 
             chr=None, 
             pv=None, 
             color=None, 
             dim=(6,4), 
             r=300, 
             ar=90, 
             gwas_sign_line=False,
             gwasp=5E-08, 
             dotsize=8, 
             markeridcol=None, 
             markernames=None, 
             gfont=8, 
             valpha=1, 
             show=False, 
             figtype='png',
             axxlabel=None, 
             axylabel=None, 
             axlabelfontsize=9, 
             axlabelfontname="Arial", 
             axtickfontsize=9, 
             figtitle='manhattan plot',
             axtickfontname="Arial", 
             ylm=None, 
             gstyle=1, 
             yskip = 1, 
             plotlabelrotation=0, 
             figname='manhattan', 
             invert=False, 
             fig=None, 
             ax=None, 
             xtickname=False):

        _x, _y = 'Chromosomes', r'$ -log_{10}(P)$'
        rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                       '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                       '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                       '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')
        '''
         rand_colors = ('#f67280', '#00a8cc', '#ffd082', '#fb8d62', '#6e5773', '#21bf73', '#d5c455', '#c9753d',
                       '#ad62aa','#d77fa1', '#a0855b', '#ffd800', '#da2d2d', '#6f9a8d', '#a8ff3e', '#b2fcff',
                       '#a0c334', '#b5525c', '#c06c84', '#3a3535', '#9b45e4', '#f6da63', '#9dab86', '#0c093c',
                       '#f6f078', '#64c4ed', '#da4302', '#5edfff', '#08ffc8', '#ca3e47', '#f7ff56', '#6c5ce7')
        '''
        # minus log10 of P-value
        if invert:
            df['tpval'] = np.log10(df[pv])
        else: 
            df['tpval'] = -np.log10(df[pv])
        df = df.sort_values(chr)
        # add indices
        df['ind'] = range(len(df))
        df_group = df.groupby(chr)
        if color is not None and len(color) == 2:
            color_1 = int(df[chr].nunique() / 2) * [color[0]]
            color_2 = int(df[chr].nunique() / 2) * [color[1]]
            if df[chr].nunique() % 2 == 0:
                color_list = list(reduce(lambda x, y: x+y, zip(color_1, color_2)))
            elif df[chr].nunique() % 2 == 1:
                color_list = list(reduce(lambda x, y: x+y, zip(color_1, color_2)))
                color_list.append(color[0])
        elif color is not None and len(color) == df[chr].nunique():
            color_list = color
        elif color is None:
            # select colors randomly from the list based in number of chr
            # color_list = sample(rand_colors, df[chr].nunique())
            color_list = rand_colors[:df[chr].nunique()]
        else:
            print("Error: in color argument")
            sys.exit(1)

        xlabels = []
        xticks = []
        if fig is None:
            fig, ax = plt.subplots(figsize=dim)
        i = 0
        for label, df1 in df.groupby(chr):
            df1.plot(kind='scatter', x='ind', y='tpval', color=color_list[i], s=dotsize, alpha=valpha, ax=ax)
            df1_max_ind = df1['ind'].iloc[-1]
            df1_min_ind = df1['ind'].iloc[0]
            xlabels.append(label)
            xticks.append((df1_max_ind - (df1_max_ind - df1_min_ind) / 2))
            i += 1

        # add GWAS significant line
        if gwas_sign_line is True:
            ax.axhline(y=-np.log10(gwasp), linestyle='--', color='#7d7d7d', linewidth=1)
        if markernames is not None:
            marker.geneplot_mhat(df, markeridcol, chr, pv, gwasp, markernames, gfont, gstyle, ax=ax, plotlabelrotation =  plotlabelrotation)
        ax.margins(x=0)
        ax.margins(y=0)
        ax.set_xticks(xticks)
        ax.set_ylim([0, max(df['tpval'] + 1)+10])
        if ylm:
            ylm = np.arange(ylm[0], ylm[1], ylm[2])
        else:
            ylm = np.arange(0, max(df['tpval']+1), yskip)
        ax.set_yticks(ylm)
        if xtickname:
            ax.set_xticklabels(map(ICDchapter_to_name,xlabels), 
                               rotation=ar, 
                               va='top',
                               ha='right',
                               fontsize=axtickfontsize)
        else: 
            ax.set_xticklabels(xlabels, rotation=ar, fontsize=axtickfontsize)
        ax.set_yticklabels(ylm, fontsize=axtickfontsize, fontname=axtickfontname, rotation=ar)
        
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        ax.set_xlabel(_x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        ax.set_ylabel(_y, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.title(figtitle, fontsize=axlabelfontsize+2)
        #general.get_figure(show, r, figtype, figname)
        return fig, ax
    
    def mhat_RE(df="dataframe", 
                 chromo=None, 
                 logp=None, 
                 color=None, 
                 dim=(10,10), 
                 rows=None,
                 columns=None,
                 nrowstop=None, # number of rows for top subplot
                 nrowsmid=None, # number of rows for middle subplot
                 topmin=None, # min y-axis value for top subplot
                 topmax=None, # max y-axis value for top subplot
                 mainmin=None, # min y-axis value for main subplot
                 mainmax=None, # max y-axis value for main subplot
                 r=300, 
                 ar=90, 
                 gwas_sign_line=False,
                 gwasp=5E-08, 
                 dotsize=8, 
                 markeridcol=None, 
                 markernames=None, 
                 gfont=8, 
                 valpha=1, 
                 show=False, 
                 figtype='pdf',
                 axxlabel=None, 
                 axylabel=None, 
                 axlabelfontsize=9, 
                 axlabelfontname="Arial", 
                 axtickfontsize=9, 
                 figtitle='manhattan plot',
                 textcolor='k', # text annotation color 
                 axtickfontname="Arial", 
                 ylm=None, 
                 gstyle=1, 
                 yskip=1, 
                 plotlabelrotation=0, 
                 figname='miami', 
                 invert=False,       
                 fig=None, 
                 ax=None, 
                 icd10_mapping=codemap3,
                 annotate=True,
                 annotatefontsize=20,
                 expand_text=(1.05, 1.2), # expand_text, expand_points, and autoalign set to adjust_text's default values
                 expand_points=(1.05, 1.2),
                 autoalign=True,
                 invisible_ticks=False,
                 suffix=None,
                 axescolors='k',
                 overlap=None): # phenotypes in overlap set will be bolded to indicate common top phenotypes among all groups
        
        # 20210823 Added icd10_mapping parameter to add in the icd-10 names

        # _y denotes what will be in y-axis
        _x, _y = 'Chromosomes', r'$ -log_{10}(P)$'
        
        # don't annotate y-axis for Fig 1, where there will be no annotation of text
        if annotate == False:
            _y = None

        # tpval1 corresponds to -log10_pvalue of first group
        df['tpval'] = df[logp]
        df = df.sort_values(chromo)

        df['ind'] = range(len(df))
        df_group = df.groupby(chromo)

        rand_colors = ('#a7414a', 
                       '#282726', 
                       '#6a8a82', 
                       '#a37c27', 
                       '#563838', 
                       '#0584f2', 
                       '#f28a30', 
                       '#f05837',
                       '#6465a5', 
                       '#00743f', 
                       '#be9063', 
                       '#de8cf0', 
                       '#888c46', 
                       '#c0334d', 
                       '#270101', 
                       '#8d2f23',
                       '#ee6c81', 
                       '#65734b', 
                       '#14325c', 
                       '#704307', 
                       '#b5b3be', 
                       '#f67280', 
                       '#ffd082', 
                       '#ffd800',
                       '#ad62aa', 
                       '#21bf73', 
                       '#a0855b', 
                       '#5edfff', 
                       '#08ffc8', 
                       '#ca3e47', 
                       '#c9753d', 
                       '#6c5ce7')
        
        #color_list = sample(rand_colors, df[chromo].nunique())
        color_list = rand_colors[:df[chromo].nunique()]

        xlabels = []
        xticks = []
        
        if fig is None:
            fig = plt.figure(figsize=dim)
            fig.tight_layout()
            
        rows = rows
        columns = columns

        # grid0 for y axis, grid1 for miami plot
        grid0 = plt.GridSpec(rows, columns, left=0.50, right=0.55) #  wspace = .25, hspace = .25, 
        grid1 = plt.GridSpec(rows, columns, hspace=0)
        
        ax0 = plt.subplot(grid0[:, 0])
        ax0.axis('off')

        i=0
        texts = []
        kwargs1 = dict(color=textcolor, fontweight='medium')
        kwargs2 = dict(color='k', fontweight='heavy')
        for label, df1 in df.groupby(chromo):
            ax = plt.subplot(grid1[(nrowstop+1):, 1])
            ax.scatter(df1['ind'], df1['tpval'], color=color_list[i], s=dotsize)
            
            
            # Add annotation (20220111)
            # source: https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe/39374693
            if annotate:
                df1 = df1.set_index('phenotype')
                for idx, row in df1.iterrows():
                    if row['annotate'+suffix] == 1:
                        if idx in overlap:
                            texts.append(ax.annotate(idx, (row['ind'], row['tpval']), fontsize=annotatefontsize, **kwargs2))
                        else:
                            texts.append(ax.annotate(idx, (row['ind'], row['tpval']), fontsize=annotatefontsize, **kwargs1))
                        
                    
            d = .007  # how big to make the diagonal lines in axes coordinates
            # arguments to pass to plot, just so we don't keep repeating them
            kwargs = dict(transform=ax.transAxes, color=axescolors, clip_on=False, linewidth=1)
            ax.plot((-d, +d), (1, 1), **kwargs)  # bottom-left diagonal
            ax.plot((1 - d, 1 + d), (1, 1), **kwargs)  # bottom-right diagonal
           
            # Add line breaks (20220112)
            # Plot the same data above on both additional axes
            
            # ax2 is top subplot
            ax2 = plt.subplot(grid1[0:nrowstop, 1])
            ax2.scatter(df1['ind'], df1['tpval'], color=color_list[i], s=dotsize)
            
            
            kwargs.update(transform=ax2.transAxes)
            ax2.plot((-d, +d), (0, 0), **kwargs)        # top-left diagonal
            ax2.plot((1 - d, 1 + d), (0, 0), **kwargs)  # top-right diagonal
           
            df1_max_ind = df1['ind'].iloc[-1]
            df1_min_ind = df1['ind'].iloc[0]
            xlabels.append(label)
            xticks.append((df1_max_ind - (df1_max_ind - df1_min_ind) / 2))
            
            i += 1
       
        adjust_text(texts,
                    arrowprops=dict(arrowstyle="-", 
                                    color='k'),
                    autoalign=autoalign,
                    expand_text=expand_text,
                    expand_points=expand_points,
                    ax=ax)
            
        ax.axhline(y=0, color='#7d7d7d', linewidth=.5, zorder=0)
 
        # 20210823 Change xlabels to ICD-10 names
        xlabels = icd10_mapping

        # add GWAS significant line
        if gwas_sign_line is True:
            ax.axhline(y=np.log10(gwasp), linestyle='--', color='#7d7d7d', linewidth=1)
            ax.axhline(y=-np.log10(gwasp), linestyle='--', color='#7d7d7d', linewidth=1)
        if markernames is not None:
            marker.geneplot_mhat_logp(df, 
                                      markeridcol, 
                                      chromo, 
                                      'tpval', 
                                      gwasp, 
                                      markernames, 
                                      gfont, 
                                      gstyle, 
                                      ax=ax, 
                                      plotlabelrotation=plotlabelrotation)

        ax.margins(x=0)
        ax.margins(y=0)
        ax.set_xticks(xticks)
        ax.set_yticks(np.arange(mainmin, (mainmax+10), 10)) # set yticks for bottom subplot here
        ax.tick_params(axis='y', labelsize=axlabelfontsize-8)
        ax.set_ylim([mainmin, mainmax]) # limit for bottom subplot here
        ax.set_xticklabels(xlabels, fontsize=axtickfontsize, rotation=ar)
        ax.spines['top'].set_visible(False)
        
        ax2.margins(x=0)
        ax2.margins(y=0)
        ax2.set_yticks(np.arange(topmin, (topmax+10), 10)) # set yticks for top subplot here
        ax2.tick_params(axis='y', labelsize=axlabelfontsize-8)
        ax2.set_ylim([topmin, topmax]) # limit for top subplot here
        ax2.spines['bottom'].set_visible(False)
        ax2.xaxis.set_visible(False)
        kwargs=dict(transform=ax2.transAxes, linespacing=1, fontweight=650)
        ax2.set_title(figtitle, fontsize=axlabelfontsize+2, **kwargs)
        
        (ymin, ymax) = (0, topmax+10)
        ax0.set_ylim([ymin, ymax])
        
        ax0.text(0.5,ymax/2,_y, fontsize=axlabelfontsize, fontname=axlabelfontname, rotation = 90, va='center')
        
        #if ylm:
          # ylm = np.arange(ylm[0], ylm[1], ylm[2])
        #else:
           # ylm = np.concatenate((np.arange(0,min(df['tpval2']-10),-yskip), np.arange(0, max(df['tpval']+10), yskip)))
            #ax.set_yticks(ylm)
        #ax.set_yticklabels(ylm.astype(int), fontsize=axtickfontsize, fontname=axtickfontname);
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        ax.set_xlabel(_x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        ax.get_yaxis().get_label().set_visible(False)
        #ax.set_ylabel(_y, fontsize=axlabelfontsize, fontname=axlabelfontname)
        
        #ax0.text(0.5,0.5,_y, fontsize=axlabelfontsize, fontname=axlabelfontname, rotation=90, va='center')
       
        # 20220301 Invisible ticks for Fig 1; also changing axes colors
        # https://www.delftstack.com/howto/matplotlib/how-to-hide-axis-text-ticks-and-or-tick-labels-in-matplotlib/
        if invisible_ticks:
            ax.yaxis.set_visible(False)
            ax.xaxis.set_visible(False)
            ax2.yaxis.set_visible(False)
            ax2.xaxis.set_visible(False)
            
            # 20220411 change spine color to specific identified race and ethnicity:
            # https://stackoverflow.com/questions/1982770/matplotlib-changing-the-color-of-an-axis
            ax.spines['bottom'].set_color(axescolors)
            ax.spines['top'].set_color(axescolors) 
            ax.spines['right'].set_color(axescolors)
            ax.spines['left'].set_color(axescolors)
        
            ax2.spines['bottom'].set_color(axescolors)
            ax2.spines['top'].set_color(axescolors) 
            ax2.spines['right'].set_color(axescolors)
            ax2.spines['left'].set_color(axescolors)
            
            # 20220411 set line width to be thicker:
            # https://stackoverflow.com/questions/2553521/setting-axes-linewidth-without-changing-the-rcparams-global-dict
            # change all spines
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(8)
                ax2.spines[axis].set_linewidth(8)
                
            # 20220411 set line width to be thicker for the cutoffs
            kwargs = dict(transform=ax.transAxes, color=axescolors, clip_on=False, linewidth=8)
            ax.plot((-d, +d), (1, 1), **kwargs)  # bottom-left diagonal
            ax.plot((1 - d, 1 + d), (1, 1), **kwargs)  # bottom-right diagonal
           
            # Add line breaks (20220112)
            # Plot the same data above on both additional axes
            
            # ax2 is top subplot
            kwargs.update(transform=ax2.transAxes)
            ax2.plot((-d, +d), (0, 0), **kwargs)        # top-left diagonal
            ax2.plot((1 - d, 1 + d), (0, 0), **kwargs)  # top-right diagonal
            
        #plt.suptitle(figtitle, fontsize=axlabelfontsize+2)
        general.get_figure(show, r, figtype, figname)
        return fig, ax
    
    
    def miami(df="dataframe", chromo=None, logp1=None, logp2=None, color=None, dim=(10,10), r=300, ar=90, gwas_sign_line=False,
             gwasp=5E-08, dotsize=8, markeridcol=None, markernames=None, gfont=8, valpha=1, show=False, figtype='png',
             axxlabel=None, axylabel=None, axlabelfontsize=9, axlabelfontname="Arial", axtickfontsize=9, figtitle = 'miami plot',
             label1='firstgroup', label2 = 'secondgroup',
             axtickfontname="Arial", ylm=None, gstyle=1, yskip = 1, plotlabelrotation = 0, figname='miami', invert = False, fig = None, ax = None, icd10_mapping=codemap3):
        
        # 20210823 Added icd10_mapping parameter to add in the icd-10 names

        _x, _y = 'Chromosomes', r'$ -log_{10}(P)$'

        df['tpval'] = df[logp1]
        df['tpval2'] = -df[logp2]
        df = df.sort_values(chromo)

        df['ind'] = range(len(df))
        df_group = df.groupby(chromo)

        rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                               '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                               '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                               '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')
        #color_list = sample(rand_colors, df[chromo].nunique())
        color_list = rand_colors[:df[chromo].nunique()]

        xlabels = []
        xticks = []
        
        if fig is None:
            #fig, ax = plt.subplots(figsize = dim)
            fig, (ax0, ax) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1,20]}, figsize = dim)
            ax0.axis('off')
            fig.tight_layout()

        i=0
        for label, df1 in df.groupby(chromo):
            df1.plot(kind='scatter', x='ind', y='tpval', color=color_list[i], s=dotsize, alpha=valpha, ax=ax)
            df1.plot(kind='scatter', x='ind', y='tpval2', color=color_list[i], s=dotsize, alpha=valpha, ax=ax)
            df1_max_ind = df1['ind'].iloc[-1]
            df1_min_ind = df1['ind'].iloc[0]
            xlabels.append(label)
            xticks.append((df1_max_ind - (df1_max_ind - df1_min_ind) / 2))
            i += 1

        ax.axhline(y=0, color='#7d7d7d', linewidth=.5, zorder = 0)
        
        # 20210823 Change xlabels to ICD-10 names
        xlabels = codemap3

        # add GWAS significant line
        if gwas_sign_line is True:
            ax.axhline(y=np.log10(gwasp), linestyle='--', color='#7d7d7d', linewidth=1)
            ax.axhline(y=-np.log10(gwasp), linestyle='--', color='#7d7d7d', linewidth=1)
        if markernames is not None:
            marker.geneplot_mhat_logp(df, markeridcol, chromo, 'tpval', gwasp, markernames, gfont, gstyle, ax=ax, plotlabelrotation=plotlabelrotation)
            marker.geneplot_mhat_logp(df, markeridcol, chromo, 'tpval2', gwasp, markernames, gfont, gstyle, ax=ax, plotlabelrotation=-plotlabelrotation, vertalign='top')

        ax.margins(x=0)
        ax.margins(y=0)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, fontsize=axlabelfontsize) # added 20210818
        (ymin, ymax) = (min(df['tpval2']-1)-10, max(df['tpval']+1)+10)
        ax.set_ylim([ymin, ymax])
        ax0.set_ylim([ymin, ymax])
        
        ax0.text(0,ymin/2,label2, fontsize=axlabelfontsize, fontname=axlabelfontname, rotation=90, va='center')
        ax0.text(0,ymax/2,label1, fontsize=axlabelfontsize, fontname=axlabelfontname, rotation=90, va='center')
        
        if ylm:
            ylm = np.arange(ylm[0], ylm[1], ylm[2])
        else:
            ylm = np.concatenate((np.arange(0,min(df['tpval2']-10),-yskip), np.arange(0, max(df['tpval']+10), yskip)))
            ax.set_yticks(ylm)
        # Test whether setting fontsize here helps with making ICD10 categories bigger (20210818)
        ax.set_xticklabels(xlabels, fontsize=axtickfontsize, rotation=ar)
        ax.set_yticklabels(ylm.astype(int), fontsize=axtickfontsize, fontname=axtickfontname);
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        ax.set_xlabel(_x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        ax.get_yaxis().get_label().set_visible(False)
        
        ax0.text(.5,0,_y, fontsize=axlabelfontsize, fontname=axlabelfontname, rotation = 90, va = 'center')
        
        plt.title(figtitle, fontsize=axlabelfontsize+2)
        general.get_figure(show, r, figtype, figname)
        return fig, ax