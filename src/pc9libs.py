import os, sys, glob
import warnings
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import scipy.cluster
import scipy.cluster.hierarchy as hier

from itertools import cycle
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
from itertools import combinations
from functools import reduce

import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

plt.style.use('ggplot')
sns.set_style("whitegrid")
plt.rcParams["axes.grid"] = True
plt.rcParams["axes.edgecolor"] = "k"
plt.rcParams["axes.linewidth"]  = 1.5


def prep_data(DATA_DIR='../../TP53loss_WGD_data/'):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drivers = pd.read_csv(os.path.join(DATA_DIR, 'drivermutations.txt'), sep='\t')
        drivers['CellLine'] = drivers['Var2']
        ploidy = pd.read_csv(os.path.join(DATA_DIR, 'PC9Ploidy.txt'), sep='\t')
        ploidy['CellLine'] = ploidy['bam'].str.replace('S_PC9|_SU', '')
        ploidy['Classed_ploidy'] = ploidy['Ploidy'].apply(lambda p : 'TRI' if abs(3 - p) < abs(6 - p) else 'HEX')
        ploidy['Parental'] = ploidy['bam'].str.contains('PC9A')
        tsg = pd.read_csv(os.path.join(DATA_DIR, 'TSG.absoluteCN.txt'), sep='\t')
        tsg['CellLine'] = tsg['Var2'].str.replace('S_PC9|_SU', '')
        tsg['Deletion'] = 0 - tsg['value']
        og = pd.read_csv(os.path.join(DATA_DIR, 'OG.absoluteCN.txt'), sep='\t')
        og['CellLine'] = og['Var2'].str.replace('S_PC9|_SU', '')
        og['Amplification'] = og['value']
        mergedf = (lambda x, y : pd.merge(left=x, right=y, how='left', on='CellLine'))
        comb = reduce(lambda x, y : mergedf(x, y), [ploidy, drivers, tsg, og])
        comb = comb.sort_values(['Classed_ploidy', 'Parental', 'T790M', 'CellLine', 'Amplification'], ascending=[False, True, True, True, True])
        comb = comb[~comb['CellLine'].str.contains('_O0|parental')]
        highamp = pd.read_csv(os.path.join(DATA_DIR, 'PC9relativeCN.txt'), sep='\t').unstack().reset_index()
        highamp.columns = ['SAMPLE', 'GENE', 'COPIES']
        highamp['CELLLINE'] = highamp['SAMPLE'].str.replace('S_PC9|_SU', '')
        highamp['PLOIDY'] = highamp['SAMPLE'].map(ploidy[['bam', 'Classed_ploidy']].set_index('bam')['Classed_ploidy'])
        resgain = {'MET', 'MAPK1', 'PIK3CA', 'CCND1', 'CCND3', 'CDK6', 'CDK4', 'MAP2K4', 'MAP3K1', 'STAT6', 'BRAF', 'ERBB2', 'NRAS', 'KRAS', 'FGFR1', 'FGFR3', 'ERBB2', 'EGFR', 'FGFR2', 'FGFR4', 'BRAF', 'NF1', 'HRAS', 'KRAS'}
        return {
            'ploidy' : ploidy, 
            'drivers' : drivers,
            'highamp' : highamp,
            'resgain' : resgain,
            'comb' : comb
        }


def mut_table(ploidy, drivers, comb, **_):
    data = pd.merge(left=ploidy, right=drivers, how='left', on='CellLine')
    data = data[~data['CellLine'].str.contains('_O0|parental')]
    data['MUT'] = (~data['value'].isna()).astype(int)
    data['Ploidy'] = data['Classed_ploidy'].apply(lambda v : '#008faa' if v=='TRI' else '#e83c5b')
    data['Parental ploidy'] = data['Parental'].apply(lambda v : '#008faa' if v==False else '#e83c5b')
    data.to_csv('PC9-mutations.tsv.gz', sep='\t', index=False)
    table = pd.pivot_table(data, index='Var1', columns='CellLine', values='MUT', aggfunc='first')
    table = table.reindex(comb['CellLine'].unique(), axis=1)
    g = sns.clustermap(table,
                    figsize=(13, 5),
                    cmap='Greys',
                    xticklabels=False,
                    yticklabels=True,
                    row_cluster=False,
                    col_cluster=False,
                    col_colors=data[['CellLine', 'Ploidy', 'Parental ploidy']].drop_duplicates().set_index('CellLine'),
                    #colors_ratio=0.03,
                    dendrogram_ratio=0.0001)
    g.cax.set_visible(False)
    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')
    g.ax_heatmap.yaxis.set_ticks_position('left')
    g.ax_col_colors.yaxis.set_ticks_position('left')
    for _, spine in g.ax_heatmap.spines.items():
        spine.set_visible(True)
    g.ax_heatmap.vlines(np.arange(*g.ax_heatmap.get_xlim()), *g.ax_heatmap.get_xlim(), colors='grey', linewidth=1, linestyles='--')
    plt.savefig('PC9drivers.pdf', bbox_inches='tight')


def all_amp(highamp, resgain, comb, **_):
    data = highamp[highamp['GENE'].isin(resgain)]
    data = data[~data['CELLLINE'].str.contains('_O0|parental')]
    data = data[data['COPIES']>=1.3]
    table = pd.pivot_table(data, index='GENE', columns='CELLLINE', values='COPIES', aggfunc='first')
    table = table.reindex(comb['CellLine'].unique(), axis=1)
    g = sns.clustermap(table,
                    figsize=(14, 5),
                    cmap='Reds',
                    vmax=3,
                    xticklabels=True,
                    yticklabels=True,
                    row_cluster=False,
                    col_cluster=False,
                    dendrogram_ratio=0.0001)
    g.cax.set_visible(False)
    for _, spine in g.ax_heatmap.spines.items():
        spine.set_visible(True)
    plt.savefig('PC9ampresist.pdf', bbox_inches='tight')


def selected_amp(ploidy, drivers, highamp, resgain, comb, **_):
    parsel = pd.merge(left=ploidy, right=drivers, how='left', on='CellLine')
    parsel = parsel[~parsel['CellLine'].str.contains('_O0|parental')]
    parsel = parsel[((parsel['Classed_ploidy']=='TRI')&(parsel['Parental']==True))|((parsel['Classed_ploidy']=='HEX')&(parsel['Parental']==False))]
    parsel['Ploidy'] = parsel['Classed_ploidy'].apply(lambda v : '#008faa' if v=='TRI' else '#e83c5b')
    parsel['Parental ploidy'] = parsel['Parental'].apply(lambda v : '#008faa' if v==False else '#e83c5b')
    data = highamp[highamp['GENE'].isin(resgain)]
    data = data[~data['CELLLINE'].str.contains('_O0|parental')]
    data = data[data['CELLLINE'].isin(parsel['CellLine'].unique())]
    data = data[data['COPIES']>=1]
    table = pd.pivot_table(data, index='GENE', columns='CELLLINE', values='COPIES', aggfunc='first')
    table = table.reindex([x for x in comb['CellLine'].unique() if x in parsel['CellLine'].unique()], axis=1)
    g = sns.clustermap(table,
                    figsize=(4, 5),
                    cmap='Reds',
                    vmax=3,
                    xticklabels=False,
                    yticklabels=True,
                    row_cluster=False,
                    col_cluster=False,
                    col_colors=parsel.rename(columns={'CellLine' : 'CELLLINE'})[['CELLLINE', 'Ploidy', 'Parental ploidy']].drop_duplicates().set_index('CELLLINE'),
                    colors_ratio=0.03,
                    dendrogram_ratio=0.0001,
                    cbar_kws={"orientation": "horizontal"})
    g.cax.set_visible(False)
    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')
    g.cax.set_position([1, 1, 0.5, 0.03])
    g.ax_heatmap.yaxis.set_ticks_position('left')
    g.ax_col_colors.yaxis.set_ticks_position('left')
    for _, spine in g.ax_heatmap.spines.items():
        spine.set_visible(True)
    plt.savefig('PC9selectedampresist.pdf', bbox_inches='tight')

