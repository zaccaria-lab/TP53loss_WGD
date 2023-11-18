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
import pingouin as pg
import itertools

from itertools import cycle
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
from itertools import combinations

from statannot import add_stat_annotation

from matplotlib.colors import LinearSegmentedColormap

from gtfparse import read_gtf

plt.style.use('ggplot')
sns.set_style("whitegrid")
plt.rcParams["axes.grid"] = True
plt.rcParams["axes.edgecolor"] = "k"
plt.rcParams["axes.linewidth"]  = 1.5



def prep_data(DATA_DIR='../../TP53loss_WGD_data/'):
    plates = pd.read_excel(os.path.join(DATA_DIR, 'SH_plates.xlsx'))
    plates['Plate number'] = plates['Plate number'].apply(lambda v : v.replace('I_D-H', 'IDH'))
    plates['PEAK'] = plates[['Plate ID', 'Plate number', 'Library numbers']].agg('_'.join, axis=1)

    descr = pd.read_csv(os.path.join(DATA_DIR, 'cell.descriptions.csv'))
    descr = descr.dropna()
    descr['cell.group'] = descr['cell.group'].apply(lambda v : v.replace('I_D-H', 'IDH'))
    descr['cell'] = descr['cell'].apply(lambda v : v.replace('I_D-H', 'IDH'))
    descr['PEAK'] = descr.apply(lambda r : '{}_{}_{}-{}'.format(*str(r['cell.group']).replace('-', '_').split('_')), axis=1)
    descr['status'] = descr['tx.status'] + '_' + descr['p53.status']

    good = pd.read_csv(os.path.join(DATA_DIR, 'good_quality_tumour_cells.tsv'))
    good['cell'] = good['tumour_cells'].apply(lambda v : v.replace('.bam', ''))

    cells = set(descr[descr['tx.status'].notna()]['cell'])
    check_uniq = (lambda L : L[0] if len(set(L)) >= 1 else len(set(L)))
    join = os.path.join
    walk = (lambda : os.walk(os.path.join(DATA_DIR, 'scpeaks/')))
    checkfile = (lambda c, n : '{}.bam'.format(c.replace('IDH', 'I_D-H')) in n)
    finddata = (lambda c : check_uniq([join(root, name) for root, _, files in walk() for name in files if checkfile(c, name)]))
    data = {c : finddata(c) for c in cells}
    data = {c : pd.read_csv(data[c], sep='\t') for c in data if isinstance(data[c], str)}
    cells = set(data.keys())

    form = (lambda p : (p[0], p[1], p[2].split('-')))
    isnum = (lambda v : v.isnumeric())
    compare = (lambda i, x, r, p : i == p[0] and x == p[1] and isnum(r[0]) and isnum(r[1]) and int(r[0]) <= int(p[2]) <= int(r[1]))
    getcells = (lambda i, x, r : filter(lambda e : compare(i, x, r, e.split('_')), cells))
    peaks = {p : list(getcells(*form(p.split('_')))) for p in set(plates['PEAK']).union(set(descr['PEAK']))}
    peaks = {p : peaks[p] for p in peaks if len(peaks[p]) > 0}
    mappk = {e : p for p in peaks for e in peaks[p]}
    ref = data[list(data.keys())[0]][['seqnames', 'start', 'end']]
    widths = data[list(data.keys())[0]]['width']
    tot = widths.sum()
    diff = (lambda i, j : np.sum(widths * (data[i]['copy.number'] != data[j]['copy.number'])) / tot if mappk[i] == mappk[j] else 0.0)
    flip = (lambda i, j, d : {(i, j, d), (j, i, d)})
    form = (lambda e : '{}_{}'.format(mappk[e],e))
    dfdiff = pd.DataFrame([{'Cell1' : form(p[0]), 'Cell2' : form(p[1]), 'Diff' : p[2]} for i, j in combinations(mappk, 2) for p in flip(i, j, diff(i, j))])
    tadiff = pd.pivot_table(dfdiff, index='Cell1', columns=['Cell2'], values='Diff')
    np.fill_diagonal(tadiff.values, 0)

    status = descr[['cell', 'status', 'numeric.ploidy']]
    status = status[status['cell'].isin(set(mappk.keys()))]
    status['WGD.status'] = status['numeric.ploidy'] >= 3
    status['isnormal'] = status.apply(lambda r : (np.sum(widths * (data[r['cell']]['copy.number'] == 2)) / tot) > 0.99, axis=1)
    status['Peak'] = status['cell'].apply(lambda c : mappk[c])
    tumour = status[status['isnormal']==False]
    tumour = tumour.sort_values('numeric.ploidy')
    ploidy = {r['PEAK'] : int(round(float(r['Sample description'].split('X')[0].split('-')[0].replace('~', '')))) for i, r in plates.iterrows()}

    props = (lambda L : L / np.sum(L))
    shaneve = (lambda P : ((-np.sum(P * np.log(P))) / np.log(P.shape[0])) if len(P) > 1 else 0.0)
    shan = []
    good_tumour = tumour[tumour['cell'].isin(set(good['cell']))].copy()
    for S, gdf in good_tumour.groupby('status'):
        for P, gtodf in gdf.groupby('Peak'):
            if len(gtodf.index) > 1:
                linkage = hier.linkage(gtodf['cell'].apply(lambda e : data[e]['copy.number']), method='single', metric='hamming', optimal_ordering=True)
                clus = hier.fcluster(linkage, t=0.2, criterion='distance')
                shan.append({'Status' : S, 'Peak' : P, 'Shannon eve.' : shaneve(props(np.array(list(Counter(clus).values())))), 'WGD Status' : gtodf['WGD.status'].value_counts().idxmax(), 'Ploidy' : ploidy[P]})
                assert not np.isnan(shan[-1]['Shannon eve.'])
    shan = pd.DataFrame(shan)
    shan['Treatment'] = shan['Status'].str.contains('ERL-Resis')
    shan['TP53 mutated'] = shan['Status'].str.contains('_FL')
    shan['ITH'] = shan['Shannon eve.']
    shan = shan.merge(right=descr[['PEAK', 'mouse', 'Tumour_ID']].drop_duplicates(), right_on='PEAK', left_on='Peak', how='left')
    shan = shan[['Status', 'Peak', 'Shannon eve.', 'Ploidy', 'Treatment', 'TP53 mutated','ITH', 'mouse', 'Tumour_ID']]
    shan = pd.merge(pd.DataFrame([(mappk[cell], cell) for cell in mappk if cell in good['cell'].unique()], columns=['Peak', 'CELL']).drop_duplicates(), shan)
    tploidy = {cell : (data[cell]['copy.number'] * data[cell]['width']).sum() / data[cell]['width'].sum() for cell in data}
    shan['CellPloidy'] = shan['CELL'].map(tploidy)
    shan['WGDPloidy'] = shan['CellPloidy'] >= 3.0
    shan['WGDCN'] = shan['CELL'].map({cell : (data[cell]['copy.number'] > 2).sum() > (data[cell]['copy.number'] <= 2).sum() for cell in data})

    return {
        'plates' : plates,
        'descr' : descr,
        'good' : good,
        'cells' : cells,
        'data' : data,
        'peaks' : peaks,
        'mappk' : mappk,
        'ref' : ref,
        'dfdiff' : dfdiff,
        'tadiff' : tadiff,
        'status' : status,
        'tumour' : tumour,
        'ploidy' : ploidy,
        'widths' : widths,
        'tot' : tot,
        'shan' : shan,
        'tploidy' : tploidy
    }


def necklace_plots(widths, data, tot, mappk, ploidy, tumour, good, descr, **_):
    diff = (lambda i, j : np.sum(widths * (data[i]['copy.number'] != data[j]['copy.number'])) / tot if mappk[i] == mappk[j] else np.nan)
    form = (lambda e : '{}+{}+{}'.format(ploidy[mappk[e]], mappk[e], e))
    rec = (lambda p : {'Cell1' : form(p[0]), 'Cell2' : form(p[1]), 'Diff' : p[2]})
    flip = (lambda i, j, d : {(i, j, d), (j, i, d)})
    pal = sns.color_palette('Reds', len(sorted(set(ploidy.values()))))
    colormap = dict(zip(sorted(set(ploidy.values())), pal))
    tumour_sel = tumour.merge(right=descr[['PEAK', 'mouse']], right_on='PEAK', left_on='Peak', how='left')
    tumour_sel = tumour_sel[(tumour_sel['status']!='ERL-Resis_WT')|(tumour_sel['mouse'].isin({'585998', '591564', '594254', '599368', '599370', '582356', '584260', '585996', '584254'}))]
    tumour_sel[tumour_sel['status']=='ERL-Resis_WT']['mouse'].unique()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for S, gdf in tumour_sel.groupby('status'):
            dfdiff = pd.DataFrame([rec(p) for i, j in combinations(set(gdf['cell']).intersection(set(good['cell'])), 2) for p in flip(i, j, diff(i, j))])
            tadiff = pd.pivot_table(dfdiff[['Cell1', 'Cell2', 'Diff']], index='Cell1', columns=['Cell2'], values='Diff')
            print('{} : {}'.format(S, len(set(dfdiff['Cell1']).union(set(dfdiff['Cell2'])))))
            np.fill_diagonal(tadiff.values, 0)
            para = {}
            para['data'] = tadiff
            para['yticklabels'] = False
            para['row_cluster'] = True
            para['xticklabels'] = False
            para['col_cluster'] = True
            para['rasterized'] = True
            para['metric'] = 'jaccard'
            para['cmap'] = 'YlGnBu'
            para['square'] = 'True'
            para['row_colors'] = pd.DataFrame({'Cell1' : tadiff.index, 'Peak' : pd.Series(tadiff.index).apply(lambda v : colormap[ploidy[v.split('+')[1]]])}).set_index('Cell1')
            g = sns.clustermap(**para)
            g.ax_row_dendrogram.set_xlim([0,0])
            g.ax_col_dendrogram.set_xlim([0,0])
            g.ax_row_dendrogram.set_ylim([0,0])
            g.ax_col_dendrogram.set_ylim([0,0])
            g.fig.suptitle(S)
            plt.savefig('{}.png'.format(S), bbox_inches='tight', dpi=600)
            dfdiff.to_csv('cell-differences-{}.tsv.gz'.format(S), sep='\t', index=False)


def cell_ploidy(shan, **_):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ddf = shan.groupby('Tumour_ID')[['Treatment', 'TP53 mutated', 'CellPloidy', 'WGDPloidy']].agg({'Treatment' : 'first', 
                                                                                                    'TP53 mutated' : 'first', 
                                                                                                    'CellPloidy' : 'mean',
                                                                                                    'WGDPloidy' : (lambda v : v.value_counts().index[0])})
        ddf['Status'] = ddf.apply(lambda r : '{} {}'.format('E' if not r['TP53 mutated'] else 'EP', 'Treated' if r['Treatment'] else 'Untreated'), axis=1)
        plt.figure(figsize=(5, 3))
        g = sns.swarmplot(data=ddf, x='Status', y='CellPloidy', order=['E Untreated', 'E Treated', 'EP Untreated', 'EP Treated'], s=9, palette=['#e2c18d', '#68c9c3', '#c3a55f', '#218f85'])
        sns.despine(ax=g, top=True, right=True)
        plt.ylim(1, 6)
        plt.plot((-1, 5), (3, 3), '--k', linewidth=2)
        add_stat_annotation(g, data=ddf, x='Status', y='CellPloidy', order=['E Untreated', 'E Treated', 'EP Untreated', 'EP Treated'],
                            box_pairs=[('E Untreated', 'EP Untreated'), ('E Treated', 'EP Untreated'), ('E Treated', 'EP Treated')],
                            test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
        ddf.to_csv('cell-ploidies.tsv.gz', sep='\t', index=False)


def shannon_eveness(shan, **_):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ddf = shan.groupby('Tumour_ID')[['Treatment', 'TP53 mutated', 'Shannon eve.']].agg({'Treatment' : 'first', 
                                                                                            'TP53 mutated' : 'first', 
                                                                                            'Shannon eve.' : 'mean'})
        ddf['Status'] = ddf.apply(lambda r : '{} {}'.format('E' if not r['TP53 mutated'] else 'EP', 'Treated' if r['Treatment'] else 'Untreated'), axis=1)
        plt.figure(figsize=(5, 3))
        g = sns.swarmplot(data=ddf, x='Status', y='Shannon eve.', order=['E Untreated', 'E Treated', 'EP Untreated', 'EP Treated'], s=9, palette=['#e2c18d', '#68c9c3', '#c3a55f', '#218f85'])
        sns.despine(ax=g, top=True, right=True)
        add_stat_annotation(g, data=ddf, x='Status', y='Shannon eve.', order=['E Untreated', 'E Treated', 'EP Untreated', 'EP Treated'],
                            box_pairs=[('E Treated', 'EP Untreated'), ('E Untreated', 'EP Untreated'), ('E Treated', 'EP Treated'), ('E Untreated', 'EP Treated')],
                            test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
        ddf.to_csv('shannon-evenness.tsv.gz', sep='\t', index=False)


def cn_heatmap_Euntreated(shan, data, ref, **_):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tt = shan.sort_values(['TP53 mutated', 'Treatment', 'mouse', 'Tumour_ID', 'CellPloidy'])
        tt = tt[(tt['TP53 mutated']==False)&(tt['Treatment']==False)]
        matrix = tt['CELL'].apply(lambda e : data[e]['copy.number'])
        row_palette = cycle(['#8dd3c7', '#ffffb3', '#bebada'])
        row_colors = {c : next(row_palette) for c in tt['Tumour_ID'].unique()}
        chr_palette = cycle(['#525252', '#969696', '#cccccc'])
        chr_colors = {c : next(chr_palette) for c in sorted(map(lambda v : int(v) if v != 'X' else 23, ref['seqnames'].unique()))}
        para = {}
        para['data'] = matrix
        para['yticklabels'] = False
        para['row_cluster'] = False
        para['xticklabels'] = False
        para['col_cluster'] = False
        para['rasterized'] = True
        para['method'] = 'single'
        para['metric'] = 'hamming'
        para['cmap'] = 'RdBu_r'
        para['center'] = 2
        para['vmin'] = 0
        para['vmax'] = 10
        para['figsize'] = (12, 10)
        para['row_colors'] = [row_colors[c] for c in tt['Tumour_ID']]
        para['col_colors'] = list(map(lambda v : chr_colors[int(v) if v != 'X' else 23], ref['seqnames']))
        g = sns.clustermap(**para)
        g.ax_row_dendrogram.set_xlim([0,0])
        g.ax_col_dendrogram.set_xlim([0,0])
        g.ax_row_dendrogram.set_ylim([0,0])
        g.ax_col_dendrogram.set_ylim([0,0])
        pd.DataFrame(matrix).rename(index=tt['CELL']).rename(columns=data[tt['CELL'].values[0]]['seqnames'].astype(str) + ':' + data[tt['CELL'].values[0]]['start'].astype(str) 
                                                                + '-' + data[tt['CELL'].values[0]]['end'].astype(str)).to_csv('heatmap_Euntreated.tsv', sep='\t')
        plt.savefig('heatmap_Euntreated.pdf', bbox_inches='tight')


def cn_heatmap_EPuntreated(shan, data, ref, **_):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tt = shan.sort_values(['TP53 mutated', 'Treatment', 'mouse', 'Tumour_ID', 'CellPloidy'])
        tt = tt[(tt['TP53 mutated']==True)&(tt['Treatment']==False)]
        matrix = tt['CELL'].apply(lambda e : data[e]['copy.number'])
        row_palette = cycle(['#8dd3c7', '#ffffb3', '#bebada'])
        row_colors = {c : next(row_palette) for c in tt['Tumour_ID'].unique()}
        chr_palette = cycle(['#525252', '#969696', '#cccccc'])
        chr_colors = {c : next(chr_palette) for c in sorted(map(lambda v : int(v) if v != 'X' else 23, ref['seqnames'].unique()))}
        para = {}
        para['data'] = matrix
        para['yticklabels'] = False
        para['row_cluster'] = False
        para['xticklabels'] = False
        para['col_cluster'] = False
        para['rasterized'] = True
        para['method'] = 'single'
        para['metric'] = 'hamming'
        para['cmap'] = 'RdBu_r'
        para['center'] = 2
        para['vmin'] = 0
        para['vmax'] = 10
        para['figsize'] = (12, 10)
        para['row_colors'] = [row_colors[c] for c in tt['Tumour_ID']]
        para['col_colors'] = list(map(lambda v : chr_colors[int(v) if v != 'X' else 23], ref['seqnames']))
        g = sns.clustermap(**para)
        g.ax_row_dendrogram.set_xlim([0,0])
        g.ax_col_dendrogram.set_xlim([0,0])
        g.ax_row_dendrogram.set_ylim([0,0])
        g.ax_col_dendrogram.set_ylim([0,0])
        pd.DataFrame(matrix).rename(index=tt['CELL']).rename(columns=data[tt['CELL'].values[0]]['seqnames'].astype(str) + ':' + data[tt['CELL'].values[0]]['start'].astype(str) 
                                                                + '-' + data[tt['CELL'].values[0]]['end'].astype(str)).to_csv('heatmap_EPuntreated.tsv', sep='\t')
        plt.savefig('heatmap_EPuntreated.pdf', bbox_inches='tight')


def cn_heatmap_Etreated(shan, data, ref, **_):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tt = shan.sort_values(['TP53 mutated', 'Treatment', 'mouse', 'Tumour_ID', 'CellPloidy'])
        tt = tt[(tt['TP53 mutated']==False)&(tt['Treatment']==True)]
        matrix = tt['CELL'].apply(lambda e : data[e]['copy.number'])
        row_palette = cycle(['#8dd3c7', '#ffffb3', '#bebada'])
        row_colors = {c : next(row_palette) for c in tt['Tumour_ID'].unique()}
        chr_palette = cycle(['#525252', '#969696', '#cccccc'])
        chr_colors = {c : next(chr_palette) for c in sorted(map(lambda v : int(v) if v != 'X' else 23, ref['seqnames'].unique()))}
        para = {}
        para['data'] = matrix
        para['yticklabels'] = False
        para['row_cluster'] = False
        para['xticklabels'] = False
        para['col_cluster'] = False
        para['rasterized'] = True
        para['method'] = 'single'
        para['metric'] = 'hamming'
        para['cmap'] = 'RdBu_r'
        para['center'] = 2
        para['vmin'] = 0
        para['vmax'] = 10
        para['figsize'] = (12, 10)
        para['row_colors'] = [row_colors[c] for c in tt['Tumour_ID']]
        para['col_colors'] = list(map(lambda v : chr_colors[int(v) if v != 'X' else 23], ref['seqnames']))
        g = sns.clustermap(**para)
        g.ax_row_dendrogram.set_xlim([0,0])
        g.ax_col_dendrogram.set_xlim([0,0])
        g.ax_row_dendrogram.set_ylim([0,0])
        g.ax_col_dendrogram.set_ylim([0,0])
        pd.DataFrame(matrix).rename(index=tt['CELL']).rename(columns=data[tt['CELL'].values[0]]['seqnames'].astype(str) + ':' + data[tt['CELL'].values[0]]['start'].astype(str) 
                                                                + '-' + data[tt['CELL'].values[0]]['end'].astype(str)).to_csv('heatmap_Etreated.tsv', sep='\t')
        plt.savefig('heatmap_Etreated.pdf', bbox_inches='tight')


def cn_heatmap_EPtreated(shan, data, ref, **_):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tt = shan.sort_values(['TP53 mutated', 'Treatment', 'mouse', 'Tumour_ID', 'CellPloidy'])
        tt = tt[(tt['TP53 mutated']==True)&(tt['Treatment']==True)]
        matrix = tt['CELL'].apply(lambda e : data[e]['copy.number'])
        row_palette = cycle(['#8dd3c7', '#ffffb3', '#bebada'])
        row_colors = {c : next(row_palette) for c in tt['Tumour_ID'].unique()}
        chr_palette = cycle(['#525252', '#969696', '#cccccc'])
        chr_colors = {c : next(chr_palette) for c in sorted(map(lambda v : int(v) if v != 'X' else 23, ref['seqnames'].unique()))}
        para = {}
        para['data'] = matrix
        para['yticklabels'] = False
        para['row_cluster'] = False
        para['xticklabels'] = False
        para['col_cluster'] = False
        para['rasterized'] = True
        para['method'] = 'single'
        para['metric'] = 'hamming'
        para['cmap'] = 'RdBu_r'
        para['center'] = 2
        para['vmin'] = 0
        para['vmax'] = 10
        para['figsize'] = (12, 10)
        para['row_colors'] = [row_colors[c] for c in tt['Tumour_ID']]
        para['col_colors'] = list(map(lambda v : chr_colors[int(v) if v != 'X' else 23], ref['seqnames']))
        g = sns.clustermap(**para)
        g.ax_row_dendrogram.set_xlim([0,0])
        g.ax_col_dendrogram.set_xlim([0,0])
        g.ax_row_dendrogram.set_ylim([0,0])
        g.ax_col_dendrogram.set_ylim([0,0])
        pd.DataFrame(matrix).rename(index=tt['CELL']).rename(columns=data[tt['CELL'].values[0]]['seqnames'].astype(str) + ':' + data[tt['CELL'].values[0]]['start'].astype(str) 
                                                                + '-' + data[tt['CELL'].values[0]]['end'].astype(str)).to_csv('heatmap_EPtreated.tsv', sep='\t')
        plt.savefig('heatmap_EPtreated.pdf', bbox_inches='tight')
