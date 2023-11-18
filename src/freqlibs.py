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

from matplotlib.colors import LinearSegmentedColormap

from gtfparse import read_gtf

plt.style.use('ggplot')
sns.set_style("whitegrid")
plt.rcParams["axes.grid"] = True
plt.rcParams["axes.edgecolor"] = "k"
plt.rcParams["axes.linewidth"]  = 1.5


def prep_data(DATA_DIR='../../TP53loss_WGD_data/'):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        humangenes = read_gtf(os.path.join(DATA_DIR, "gencode.v38.annotation.gtf.gz"))
        humangenes = humangenes[humangenes['feature']=='gene']
        humangenes['GENE_NAME'] = humangenes['gene_name'].str.upper()
        humangenes = humangenes[~humangenes['seqname'].isna()]
        humangenes['GENOME'] = tuple(zip(humangenes['seqname'], humangenes['start'], humangenes['end']))
        humangenes = humangenes[['GENE_NAME', 'seqname', 'start', 'end']]

        pc9 = os.path.join(DATA_DIR, 'PC9/*readdepth.txt.gz')
        pc9 = {os.path.basename(f).replace('readdepth.txt.gz', '') : f for f in glob.glob(pc9)}
        parents = {p : [f for f in pc9.keys() if p == '_'.join(f.split('_')[:2])] for p in {f.split('_SU_')[0] for f in pc9}}
        assert set(pc9.keys()) == set(f for p in parents for f in parents[p])
        mapp = {f : p for p in parents for f in parents[p]}
        info = pd.read_excel(os.path.join(DATA_DIR, 'PC9_Corrected_ploidies_4.11.2020.xlsx'))
        info = info.loc[:, ~info.columns.str.contains('^Unnamed')].dropna().copy()
        info['name'] = info['bam'].str.replace('.bam', '')
        extract = (lambda D : D[['name', 'Corrected Ploidy']].set_index('name').to_dict()['Corrected Ploidy'])
        pploidy = info[info['name'].str.contains('SU_O0')].copy()
        pploidy['name'] = pploidy['name'].str.replace('_SU_O0', '')
        pploidy = extract(pploidy)
        assert set(pploidy.keys()) == set(parents.keys())
        sel = (lambda D, p : (D != p) & (D.str.contains('_'.join(p.split('_')[:3]))))
        dploidy = {p : extract(info[sel(info['name'], '{}_SU_O0'.format(p))]) for p in pploidy}
        assert set(f for p in dploidy for f in dploidy[p]) == set(f for p in parents for f in parents[p])
        getp = (lambda p : {'PARENT_LINE' : mapp[p], 'CELLLINE' : p, 'PARENTAL_PLOIDY' : pploidy[mapp[p]], 'PLOIDY' : dploidy[mapp[p]][p]})
        cols = ['CHR', 'START', 'STOP', 'read.pvalue', 'gain.vs.loss', 'final.status', 'read.pvalue.group']
        data = pd.concat((pd.read_csv(pc9[p], sep='\t')[cols].assign(**getp(p)) for p in pc9)).drop_duplicates()

        def getsegs(cellline):
            view = humangenes.merge(right=cellline, left_on='seqname', right_on='CHR', how='left')
            view['OVERLAP'] = np.maximum(view['START'], view['start']) - np.minimum(view['STOP'], view['end'])
            view = view[view['OVERLAP'] > 0].sort_values(['GENE_NAME', 'OVERLAP'], ascending=[True, False])
            return view.groupby('GENE_NAME', sort=False).first().reset_index()
        humancombo = pd.concat((getsegs(cellline) for _, cellline in data.groupby('CELLLINE', sort=False)))

        mousegenes = read_gtf(os.path.join(DATA_DIR, "Mus_musculus.GRCm38.102.gtf.gz"))
        mousegenes = mousegenes[mousegenes['feature']=='gene']
        mousegenes['GENE_NAME'] = mousegenes['gene_name'].str.upper()
        mousegenes = mousegenes[~mousegenes['seqname'].isna()]
        mousegenes['GENOME'] = tuple(zip(mousegenes['seqname'], mousegenes['start'], mousegenes['end']))
        mousegenes = mousegenes[['GENE_NAME', 'seqname', 'start', 'end']]

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

        cells = set(descr[descr['tx.status'].notna()]['cell']).intersection(good['cell'])
        annot = descr[['cell', 'PEAK', 'tx.status', 'p53.status']].set_index('cell').to_dict('index')
        check_uniq = (lambda L : L[0] if len(set(L)) >= 1 else len(set(L)))
        join = os.path.join
        walk = (lambda : os.walk(os.path.join(DATA_DIR, 'scpeaks/')))
        checkfile = (lambda c, n : '{}.bam'.format(c.replace('IDH', 'I_D-H')) in n)
        finddata = (lambda c : check_uniq([join(root, name) for root, _, files in walk() for name in files if checkfile(c, name)]))
        data = {c : finddata(c) for c in cells}
        form = (lambda p : (p[0], p[1], p[2].split('-')))
        isnum = (lambda v : v.isnumeric())
        compare = (lambda i, x, r, p : i == p[0] and x == p[1] and isnum(r[0]) and isnum(r[1]) and int(r[0]) <= int(p[2]) <= int(r[1]))
        getcells = (lambda i, x, r : filter(lambda e : compare(i, x, r, e.split('_')), cells))
        peaks = {p : list(getcells(*form(p.split('_')))) for p in set(plates['PEAK']).union(set(descr['PEAK']))}
        peaks = {p : peaks[p] for p in peaks if len(peaks[p]) > 0}
        mappk = {e : p for p in peaks for e in peaks[p]}
        data = pd.concat((pd.read_csv(data[c], sep='\t').assign(CELL=c, **annot[c]) for c in data if isinstance(data[c], str) 
                                                                                    if c in set(good['cell']) and c in mappk))
        data = data[['seqnames', 'start', 'end', 'width', 'copy.number', 'CELL', 'PEAK', 'tx.status', 'p53.status']]
        data = data.merge(right=descr[['PEAK', 'Tumour_ID']].drop_duplicates(), left_on='PEAK', right_on='PEAK', how='left')
        func = {'copy.number' : 'median', 'tx.status' : 'first', 'p53.status' : 'first'}
        data = data.groupby(['Tumour_ID', 'seqnames', 'start', 'end'])[['copy.number', 'tx.status', 'p53.status']].agg(func).reset_index()
        getploidy = (lambda r : np.sum(r['copy.number'] * (r['end'] - r['start'])) / np.sum(r['end'] - r['start']))
        data['PLOIDY'] = data['Tumour_ID'].map(data.groupby('Tumour_ID')[['start', 'end', 'copy.number']].apply(getploidy))
        data.columns = ['Tumour_ID', 'CHR', 'START', 'STOP', 'copy.number', 'tx.status', 'p53.status', 'PLOIDY']

        def getsegs(peak):
            view = mousegenes.merge(right=peak, left_on='seqname', right_on='CHR', how='left')
            view['OVERLAP'] = np.maximum(view['start'], view['START']) - np.minimum(view['end'], view['STOP'])
            view = view[view['OVERLAP'] > 0].sort_values(['GENE_NAME', 'OVERLAP'], ascending=[True, False])
            return view.groupby('GENE_NAME', sort=False).first().reset_index()
        mousecombo = pd.concat((getsegs(peak) for _, peak in data.groupby('Tumour_ID', sort=False)))

        humancombo['DESC_PLOIDY'] = humancombo['PLOIDY'].apply(lambda v : 3 if abs(v - 3) < abs(v - 6) else 6)
        humancombo['SIGN_AMP'] = (humancombo['read.pvalue'] <= 0.05) & (humancombo['gain.vs.loss'] == 'gain')
        humancombo['SIGN_DEL'] = (humancombo['read.pvalue'] <= 0.05) & (humancombo['gain.vs.loss'] == 'loss')
        count = (lambda C, V : C[V] / C.sum() if V in C else 0.0)
        form = {'FREQ_AMP' : ('SIGN_AMP', (lambda c : count(c.value_counts(), True))), 
                'FREQ_DEL' : ('SIGN_DEL', (lambda c : count(c.value_counts(), True)))}
        humanfreq = humancombo.groupby(['GENE_NAME'])[[form[k][0] for k in form]].agg(**form).reset_index()
        cosmic = pd.read_csv(os.path.join(DATA_DIR, 'COSMIC_19_35_00_2021.tsv'), sep='\t')['Gene Symbol'].str.upper()

        humanratiod = humancombo.groupby('CELLLINE')['SIGN_DEL'].apply(lambda v : v.value_counts()[True] / len(v) if True in v.value_counts() else 0.0).to_numpy()
        humanratiod = sorted(np.sum(scipy.stats.bernoulli.rvs(humanratiod)) / len(humanratiod) for x in range(100000))
        humanratiod = min(humanratiod[-int(len(humanratiod) * 0.05):])
        humanfreq['SEL_DEL'] = humanfreq['FREQ_DEL'] >= humanratiod

        humanratiog = humancombo.groupby('CELLLINE')['SIGN_AMP'].apply(lambda v : v.value_counts()[True] / len(v) if True in v.value_counts() else 0.0).to_numpy()
        humanratiog = sorted(np.sum(scipy.stats.bernoulli.rvs(humanratiog)) / len(humanratiog) for x in range(100000))
        humanratiog = min(humanratiog[-int(len(humanratiog) * 0.05):])
        humanfreq['SEL_AMP'] = humanfreq['FREQ_AMP'] >= humanratiog

        mousecombo['SIGN_AMP'] = mousecombo['copy.number'] > np.round(mousecombo['PLOIDY']).astype(int)
        mousecombo['SIGN_DEL'] = mousecombo['copy.number'] < np.round(mousecombo['PLOIDY']).astype(int)
        count = (lambda C, V : C[V] / C.sum() if V in C else 0.0)
        form = {'FREQ_AMP' : ('SIGN_AMP', (lambda c : count(c.value_counts(), True))), 
                'FREQ_DEL' : ('SIGN_DEL', (lambda c : count(c.value_counts(), True)))}
        mousefreq = mousecombo[(mousecombo['tx.status']=='ERL-Resis')
                            &(mousecombo['p53.status']=='FL')].groupby(['GENE_NAME'])[[form[k][0] for k in form]].agg(**form).reset_index()

        mouseratiod = mousecombo.groupby('Tumour_ID')['SIGN_DEL'].apply(lambda v : v.value_counts()[True] / len(v) if True in v.value_counts() else 0.0).to_numpy()
        mouseratiod = sorted(np.sum(scipy.stats.bernoulli.rvs(mouseratiod)) / len(mouseratiod) for x in range(100000))
        mouseratiod = min(mouseratiod[-int(len(mouseratiod) * 0.05):])
        mousefreq['SEL_DEL'] = mousefreq['FREQ_DEL'] >= mouseratiod

        mouseratiog = mousecombo.groupby('Tumour_ID')['SIGN_AMP'].apply(lambda v : v.value_counts()[True] / len(v) if True in v.value_counts() else 0.0).to_numpy()
        mouseratiog = sorted(np.sum(scipy.stats.bernoulli.rvs(mouseratiog)) / len(mouseratiog) for x in range(100000))
        mouseratiog = min(mouseratiog[-int(len(mouseratiog) * 0.05):])
        mousefreq['SEL_AMP'] = mousefreq['FREQ_AMP'] >= mouseratiog

        humancombo['DESC_PLOIDY'] = humancombo['PLOIDY'].apply(lambda v : 3 if abs(v - 3) < abs(v - 6) else 6)
        humancombo['SIGN_AMP'] = (humancombo['read.pvalue'] <= 0.05) & (humancombo['gain.vs.loss'] == 'gain')
        humancombo['SIGN_DEL'] = (humancombo['read.pvalue'] <= 0.05) & (humancombo['gain.vs.loss'] == 'loss')
        humancombo1 = humancombo[(humancombo['PARENTAL_PLOIDY']>4)&(humancombo['PLOIDY']>4)]
        count = (lambda C, V : C[V] / C.sum() if V in C else 0.0)
        form = {'FREQ_AMP' : ('SIGN_AMP', (lambda c : count(c.value_counts(), True))), 
                'FREQ_DEL' : ('SIGN_DEL', (lambda c : count(c.value_counts(), True)))}
        humanfreq = humancombo1.groupby(['GENE_NAME'])[[form[k][0] for k in form]].agg(**form).reset_index()

        humanratiod = humancombo1.groupby('CELLLINE')['SIGN_DEL'].apply(lambda v : v.value_counts()[True] / len(v) if True in v.value_counts() else 0.0).to_numpy()
        humanratiod = sorted(np.sum(scipy.stats.bernoulli.rvs(humanratiod)) / len(humanratiod) for x in range(100000))
        humanratiod = min(humanratiod[-int(len(humanratiod) * 0.05):])
        humanfreq['SEL_DEL'] = humanfreq['FREQ_DEL'] >= humanratiod

        humanratiog = humancombo1.groupby('CELLLINE')['SIGN_AMP'].apply(lambda v : v.value_counts()[True] / len(v) if True in v.value_counts() else 0.0).to_numpy()
        humanratiog = sorted(np.sum(scipy.stats.bernoulli.rvs(humanratiog)) / len(humanratiog) for x in range(100000))
        humanratiog = min(humanratiog[-int(len(humanratiog) * 0.05):])
        humanfreq['SEL_AMP'] = humanfreq['FREQ_AMP'] >= humanratiog

        freq = humanfreq.set_index('GENE_NAME').join(mousefreq.set_index('GENE_NAME'), lsuffix='_human', rsuffix='_mouse').reset_index().dropna()
        freqgen = freq.merge(right=humancombo[['GENE_NAME', 'seqname', 'start', 'end']], on='GENE_NAME', how='left').drop_duplicates()
        freqgen['CHR'] = freqgen['seqname'].str.replace('chr', '').astype(int)
        freqgen = freqgen.sort_values(['CHR', 'start', 'end'])
        freqgen['end_chr'] = freqgen['CHR'].map(freqgen.groupby('CHR')['end'].max().cumsum() - freqgen.groupby('CHR')['end'].max())
        freqgen['START'] = freqgen['start'] + freqgen['end_chr']
        freqgen['END'] = freqgen['end'] + freqgen['end_chr']
        freqgen['SEL_DEL'] = freqgen['SEL_DEL_human'] & freqgen['SEL_DEL_mouse'] & freqgen['GENE_NAME'].isin(cosmic)
        freqgen['SEL_AMP'] = freqgen['SEL_AMP_human'] & freqgen['SEL_AMP_mouse'] & freqgen['GENE_NAME'].isin(cosmic)

        return {
            'freqgen' : freqgen,
            'humanratiod' : humanratiod,
            'mouseratiod' : mouseratiod,
            'humanratiog' : humanratiog,
            'mouseratiog' : mouseratiog
        }


def freq_losses(freqgen, humanratiod, mouseratiod, **_):
    fig, ax = plt.subplots(2, 1, figsize=(20, 12), sharex=True, sharey=True)

    nonselgen = freqgen[~freqgen['SEL_DEL']].reset_index(drop=True)
    ax[0].plot((nonselgen['START'], nonselgen['END']), (nonselgen['FREQ_DEL_human'], nonselgen['FREQ_DEL_human']), c='k', linewidth=3)
    for x1, x2, y in zip(nonselgen['START'], nonselgen['END'], nonselgen['FREQ_DEL_human']):
        ax[0].fill_between((x1, x2), (y, y), (0, 0), color='k', alpha=0.1)
    ax[1].plot((nonselgen['START'], nonselgen['END']), (nonselgen['FREQ_DEL_mouse'], nonselgen['FREQ_DEL_mouse']), c='k', linewidth=3)
    for x1, x2, y in zip(nonselgen['START'], nonselgen['END'], nonselgen['FREQ_DEL_mouse']):
        ax[1].fill_between((x1, x2), (y, y), (0, 0), color='k', alpha=0.1)
        
    selqgen = freqgen[freqgen['SEL_DEL']].reset_index(drop=True)
    ax[0].plot((selqgen['START'], selqgen['END']), (selqgen['FREQ_DEL_human'], selqgen['FREQ_DEL_human']), c='blue', linewidth=8)
    for x1, x2, y in zip(selqgen['START'], selqgen['END'], selqgen['FREQ_DEL_human']):
        ax[0].fill_between((x1, x2), (y, y), (0, 0), color='blue', alpha=0.6)

    c1pos = 0
    start = -1
    prey = -1
    schr = -1
    text = ""
    for c, x1, x2, y, n in zip(selqgen['CHR'], selqgen['START'], selqgen['END'], selqgen['FREQ_DEL_human'], selqgen['GENE_NAME']):
        x = (x1 + x2) / 2
        if start == -1:
            start = x
            schr = c
            prey = y
        if x - start > 20000000 or c != schr:
            if schr in [2, 3, 13, 14, 15, 16, 19]:
                c1pos += 1
            cox = 0
            coy = 25
            if c1pos == 1 and schr == 2:
                cox = -80
                coy = 25
            elif c1pos == 2 and schr == 2:
                cox = -25
                coy = 25
            elif c1pos == 3 and schr == 2:
                cox = 20
                coy = 25
            elif c1pos == 4 and schr == 2:
                cox = 50
                coy = 25
            elif c1pos == 5 and schr == 14:
                cox = -70
                coy = 40
            elif c1pos == 6 and schr == 14:
                cox = -50
                coy = 65
            elif c1pos == 7 and schr == 15:
                cox = -50
                coy = 90
            elif c1pos == 8 and schr == 15:
                cox = 0
                coy = 70
            elif c1pos == 9 and schr == 16:
                cox = 0
                coy = 30
            elif c1pos == 10 and schr == 19:
                cox = -50
                coy = 25
            elif c1pos == 11 and schr == 19:
                cox = 0
                coy = 25
            elif c1pos == 12 and schr == 19:
                cox = 18
                coy = -5
            ax[0].annotate(text.strip(),
                    xy=(start, prey), xycoords='data',
                    xytext=(cox, coy), textcoords='offset points',
                    arrowprops=dict(facecolor='blue', shrink=0.05, linewidth=1),
                    horizontalalignment='left', verticalalignment='bottom')
            text = "{}\n".format(n)
            start = x
            schr = c
            prey = y
        else:
            text += "{}\n".format(n)
    ax[0].annotate(text,
        xy=(start, prey), xycoords='data',
        xytext=(cox, coy), textcoords='offset points',
        arrowprops=dict(facecolor='blue', shrink=0.05, linewidth=1),
        horizontalalignment='left', verticalalignment='bottom')
            
    ax[1].plot((selqgen['START'], selqgen['END']), (selqgen['FREQ_DEL_mouse'], selqgen['FREQ_DEL_mouse']), c='blue', linewidth=8)
    for x1, x2, y in zip(selqgen['START'], selqgen['END'], selqgen['FREQ_DEL_mouse']):
        ax[1].fill_between((x1, x2), (y, y), (0, 0), color='blue', alpha=0.6)

    ax[0].set_ylabel('Frequency of losses')
    ax[0].set_title('Lost genes in hexaploid PC9 or with hexaploid parents')
    ax[1].set_ylabel('Frequency of losses')
    ax[1].set_title('Lost genes in EP treated single mouse cells')
    ticks = freqgen[freqgen['CHR'] != freqgen['CHR'].shift()]
    ax[0].plot((ticks['START'].min(), ticks['START'].max()), (humanratiod, humanratiod), '--g')
    ax[1].plot((ticks['START'].min(), ticks['START'].max()), (mouseratiod, mouseratiod), '--g')
    plt.xticks(ticks=ticks['START'], labels=ticks['CHR'])
    plt.xlim(ticks['START'].min(), ticks['START'].max())
    plt.ylim(ymin=0, ymax=1)
    plt.savefig('PC9_freq_losses.png', bbox_inches='tight', dpi=300)
    selqgen.to_csv('PC9-losses-frequencies.tsv.gz', sep='\t', index=False)


def freq_gains(freqgen, humanratiog, mouseratiog, **_):
    fig, ax = plt.subplots(2, 1, figsize=(20, 12), sharex=True, sharey=True)

    nonselgen = freqgen[~freqgen['SEL_AMP']].reset_index(drop=True)
    ax[0].plot((nonselgen['START'], nonselgen['END']), (nonselgen['FREQ_AMP_human'], nonselgen['FREQ_AMP_human']), c='k', linewidth=3)
    for x1, x2, y in zip(nonselgen['START'], nonselgen['END'], nonselgen['FREQ_AMP_human']):
        ax[0].fill_between((x1, x2), (y, y), (0, 0), color='k', alpha=0.1)
    ax[1].plot((nonselgen['START'], nonselgen['END']), (nonselgen['FREQ_AMP_mouse'], nonselgen['FREQ_AMP_mouse']), c='k', linewidth=3)
    for x1, x2, y in zip(nonselgen['START'], nonselgen['END'], nonselgen['FREQ_AMP_mouse']):
        ax[1].fill_between((x1, x2), (y, y), (0, 0), color='k', alpha=0.1)
        
    selqgen = freqgen[freqgen['SEL_AMP']].reset_index(drop=True)
    ax[0].plot((selqgen['START'], selqgen['END']), (selqgen['FREQ_AMP_human'], selqgen['FREQ_AMP_human']), c='red', linewidth=8)
    for x1, x2, y in zip(selqgen['START'], selqgen['END'], selqgen['FREQ_AMP_human']):
        ax[0].fill_between((x1, x2), (y, y), (0, 0), color='red', alpha=0.6)

    c1pos = 0
    start = -1
    prey = -1
    schr = -1
    text = ""
    for c, x1, x2, y, n in zip(selqgen['CHR'], selqgen['START'], selqgen['END'], selqgen['FREQ_AMP_human'], selqgen['GENE_NAME']):
        x = (x1 + x2) / 2
        if start == -1:
            start = x
            schr = c
            prey = y
        if x - start > 10000000 or c != schr:
            if schr in [3, 7, 10, 12, 20]:
                c1pos += 1
            cox = 0
            coy = 25
            if c1pos == 1 and schr == 3:
                cox = -80
                coy = 10
            if c1pos == 2 and schr == 3:
                cox = -100
                coy = 50
            if c1pos == 3 and schr == 3:
                cox = -50
                coy = 40
            if c1pos == 4 and schr == 3:
                cox = -20
                coy = 80
            if c1pos == 5 and schr == 3:
                cox = 40
                coy = 30
            if c1pos == 6 and schr == 3:
                cox = 100
                coy = 15
            if c1pos == 7 and schr == 7:
                cox = -80
                coy = 35
            if c1pos == 8 and schr == 7:
                cox = -40
                coy = 45
            if c1pos == 9 and schr == 7:
                cox = 0
                coy = 60
            if c1pos == 10 and schr == 7:
                cox = 30
                coy = 25
            if c1pos == 11 and schr == 10:
                cox = -50
                coy = 45
            if c1pos == 12 and schr == 10:
                cox = 0
                coy = 45
            if c1pos == 13 and schr == 10:
                cox = 20
                coy = 25
            if c1pos == 14 and schr == 12:
                cox = -50
                coy = 55
            if c1pos == 15 and schr == 12:
                cox = 0
                coy = 100
            if c1pos == 16 and schr == 12:
                cox = 30
                coy = 40
            if c1pos == 17 and schr == 12:
                cox = 80
                coy = 30
            if c1pos == 18 and schr == 20:
                cox = -100
                coy = 40
            if c1pos == 19 and schr == 20:
                cox = -60
                coy = 60
            if c1pos == 20 and schr == 20:
                cox = -20
                coy = 80
            ax[0].annotate(text.strip(),
                    xy=(start, prey), xycoords='data',
                    xytext=(cox, coy), textcoords='offset points',
                    arrowprops=dict(facecolor='red', shrink=0.05, linewidth=1),
                    horizontalalignment='left', verticalalignment='bottom')
            text = "{}\n".format(n)
            start = x
            schr = c
            prey = y
        else:
            text += "{}\n".format(n)
    ax[0].annotate(text,
        xy=(start, prey), xycoords='data',
        xytext=(cox, coy), textcoords='offset points',
        arrowprops=dict(facecolor='red', shrink=0.05, linewidth=1),
        horizontalalignment='left', verticalalignment='bottom')
            
    ax[1].plot((selqgen['START'], selqgen['END']), (selqgen['FREQ_AMP_mouse'], selqgen['FREQ_AMP_mouse']), c='red', linewidth=8)
    for x1, x2, y in zip(selqgen['START'], selqgen['END'], selqgen['FREQ_AMP_mouse']):
        ax[1].fill_between((x1, x2), (y, y), (0, 0), color='red', alpha=0.6)

    ax[0].set_ylabel('Frequency of gains')
    ax[0].set_title('Gained genes in hexaploid PC9 or with hexaploid parents')
    ax[1].set_ylabel('Frequency of gains')
    ax[1].set_title('Gained genes in EP treated single mouse cells')
    ticks = freqgen[freqgen['CHR'] != freqgen['CHR'].shift()]
    ax[0].plot((ticks['START'].min(), ticks['START'].max()), (humanratiog, humanratiog), '--g')
    ax[1].plot((ticks['START'].min(), ticks['START'].max()), (mouseratiog, mouseratiog), '--g')
    plt.xticks(ticks=ticks['START'], labels=ticks['CHR'])
    plt.xlim(ticks['START'].min(), ticks['START'].max())
    plt.ylim(ymin=0, ymax=1)
    plt.savefig('PC9_freq_gains.png', bbox_inches='tight', dpi=300)
    selqgen.to_csv('PC9-gains-frequencies.tsv.gz', sep='\t', index=False)


