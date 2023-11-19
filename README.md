# TP53loss_WGD

This repository includes the code needed to reproduce the analysis and figures in the manuscript:

[TP53 loss with whole genome doubling mediates intra - patient heterogeneous 2 therapy response s through Chromosomal Instability]()

Pre-executed notebooks and demos for the different analyses are available for:
- [Single-cell resistance analysis](scripts/sc_analysis.ipynb)
- [PC9 bulk analysis](scripts/PC9_analysis.ipynb)
- [PC9 human-mouse gene frequency analysis](scripts/freq_analysis.ipynb)

In addition the CNA calling pipeline from PC9 cell lines is [available](cnacalling_PC9/copy.number.detection.cellline.R).

## Requirements

Executing the code in this repository requires to first download the publicly available pre-process data:

[Zenodo link to appear soon]()

Moreover, the code is implemented in Python 3 and a `conda` enviroment ([miniconda](https://docs.conda.io/projects/miniconda/en/latest/) reccomended) with all the requirements can be created with the following command (use of [mamba](https://mamba.readthedocs.io/en/latest/) is reccomended):

```
conda create -n tp53loss_wgd python \
                             numpy \
                             pandas=1.3.5 \
                             matplotlib-base \
                             seaborn=0.11 \
                             scipy \
                             pingouin \
                             statannot \
                             gtfparse=1.2.1 \
                             jupyterlab \
                             openpyxl
```

The environment has to be activated before each use with
```
conda activate tp53loss_wgd
```
