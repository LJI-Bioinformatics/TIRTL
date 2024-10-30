# TIRTL-seq data analysis pipeline
Code for paired TCR data analysis from TIRTL-seq experiments. Please see the [preprint](https://www.biorxiv.org/content/10.1101/2024.09.16.613345v1) for details. Data for the manuscript is available at Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14010377.svg)](https://doi.org/10.5281/zenodo.14010377).

## Installation
### Requirements

For Nvidia GPU functionality you will need to have CUDA installed. You can install it from [here](https://developer.nvidia.com/cuda-downloads).
For Apple Silicon GPU install Apple MLX library from [here](https://github.com/ml-explore/mlx).



R dependencies are _Matrix_ and _data.table_ packages, available by default with most distributions. If not, you can install them by running from R console:
```R
install.packages("Matrix")
install.packages("data.table")
```

Python3 dependencies are _numpy_, _pandas_, _cupy_ (for Nvidia GPU runs, optional), _mlx_ (for Apple GPU runs, optional). You can install them by running:
```bash
pip install numpy pandas
pip install cupy
pip install mlx
```
## Quick Start
Install dependencies, download repository and run run the following commands in R:
```R
source("TIRTL_functions.R")
# to run on cpu
test<-run_single_point_analysis_sub_gpu("data/")
# to run on Nvidia gpu
test_gpu_nvidia<-run_single_point_analysis_sub_gpu("data/",backend="cupy")
# to run on Apple Silicon gpu
test_gpu_apple<-run_single_point_analysis_sub_gpu("data/",backend="mlx")
# to run on 192 wells, first 12 columns of 384-well plate
test_half_plate<-run_single_point_analysis_sub_gpu("data/",wellset1=get_well_subset(1:16,1:12))
# now lets run on the small dataset from the manuscript (download and unpack data from zenodo https://doi.org/10.5281/zenodo.14010377 first!)
cd8_tp2<-run_single_point_analysis_sub_gpu("exp3_clones/TCR_clones_ID03/",wellset1=get_well_subset(1:16,1:12),backend="numpy") # CD8 repertoire for time point 2 for COVID patient (left half of 384-well plate)
print(cd8_tp2) # we got alpha-beta pairs!
table(cd8_tp2$method) # note that same alpha-beta pairs can be called by different methods, so there are duplicates
```

### Input for pairing pipeline
The input data should be a folder with TCRalpha and TCRbeta repertoires from one plate in separate tab separated files from standard _mixcr_ output.

Other parameters are optional and can be set to change the default behavior of the pipeline:

_prefix_ (default="tmp") - prefix for the output files

_well_filter_thres_ (default=0.5 of mean number of unique clones across alpha wells) - threshold for the number unique clones in a well to be considered as a valid well. Decreasing this threshold will increase the number of wells considered as valid, but may also increase the number of false positives.

_min_reads_ (default=0) - minimum number of reads a clone to be included. By default all clones are included. Increase this number if you suspect cross-contamination between wells or high number of sequencing errors. 

_well_pos_ (default=3) - position of the well id in the file name if split by underbar. For example, if the file name is "3123233_MVP093_TIRTLseq_**A1**_S8.alpha_old_ID10.clones_TRA.tsv", the well_pos should be 4.

_wellset1_ (default=get_well_subset(1:16,1:24)) - subset of wells to be analyzed (if you loaded your sample to part of the plate). By default all wells in 384 well plate are analyzed. For example, to analyze the first 12 columns of a 384-well plate, you can use get_well_subset(1:16,1:12).

_compute_ (default=T) - whether to run backend script for pairing. Only set F if you want to run the python backend script separately.

_backend_ (default="numpy") - backend for the pairing. Can be "numpy" for CPU, "cupy" for Nvidia GPU, or "mlx" for Apple Silicon GPU.

### Output for pairing pipeline
This will take TCRalpha and TCRbeta repertoires from the folder "data/" and pair them using the default parameters. The output will be saved in the working directory under _tmp_TIRTLoutput.tsv_ and also returned to R as _data.table_.

_wi_ - number of wells where alpha is found, but not beta

_wj_ - number of wells where beta is found, but not alpha

_wij_ - number of wells where both alpha and beta are found

_wa_ - number of wells where alpha is found

_wb_ -number of wells where beta is found

_alpha_nuc_seq_, _alpha_nuc_ - alpha CDR3nuc sequence

_beta_nuc_seq_, _beta_nuc_ - beta CDR3nuc sequence

_alpha_beta_ - alpha-beta pair concatenated nucleotide sequence

_method_ - method used for pairing t-shell, or MAD-HYPE

_r_ - for t-shell, the Pearson correlation coefficient between alpha and beta frequencies

_ts_ - for t-shell, the t-statistic of the correlation

_pval_ - for t-shell, the p-value of the correlation

_pval_adj_ - for t-shell, the adjusted p-value of the correlation

_loss_a_frac_ - the fraction of wells with lost alpha chain

_loss_b_frac_ - the fraction of wells with lost beta chain

_score_ - for MAD-HYPE, the score of the pairing (higher is better)

_cdr3a_ - alpha CDR3aa sequence

_va_ - alpha V gene

_ja_ - alpha J gene

_cdr3b_ - beta CDR3aa sequence

_vb_ - beta V gene

_jb_ - beta J gene



