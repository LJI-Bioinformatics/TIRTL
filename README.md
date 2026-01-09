# TIRTL-seq data analysis pipeline
Code for paired TCR data analysis from TIRTL-seq experiments. Please see the [paper](https://www.nature.com/articles/s41592-025-02907-9) for details. Data for the manuscript is available at Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14010377.svg)](https://doi.org/10.5281/zenodo.14010377). More functions for TIRTL-seq data analysis and visualization are available in actively developing [TIRTLtools](https://github.com/NicholasClark/TIRTLtools) R package.

## Installation

Install from Github with:
```{r}
devtools::install_github("LJI-Bioinformatics/TIRTL")
```

### Additional requirements

For Nvidia GPU functionality you will need to have CUDA installed. You can install it from [here](https://developer.nvidia.com/cuda-downloads).
For Apple Silicon GPU install Apple MLX library from [here](https://github.com/ml-explore/mlx).

Python3 dependencies are _numpy_ (v.>1.26), _pandas_ (v.>2.2), _cupy_ (v. >13.2 for Nvidia GPU runs, optional), _mlx_ (v. >=0.13 for Apple GPU runs, optional). You can install them by running:
```bash
pip install numpy pandas
pip install cupy
pip install mlx
```

The code was tested on macOS 14-15 (CPU/GPU), Red Hat Enterprise Linux 8.8 (CPU/Nvidia GPU) and Windows 11 Enterprise (tested on CPU only).

## Quick Start

Install dependencies, download repository and run the following commands in R:

```R
library(TIRTL)
# to run on cpu
test<-run_single_point_analysis_sub_gpu("data/")
# to run on Nvidia gpu
test_gpu_nvidia<-run_single_point_analysis_sub_gpu("data/",backend="cupy")
# to run on Apple Silicon gpu
test_gpu_apple<-run_single_point_analysis_sub_gpu("data/",backend="mlx")
# to run on 192 wells, first 12 columns of 384-well plate
test_half_plate<-run_single_point_analysis_sub_gpu("data/",wellset1=get_well_subset(1:16,1:12))
# DEMO! lets run on the small dataset from the manuscript (download and unpack data from zenodo https://doi.org/10.5281/zenodo.14010377 first!)
# This test takes ~6 minutes on CPU on a 2021 Macbook Pro (M1 Pro 32 GB RAM)
cd8_tp2<-run_single_point_analysis_sub_gpu("exp3_clones/TCR_clones_ID03/",wellset1=get_well_subset(1:16,1:12),backend="numpy") # CD8 repertoire for time point 2 for COVID patient (left half of 384-well plate)
print(cd8_tp2) # we got alpha-beta pairs!
table(cd8_tp2$method) # note that same alpha-beta pairs can be called by different methods, so there are duplicates
#expected output: 
#madhype  tshell 
#  16369    9891 
```

### Input for pairing pipeline
The input data should be a folder with TCRalpha and TCRbeta repertoires from one plate in separate tab separated files from standard _mixcr_ output. Filenames should contain alphanumeric well id in the format "A1", "B2", etc separated by underbars.

Other parameters are optional and can be set to change the default behavior of the pipeline:

_prefix_ (default="tmp") - prefix for the output files

_well_filter_thres_ (default=0.5 of mean number of unique clones across alpha wells) - threshold for the number unique clones in a well to be considered as a valid well. Decreasing this threshold will increase the number of wells considered as valid, but may also increase the number of false positives.

_min_reads_ (default=0) - minimum number of reads a clone to be included. By default all clones are included. Increase this number if you suspect cross-contamination between wells or high number of sequencing errors. 

_well_pos_ (default=3) - position of the well id in the file name if split by underbar. For example, if the file name is "3123233_MVP093_TIRTLseq_**A1**_S8.alpha_old_ID10.clones_TRA.tsv", the well_pos should be 4.

_wellset1_ (default=get_well_subset(1:16,1:24)) - subset of wells to be analyzed (if you loaded your sample to part of the plate). By default all wells in 384 well plate are analyzed. For example, to analyze the first 12 columns of a 384-well plate, you can use get_well_subset(1:16,1:12).

_compute_ (default=T) - whether to run backend script for pairing. Only set F if you want to run the python backend script separately.

_backend_ (default="numpy") - backend for the pairing. Can be "numpy" for CPU, "cupy" for Nvidia GPU, or "mlx" for Apple Silicon GPU.

_wij_thres_tshell_ (default=2) - threshold for the number of wells with both alpha and beta chains for T-SHELL method. Increasing this number will increase the stringency of the pairing, but may also decrease the number of pairs found. Recommended value for 96-well experiments is 3.

_pval_thres_tshell_ (default=1e-10) - p-value threshold for the T-SHELL method. Recommended value for 96-well experiments is 1e-3.

_pseudobulk_ (default=F) - this option generates pseudobulk files, a single file containing all merged TCRalpha and TCRbeta singlechain data aggregated across wells. Resulting files are written out with _pseudobulk_TRA.tsv and _pseudobulk_TRB.tsv suffixes.

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

### Running on your own data
To run the pipeline on your own data, you need to have TCRalpha and TCRbeta repertoires processed with _mixcr_ (version 4.6) separate for each well. 

We use the following mixcr pipeline for each well: 
```bash
$MIXCR analyze generic-amplicon --species hsa --rna \
--tag-pattern "^N{0:1}(SAMPLE:N{7})(R1:*)\^(R2:*)" \
--floating-left-alignment-boundary \
--floating-right-alignment-boundary C \
--split-by-sample \
--sample-sheet plate_index.tsv \
${fp}_L{{n}}_R1.fastq.gz ${fp}_L{{n}}_R2.fastq.gz \
$OUTDIR/${sample} \
--use-local-temp
```
where _plate_index.tsv_ is a tab separated file containing you PCR I plate indices (see example in this repo!). 

_MIXCR_ variable contains path to _mixcr_ executable.

_fp_ is the path to the fastq files for given wells without the _L{n}_R1.fastq.gz and _L{n}_R2.fastq.gz suffixes.

_OUTDIR_ is the output directory for the mixcr output.

_sample_ is the file name for the output for the given well (usually only filename from _fp_).


