# Synthetic Samples

The relevant configs and scripts have been copied here. Note, though, that none of the filepaths
will be valid.

General usage is 
```
python synthetic_samples_methatlas.py <config.yaml>
```

where a config file specifies the different bamfiles to sample and with what weightings
and ratios

Additionally, beta files can now be directly converted to Illumina Indices with WGBStools so the
- Beta to Blocks
- Beta/Blocks to Coverage
- CpG Coverage to Illumina Coverage
sections are no longer necessary. Previously, this used a mapping between hg38 chromosomal loci
and Illumina indices. Something along the lines of 

```
os.system(f"wgbstools beta_to_450k subsampled.sorted.beta")
```
can be used instead, with some post-processing for appropriate formating.

The mapping file used is included here. Note that it only includes mapping for the (~8000) sites used by
methatlas, not all Illumina sites.
