Figure 1A, 2A - Biorender
Figure 1B, 1D, 1E, 2C, 3B - methatlas + edits for sample labels and 'other' category
Figure 1C, Supplementary Figure 1A-F - synth_downsample_plots.py
Figure 1F - Excel
Figure 1G, 2E, 3D, Supp.2A-E, Supp.3A-C - outliers.R + edits for sample labels
Figure 2B, 3A - IGV + edits for variant annotation
Figure 2D, 3C - plot_deconv_proportions.py
Figure 3E - plot_contrasting_bioheart_proportions.py



#### outliers.R: ####

Input data should be a matrix with 
- Columns of samples
- Rows of loci
- Entries of methylation frequencies at the loci


#### plot_contrasting_bioheart_proportions ####

Input data should be 2 files as described in the lead comment

- One of columns samples, rows cell types, with entries cell type contributions
- The other columns metadata, rows sample, celltype, and contribution

- edit the lines corresponding to mapping the cell type names between
the files as appropriate

#### plot_deconv_proportions ####

Input data should essentially be a methatlas deconv output file
- columns of samples
- rows of celltypes
- entries of cell type contributions to the sample

#### synth_downsample_plots #### 

- Input data: a series of files, each a repeat of the downsample experiment,
where each file is a methatlas output csv
- columns of samples
- rows of celltypes
- entries of cell type contributions
- column and row names should match between files
