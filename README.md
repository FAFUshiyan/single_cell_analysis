# single cell analysis pipeline(Seurat)
## annotate_clusters.R

Run this script to annotate Seurat clusters with user-defined labels and generate gene expression plots.

### Usage
`Rscript annotate_clusters.R \`

`--anno_file cluster_map_B.txt \`

`--rds_file  Result_V1/results/rds/Bsub_recluster.rds \`

`--out_dir   Result_V1/results/Bsub_out \`

`--prefix    Bcells \`

`--genes     SOD1`

### Description
This command will:

Read a Seurat object from the provided .rds file.

Apply cluster annotations based on your mapping file.

Produce per-cluster violin plots and UMAP feature plots for one or more genes.

Write all outputs (plots and annotated metadata) into the specified output directory, using your given prefix.

### Arguments
--anno_file
Path to a text file mapping cluster IDs to human-readable labels. Each line must be:

php-template

`<cluster_id>,<annotation_label>`
For example:

`29,CD8 MAIT`
`34,CD4 T NAIVE`

--rds_file
Path to the input Seurat object (saved via saveRDS()).

--out_dir
Directory where annotated metadata and plots will be saved.

--prefix
Filename prefix for all outputs (e.g., Bcells will produce files like Bcells_violin_PLAAT4.png).

--genes
Comma-separated list of gene symbols to visualize.  
  - **Multiple genes** (e.g. `--genes CD8A,SOD1`): generates a single multi-gene violin plot and one multi-panel UMAP feature plot  with automatic layout.  
  - **Single gene**: falls back to the original behavior, producing individual gene plots (with `split.by` support).
A violin plot (VlnPlot) showing expression distributions across annotated clusters 

A UMAP feature plot (FeaturePlot) mapping expression onto the UMAP embedding 

You can request multiple genes at once, e.g. --genes PLAAT4,CD8A,SOD1.

### Requirements
Rscript: the command-line front end for running R scripts without interactive startup 

Seurat (v5 or later): for single-cell object handling and plotting.