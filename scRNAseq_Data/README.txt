Each folder in the scRNA-seq_Data directory contains the 3 files outputted by CellRanger's pipeline: barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz. 
These files can be downloaded to a local machine separately or as a folder.
It is recommended to download the each feature_bc_matrix folder, so it can be imported into R as a Seurat object.
For example, run Read10x(data.dir= "~/22-0160_feature_bc_matrix") where ~ is the local location of the folder 
