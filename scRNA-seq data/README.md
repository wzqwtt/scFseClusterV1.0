# scRNA-seq data

The directory contains the dataset for the scFseCluster test. In our experiments, each dataset contains three files.

- `*_cell_label.csv`: Clustering labels for each cell. 
- `raw_data.csv`: scRNA-seq expression matrix; rows represent cells, columns represent genes. 
- `T2000_expression.txt`: Top 2000 highly variable genes (HVGs). Not necessary, scFseCluster can help you generate HVGs. rows represent genes, columns represent cells. 

In use, you can enter only the expression matrix of scRNA-seq. 

More datasets are available for download on our website at http://cdsic.njau.edu.cn/data/scFseCluster.
