# scFseCluster code

scFseCluster provides four python files and one R file. We have added good comments to the code to help those who use our proposed method. Each of the five files performs the following functions.

- `scFseCluster.py`: This file is the entry point for the execution of the scFseCluster program.
- `Preprocessing.py`: In this file, scFseCluster will help us to read and preprocess the scRNA-seq data.
- `FSQSSA.py`: Core program. Here the FSQSSA algorithm will be iterated `max_iter` times (default 100 times). It will work in concert with `Cluster.py`.
- `Cluster.py`: The feature subset corresponding to each squirrel is obtained and the fitness value is calculated.
- `plot_scFseCluster.R`: Visualize cell clustering results.