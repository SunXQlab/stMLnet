# stMLnet

## What is stMLnet?
<p align="center">
  <img src="https://github.com/SunXQlab/stMLnet/blob/main/fig1.framework.png">
</p>
`stMLnet` is a tool to infer spatial intercellular communication and multilayer signaling regulations from `spatial transcriptomic data (ST)` by quantifying distance-weighted ligand–receptor (LR) signaling activity based on diffusion and mass action models and mapping it to intracellular targets. stMLnet can infer, quifity, and visualize both intercelluar communitations and the intracelluar gene regulatory network from ST data. stMLnet allows:<br>
* to construct a multilayer signaling network, infer LR signling activate, and predicte LR-target gene regulation <br>
* to leverage spatial information in the ST data to quantify intercellular signaling activity and connect extracellular signals to intracellular gene expression<br>
* to visualize inter- and intra-cellular signaling networks and functions associated with cellular communications and molecular regulations <br>

We also provide the R code used for collection and integration of prior databases, detailed in `stMLnet-AnalysisCode` repository

## Package Structure
The repository is centered around the `R` module:
* `creat_multilayer_network` contains the scripts to create mulitlayer signling network <br>
* `calculate_signal_activity` contains the scripts to obtain the upstream paired signaling activity and downstream target gene expression based on the mulitlayer signling network <br>
* `calculate_signal_importance` contains the scripts to calculate the upstream signal pairs or signals importance in the multilay signal network of cell communication <br>
* `visualize_cell_communication` contains the scripts to visualize cell-cell interations <br>
* `check_feedback_loop` contains the scripts to select multicellular feedback circuits <br>

## Usage

To install and use the stMLnet package, please make sure you have already installed related dependencies:

       # Check if the following dependencies are installed.
       pkgs <- c('Seurat','Giotto','reshape2','stringr','dplyr', # for data preprocessing
                        'caret','doParallel','snow','foreach', # for quantitative model
                         'ggplot2','ggsci','org.Hs.eg.db', 'plotrix','ggalluvial','ggraph','igraph' # for visualization
                         )
       for (pkg in pkgs) {
         if (!requireNamespace(pkg)) { cat(paste0('please install and library the package: ',pkg,'\n')) }
       }
       
       # Installing related dependencies.
       pkgs <- c( 'caret','doParallel','snow','foreach','ggplot2','ggsci','org.Hs.eg.db','plotrix','ggalluvial','ggraph','igraph')
       for (pkg in pkgs) {install.packages(pkg, repos = 'https://cloud.r-project.org')}
       
       devtools::install_version("Seurat", version = "4.2.0", repos="https://cloud.r-project.org/")
       remotes::install_github("satijalab/seurat-wrappers")
       remotes::install_github("drieslab/Giotto",  ref="v1.1.0")

After building dependent environment, you can download stMLnet from github:

       git clone https://github.com/SunXQlab/stMLnet.git
       
and then install stMLnet from local:
 
       install.packages("path/to/stMLnet/stMLnet_0.1.2.tar.gz", repos = NULL, type = "source")
       library(stMLnet)
       
or you can directly install stMLnet from github:

       devtools::install_github("SunXQlab/stMLnet")
       library(stMLnet)

If you have problems installing the environment manually, you can also choose to install the dependent environment via dockfile:

       # Bash
       # built a docker image
       # ensure that dockerfile and postInstall are in the same path
       docker build -f Dockerfile -t stMLnetEnv:0.1 .
       # Run docker image
       docker run -it stMLnetEnv:0.1 /bin/bash

To learn how to use this tool, check [Tutorial of stMLnet.Rmd](https://github.com/SunXQlab/stMLnet/blob/main/Tutorial%20of%20stMLnet.Rmd). This tutorial shows the installation and application of stMLnet in the demo dataset, which can be download from [here](https://www.zenodo.org/record/8024150). It will take about 15 mins to run this demo (excluding environment installation) mainly depending on the parameter setting in the quantitative analysis step.

## Examples and Reproducibility

All the examples and the reproducibility codes for the plots in the paper could be found in the `stMLnet-AnalysisCode` repository which includes:

* `prior_knowledge` contains the code used for collection and integration of prior databases <br>
* `apply_in_simu` contains the code to reproduce the simulation study of stMLnet <br>
* `apply_in_scST` contains the code to reproduce the plot and detailed analysis of the three single-cell resolution ST datasets <br>
* `apply_in_COVID19` contains the code to reproduce plots and detailed analysis of the COVID-19 ST dataset <br>
* `benchmark` contains the code to reproduce plots and benchmarking <br>
* `code` contains all functions of stMLnet to analysis cell-cell interactions <br>

See detials therein.

We also provide a web-based application to demonstrate the functionality and visualization of stMLnet, available at www.stmlnet.top/net.

## Contact
If you have questions or suggestions for imrpoving stMLnet, please contact Xiaoqiang Sun via sunxq6@mail.sysu.edu.cn.





