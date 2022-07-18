# stMLnet

## What is stMLnet?
![image](https://github.com/SunXQlab/stMLnet/blob/main/overview_stMLnet.png)
`stMLnet` is a tool to infer spatial intercellular communication and multilayer signaling regulations from `spatial transcriptomic data (ST)` by quantifying distance-weighted ligand–receptor (LR) signaling activity based on diffusion and mass action models and mapping it to intracellular targets. stMLnet can infer, quifity, and visualize both intercelluar communitations and the intracelluar gene regulatory network from ST data. stMLnet allows:
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

## Usage
To learn how to use this tool, chek `Tutorial of stMLnet.Rmd`.

## Examples and Reproducibility
All the examples and the reproducibility codes for the plots in the paper could be found in the `stMLnet-AnalysisCode` repository which includes:

* `prior_knowledge` contains the code used for collection and integration of prior databases <br>
* `apply_in_simu` contains the code to reproduce the simulation study of stMLnet <br>
* `apply_in_stBC` contains the code to reproduce plots and benchmarking of the breast cancer dataset <br>
* `apply_in_COVID19` contains the code to reproduce plots and detailed analysis of the COVID-19 ST dataset <br>
* `apply_in_scGNM` contains the code to reproduce plots and detailed analysis for appling stMLnet on the scRNA-seq dataset of glioma <br>
* `code` contains all functions of stMLnet to analysis cell-cell interactions <br>

See detials therein.


## Support and contribute
If you have questions or suggestions for imrpoving stMLnet, please contact Xiaoqiang Sun via sunxq6@mail.sysu.edu.cn.





