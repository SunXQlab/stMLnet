# stMLnet

## what is stMLnet?
![image](https://github.com/SunXQlab/stMLnet/blob/main/overview_stMLnet.png)
`stMLnet` is a tool to infer spatial intercellular communication and multilayer signaling regulations from `spatial transcriptomic data (ST)` by quantifying distance-weighted ligand–receptor (LR) signaling activity based on diffusion and mass action models and mapping it to intracellular targets. stMLnet can infer, quifity, and visualize both intercelluar communitations and the intracelluar gene regulatory work from ST data. stMLnet allows:
* constructs a multilayer signaling network, infer LR signling activate, and predicte LR-target gene regulation <br>
* leverages spatial information in the ST data to quantify intercellular signaling activity and connect extracellular signals to intracellular gene expression<br>
* provides various visualizations to vividly depict mechanisms underlying cellular communications and molecular regulations <br>

## Package Structure
The repository is centered around the `runMLnet` module:
* `runMLnet` contains the scripts to create mulitlayer signling network <br>
* `runImputation` contains the scripts to imputate the ST normalized expression matrix <br>
*  `getSiganlActivity` contains the scripts to obtain the upstream paired signaling activity and downstream target gene expression based on the mulitlayer signling network <br>
*  


