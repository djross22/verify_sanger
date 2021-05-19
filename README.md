Python package for automated analysis of Sanger sequencing data to verify construct insters into plasmids.

Written for Python 3.

Only uses the built-in BioPython pairwise sequence aligner so that it is usable across platforms without having to install additional alignment tools.

Still a work in progress.

An example of the analysis and the plots generated is in the Jupyter notebook in the 'example data' folder.
Before running the notebook, unzip the contents of the 'Sanger data.zip' file contained in the 'example data\data' folder.

### Installation
Download the verify_sanger package from GitHub or GitLab and save the source code in a local directory on your computer (e.g., on a Windows PC, C:\Users\username\Documents\Python Scripts\verify_sanger).

Open a Command Prompt (Windows) or the Terminal applicaiton (Mac) and navigate to the directory where you saved the FlowcCalNIST source code. Then install the FlowGateNIST package in editable mode using the command, "pip install -e ."
```
C:\Users\username\Documents\Python Scripts\verify_sanger>pip install -e .
```
Note that the "." at the end is important, as it indicates that the package should be installed from the current directory.