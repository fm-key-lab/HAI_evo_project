# Protein prediction and visualization

Before running the analysis, please check that you have: 
- access to the [AlphaFold 3 server](https://alphafoldserver.com/welcome) or a local installation
- installed [ABPS](#abps-installation)
- installed [ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html)


## ABPS installation

To run ABPS, follow the steps below:
```
conda create -n pqr python=3.12
conda activate pqr
cd ~/library/apbs
## install propka (dependency of pdb2pqr)
wget https://github.com/jensengroup/propka/archive/refs/tags/v3.5.1.tar.gz
untar v3.5.1.tar.gz
cd propka-3.5.1/ ; pip install . ; cd ../
pip install pdb2pqr ## requires python >= 3.8
conda install -c conda-forge gemmi
## install apbs
wget https://github.com/Electrostatics/apbs/releases/download/v3.4.1/APBS-3.4.1.Linux.zip
unzip APBS-3.4.1.Linux.zip
```

Lastly, modify `psize.py`
    - open file `psize.py` ( APBS-3.4.1/share/apbs/tools/manip/psize.py )
    - change line 174 from `if nsmall <= 0:` to `if nsmall[i] <= 0:`


## Visualize FimZ

1. Run AlphaFold 3 either via a local installation or on the server using the provided json file (`AF3_input/fimZ.json`).
    - NOTE: If using a local installation, one can implement this into the `run_AF3_MPCDF.sh` file
2. Update the `parameters.inc` and `visualization/chimerax_load_and_snapshot_fimZ.cxc` for the respective paths
3. Run `bash run_AF3_MPCDF.sh` 
4. Visualize the protein using the `visualization/chimerax_load_and_snapshot_fimZ.cxc` by running `/Applications/ChimeraX-1.6.1.app/Contents/bin/ChimeraX --script chimerax_load_and_snapshot_fimZ.cxc` within the respective directory