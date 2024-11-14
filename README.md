# Deep learning models of dehydroamino acid reductases from 4-track mRNA display profiling data
 
## About

This repository holds the code and a demo to reproduce the results of the LanJ<sub>C</sub> mRNA display profiling study (reference tbd). This a fork of the previously published [LazBF and LazDEF](https://github.com/avngrdv/mRNA-display-deep-learning) repository. The repository holds python scripts to reproduce the entire pipeline: from loading NGS output files to training and evaluating tensorflow-based models.

1. All metaparameters which include DNA and peptide library designs, a genetic code table, model hyperparameters, etc are specified in ```./srs/config.py```
2. Primary code is in ```./srs/clibas```
3. Model architecture is defined in ```./srs/tf/cnn_model.py``` (```cnn_vm_J2``` was used to build both the denoiser and LanJ<sub>C</sub> models).
5. Fully trained model weights stored as .h5 files are in ```./tf_trained_models```
6. The associated NGS sequencing data is uploaded to DDBJ (accession number: [DRA018924](https://ddbj.nig.ac.jp/search/entry/sra-submission/DRA018924))
\
\
\
For further details, refer to the accompanying paper: [_J. Am. Chem. Soc._ **2024**, 146, 31124−31136](https://pubs.acs.org/doi/full/10.1021/jacs.4c11013). Please cite it if you use this code.

## Dependencies

The scripts were written and tested for 

python 3.8.5 \
numpy 1.19.5 \
pandas 1.2.4 \
rdkit 2021.03.3 \
tensorflow 2.4.1 \
tensorflow-io 0.17.1 \
h5py 2.10.0 \
matplotlib 3.3.4 \
seaborn 0.11.1

## Usage examples

The usage of this code to build the pipeline from the manuscript is described in  ```./srs/4-track-pipeline.ipynb```

## License

_The code is released under the GNU General Public License v3.0 without any explicit or implied warranty. Use it at your own discretion._
