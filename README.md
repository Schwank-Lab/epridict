# ePRIDICT: epigenetic PRIme editing preDICTion

## Overview

This repository is part of our study in Nature Biotechnology ["Machine learning prediction of prime editing efficiency across diverse chromatin contexts"](https://rdcu.be/dLu0f) and the initial [BioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.10.09.561414v1).

Predict prime editing efficiency based on chromatin context of a genomic location in K562 (ATCC CCL-243) cells.
Repository containing `Python` package for running trained `ePRIDICT` (epigenetic PRIme editing preDICTion) models. 
Models were trained in K562 cells. Prediction performance may vary in different cellular contexts. Check out our [publication](https://rdcu.be/dLu0f) for further details.

## Complementary Model

- **PRIDICT2.0**: This model focuses on the sequence-context based prime efficiency prediction and pegRNA design. We recommend to first select the most suitable pegRNA with PRIDICT2.0 and then assess its overall endogenous targetability with ePRIDICT. [Access PRIDICT2.0 GitHub Repository](https://github.com/uzh-dqbm-cmi/PRIDICT2)

## Resources

- **Supplementary Files**: [Access Here](https://github.com/Schwank-Lab/epridict/tree/supplementary_files)
- **Web Application**: For an online version of ePRIDICT, visit our [webapp](https://pridict.it/epridict)*.

*Default model for this repository and online webapp is `ePRIDICT-light`. For running the full `ePRIDICT` model, check the description below.

## Contact

For questions or suggestions, please either:
- Email us at [nicolas.mathis@pharma.uzh.ch](mailto:nicolas.mathis@pharma.uzh.ch)
- Open a GitHub issue

## Citation

If find our work useful for your research please cite:
- [Mathis et al., Nature Biotechnology, 2024](https://rdcu.be/dLu0f) (ePRIDICT and PRIDICT2.0):
```bibtex
@article{mathis2024machine,
  title={Machine learning prediction of prime editing efficiency across diverse chromatin contexts},
  author={Mathis, Nicolas and Allam, Ahmed and Tálas, András, and Kissling, Lucas and Benvenuto, Elena and Schmidheini, Lukas and Schep, Ruben and Damodharan, Tanav and Balázs, Zsolt and Janjuha, Sharan and others},
  journal={Nature Biotechnology},
  year={2024},
  publisher={Nature Publishing Group},
  doi={10.1038/s41587-024-02268-2}
}
```



## Getting Started

### Installation using Anaconda (Linux and Mac OS) 🐍
📣 `ePRIDICT` can only be installed on `Linux` and `Mac OS` since `pybigwig` package is not available for `Windows` 📣

The easiest way to install and manage Python packages on various OS platforms is through [Anaconda](https://docs.anaconda.com/anaconda/install/). Once installed, any package (even if not available on Anaconda channel) could be installed using pip. 

* Install [Anaconda](https://docs.anaconda.com/anaconda/install/).
* Start a terminal and run:
    ```shell
    # clone ePRIDICT repository
    git clone https://github.com/Schwank-Lab/epridict.git
    # navigate into repository
    cd epridict
    # create conda environment and install dependencies for ePRIDICT
    # (only has to be done before first run/install)
    conda env create -f epridict_env.yml

    # activate the created environment
    conda activate epridict
    ```

* Next, downloading ENCODE datasets is needed for prediction with ePRIDICT. Files will be downloaded in `bigwig` folder.
  Note: For running the `full` model, 455 datasets will be downloaded, requiring **624 GB** of storage space!
        For the `light` model, with slightly lower accuracy, 6 datasets will be downloaded, requiring **5.3 GB** of storage space.
  ```shell
  # make download script executable:
  chmod +x epridict_download_encode.sh
  # run download script:
  ./epridict_download_encode.sh light # or ./epridict_download_encode.sh full
  ```
  
--------------------------

### Running ePRIDICT in 'single' mode:
  #### Required:
  -  `--chromosome`: Chromosome of desired location. Format: "chr1", "chr2", ... (human chromosomes; Y-chromosome not supported)
  -  `--position_hg38`: Position within chromosome (hg38). Example: "1192940"
  #### Optional:
  -  `--use_full_model`: Use `full` model (455 ENCODE datasets) for prediction. Only possible when downloaded all datasets with `./epridict_download_encode.sh full`. Default: `light` model

```shell

python epridict_prediction.py single --chromosome chr3 --position_hg38 44843504
# for full model:
# python epridict_prediction.py single --chromosome chr3 --position_hg38 44843504 --use_full_model
```

--------------------------

### Running in 'batch' mode:
  ####  Required:
  -  `input_filename`: Input file name - name of .csv file that has two columns [`chromosome`, `position_hg38`]. See `sample_epridict_batch.csv` in the `./input` folder.
  #### Optional:
  -  `--output-fname`: Alternative output filename. Default is `input_filename_output.csv`. (e.g. in the example below the default output will be `sample_epridict_batch_output.csv`)
  -  `--use_full_model`: Use `full` model (455 ENCODE datasets) for prediction. Only possible when downloaded all datasets with `./epridict_download_encode.sh full`. Default: `light` model
  
```shell
# python epridict_prediction.py batch input_filename
python epridict_prediction.py batch sample_epridict_batch.csv
# for full model and alternative output filename:
# python epridict_prediction.py batch sample_epridict_batch.csv --output-fname alternative_output_filename_batch.csv --use_full_model
```
--------------------------
