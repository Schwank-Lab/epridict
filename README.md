# ePRIDICT: epigenetic PRIme editing preDICTion

This repository is part of the study [Predicting prime editing efficiency across diverse edit types and chromatin contexts with machine learning](https://pridict.it/epridict).

Predict prime editing efficiency based on chromatin context of a genomic location in K562 cells.
Repository containing `Python` package for running trained `ePRIDICT` (epigenetic PRIme editing preDICTion) models. 
Models were trained in K562 cells. Prediction performance may vary in different cellular contexts. Check out the publication for further details.

To run `ePRIDICT` online, see our webapp [pridict.it/epridict](https://pridict.it/epridict).*

*Default model for this repository and online webapp is `ePRIDICT-light`. For running the full `ePRIDICT` model, check the description below.

--------------------------

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
    # create conda environment and install dependencies for ePRIDICT (only has to be done before first run/install)
    # use epridict_linux for linux machine or epridict_mac for a macbook
    conda env create -f epridict_linux.yml # epridict_mac.yml for macbook

    # activate the created environment
    conda activate epridict
    ```

* Next, downloading ENCODE datasets is needed for prediction with ePRIDICT. Files will be downloaded in `bigwig` folder.
  Note: For running the `full` model, 455 datasets will be downloaded, requiring **624 GB** of storage space!
        For the `light` model, with near on-par performance, 6 datasets will be downloaded, requiring **5.3 GB** of storage space.
  ```shell
  ./epridict_download_encode.sh light # or ./epridict_download_encode.sh full
  ```
  
--------------------------

### Running ePRIDICT in 'manual' mode:
  #### Required:
  -  `--chromosome`: Chromosome of desired location. Format: "chr1", "chr2", ...
  -  `--position_hg38`: Position within chromosome (hg38). Example: "1192940"
  #### Optional:
  -  `--use_full_model`: Use `full` model (455 ENCODE datasets) for prediction. Only possible when downloaded all datasets with `./epridict_download_encode.sh full`. Default: `light` model

```shell

python epridict_prediction.py manual --chromosome chr3 --position_hg38 44843504
# for full model:
# python epridict_prediction.py manual --chromosome chr3 --position_hg38 44843504 --use_full_model
```

--------------------------

### Running in 'batch' mode:
  ####  Required:
  -  `input_filename`: Input file name - name of .csv file that has two columns [`chromosome`, `position_hg38`]. See `sample_epridict_batch.csv` in the `./input` folder.
  #### Optional:
  -  `--output-fname`: Alternative output filename. Default is `input_filename_output.csv`.
  -  `--use_full_model`: Use `full` model (455 ENCODE datasets) for prediction. Only possible when downloaded all datasets with `./epridict_download_encode.sh full`. Default: `light` model
  
```shell

python epridict_prediction.py batch sample_epridict_batch.csv
# for full model and alternative output filename:
# python epridict_prediction.py batch sample_epridict_batch.csv --output-fname alternative_output_filename_batch.csv --use_full_model
```

### Citation
If you find our work is useful in your research, please cite the following paper: