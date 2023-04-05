# CH-ACWE
This is an implementation of active contours without edges (ACWE) on solar extreme ultraviolet (EUV) images from the Atmospheric Imager Assembly aboard the Solar Dynamics Observatory (SDO-AIA). This implementation includes the ability to generate Confidence Maps based on the homogeneity of each region with respect to the core of the observed region.

Requirements: [environment.yml](environment.yml)  
This environment file specifies the packages necessary to implement all code contained within this repository, including packages necessary for downloading the dataset, generating segmentations (including confidence maps), and analyzing segmentation results.

## Downloading the Dataset
The dataset can be downloaded, and new datasets can be created, using the the tools provided in the `DatasetTools` folder.

- `DownloadLists`: This folder contains an organized lists of the dataset. These lists are organized into four `.csv` files, one for each Carrington Rotation (CR).
- `Carrington Rotation Start Dates.csv`: This file is a list of the start dates for each Carrington Rotation from CR -10 through CR 2300. This file is used by `DownloadByRotation.py` for both downloading and organizing the dataset.
- `DataManagmentTools.py`: Tools/functions for formatting dates to allow for request of data from [jsoc.stanford.edu](jsoc.stanford.edu).
- `DownloadByRotation.py`: Download `aia.lev1_euv_12s` and `hmi.M_720s` images at a 1 hour cadence for the specified Carrington rotation(s). 
  - This script will attempt to omit any time frame wherein at least one file is missing or does not and have a 'QUALITY' key of '0'. 
  - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
  - User will need to register their email at [http://jsoc.stanford.edu/ajax/register_email.html](http://jsoc.stanford.edu/ajax/register_email.html) and add that email to the appropriate variable in the `Key Variables` cell.
  - This script can be used to speed up the process of rebuilding the dataset. This is achieved by
    - Creating a temporary subfolder within the `DatasetTools` folder
    - Ensuring that the variable `traceFolder` points to that temporary subfolder
    - Setting the remaining variables in the `Key Variables` cell to ensure the correct CR is downloaded and saved where the user wishes
    - Running `DownloadByRotation.py`
    - Deleting the temporary subfolder
    - Running `RebuildDataset.py` with `traceFolder = 'DownloadLists/'` to download any missing files
- `GapCheck.py`: Inform the user as to the largest hour gap between entries in the specified CR within the dataset.
- `RebuildDataset.py`: Find and download any file that is missing from the dataset folder.
  - This script will rebuild the dataset.
  - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
  - User will need to register their email at [http://jsoc.stanford.edu/ajax/register_email.html](http://jsoc.stanford.edu/ajax/register_email.html) and add that email to the appropriate variable in the `Key Variables` cell.

## Generating ACWE Segmentations
### General Tools
The folder `ACWE_python_spring_2023` contains functions for running ACWE and saving the results.

- `ACWE_python_v3`: This folder contains the original ACWE functions, updated to be operate on python version 3.0 or greater.
- `acweConfidenceMapTools_v3.py`: Tools/functions for combining a segmentation group (collection of segmentations from the same EUV observation) in order to generate a confidence map.
- `acweFunctions_v6.py`: Tools/functions for preprocessing an EUV image, generating an initial mask, and running ACWE for both single output/segmentation and for a confidence map. 
  - The function `run_acwe` performs all processing and returns the final segmentation and initial mask. 
  - The function `run_acwe_confidenceMap` performs all processing and returns the final confidence map as a series of segmentations and initial mask.
  - Additional functions are also provided to perform each step separately.
- `acweRestoreScale.py`: Tools/functions for resizing a segmentation to match the spatial resolution of the input image.
  - Upscale a confidence map using `upscaleConMap`
  - Upscale a single segmentation using `upscale` 
  - Both functions take in the ACWE header and the segmentation or confidence map and return the same segmentation or confidence map, upscaled to match the resolution of the original EUV image.
- `acweSaveSeg_v5.py`: Tools/functions for saving and opening segmentations. 
  - The function `saveSeg` takes in the header of the original EUV image, the final segmentation(s), and the list of ACWE parameters to generate an .npz file which saves the final segmentation with a header outlining the ACWE parameters and a copy of the header for the original EUV image. 
  - The function `openSeg` opens and returns the header of the original EUV image, as a dictionary, the header outlining the options used to generate the ACWE segmentation, organized as a dictionary, and the final ACWE segmentation(s).
  - Both functions work for both single segmentations and for confidence maps.
