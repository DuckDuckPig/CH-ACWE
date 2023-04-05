# CH-ACWE
This is an implementation of active contours without edges (ACWE) on solar extreme ultraviolet (EUV) images from the Atmospheric Imager Assembly aboard the Solar Dynamics Observatory (SDO-AIA). This implementation includes the ability to generate Confidence Maps based on the homogeneity of each region with respect to the core of the observed region.

Requirements: [environment.yml](environment.yml)  
This environment file specifies the packages necessary to implement all code contained within this repository, including packages necessary for downloading the dataset, generating segmentations (including confidence maps), and analyzing segmentation results.

## Downloading the Dataset
The dataset can be downloaded, and new datasets can be created, using the the tools provided in the `DatasetTools` folder.

- `DownloadLists`: This folder contains an organized lists of the dataset. These lists are organized into four `.csv` files, one for each Carrington Rotation (CR).
- `Carrington Rotation Start Dates.csv`: This file is a list of the start dates for each Carrington Rotation from CR -10 through CR 2300. This file is used by `DownloadByRotation.py` for both downloading and organizing the dataset.
- `DataManagmentTools.py`: Tools/Functions for formatting dates to allow for request of data from [jsoc.stanford.edu](jsoc.stanford.edu).
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
