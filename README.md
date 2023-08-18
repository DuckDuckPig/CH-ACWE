# CH-ACWE
This is an implementation of active contours without edges (ACWE) on solar extreme ultraviolet (EUV) images from the Atmospheric Imager Assembly aboard the Solar Dynamics Observatory (SDO-AIA). This implementation includes the ability to generate Confidence Maps based on the homogeneity of each region with respect to the core of the observed region. This code is related the paper [https://doi.org/10.48550/arXiv.2308.05679](https://doi.org/10.48550/arXiv.2308.05679) (submited to Solar Physics).

## Requirements: [environment.yml](environment.yml)
This environment file specifies the packages necessary to implement all code contained within this repository, including packages necessary for downloading the dataset, generating segmentations (including confidence maps), and analyzing segmentation results.

## Downloading the Dataset
The dataset can be downloaded, and new datasets can be created, using the the tools provided in the `DatasetTools` folder.

- `DownloadLists`: This folder contains an organized lists of the dataset. These lists are organized into four `.csv` files, one for each Carrington Rotation (CR).
- `Carrington Rotation Start Dates.csv`: This file is a list of the start dates for each Carrington Rotation from CR -10 through CR 2300. This file is used by `DownloadByRotation.py` for both downloading and organizing the dataset.
- `DataManagmentTools.py`: Tools/functions for formatting dates to allow for request of data from [jsoc.stanford.edu](jsoc.stanford.edu).
- `RebuildDataset.py`: Find and download any file that is missing from the dataset folder.
  - This script will rebuild the dataset direcly from the specified file present in the `DownloadLists` folder.
  - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
  - User will need to register their email at [http://jsoc.stanford.edu/ajax/register_email.html](http://jsoc.stanford.edu/ajax/register_email.html) and add that email to the appropriate variable in the `Key Variables` cell.
- `DownloadByRotation.py`: Download `aia.lev1_euv_12s` and `hmi.M_720s` images at a 1 hour cadence for the specified Carrington rotation(s). 
  - This script will attempt to omit any time frame wherein at least one file is missing or does not and have a `QUALITY` key of `0`. 
  - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
  - User will need to register their email at [http://jsoc.stanford.edu/ajax/register_email.html](http://jsoc.stanford.edu/ajax/register_email.html) and add that email to the appropriate variable in the `Key Variables` cell.
  - This script can be used to speed up the process of rebuilding the dataset. This is achieved by
    1. Creating a temporary subfolder within the `DatasetTools` folder
    2. Ensuring that the variable `traceFolder` points to that temporary subfolder
    3. Setting the remaining variables in the `Key Variables` cell to ensure the correct CR is downloaded and saved where the user wishes
    4. Running `DownloadByRotation.py`
    5. Deleting the temporary subfolder
    6. Running `RebuildDataset.py` with `traceFolder = 'DownloadLists/'` to download any missing files
- `GapCheck.py`: Inform the user as to the largest hour gap between entries in the specified CR within the dataset.

## Generating ACWE Segmentations

![ACWE Segmentation Sample](/Figures/ReadMe/aia.lev1_euv_12s.2010-09-18T000008Z.193.image_lev1.fits.png?raw=true "ACWE Segmentation Sample")

### General Tools
The folder `ACWE_python_spring_2023` contains functions for running ACWE and saving the results.

- `ACWE_python_v3`: This folder contains the original ACWE functions, updated to operate on python version 3.0 or greater.
- `acweConfidenceMapTools_v3.py`: Tools/functions for combining a segmentation group (collection of segmentations from the same EUV observation) in order to generate a confidence map.
- `acweFunctions_v6.py`: Tools/functions for preprocessing an EUV image, generating an initial mask, and running ACWE for both single output/segmentation and for a confidence map. 
  - The function `run_acwe` performs all processing and returns the final segmentation and initial mask. 
  - The function `run_acwe_confidenceMap` performs all processing and returns the final confidence map as a series of segmentations and initial mask.
  - Additional functions are also provided to perform each step separately.
  - These functions will work for both AIA and Solar Terrestrial RElations Observatory (STEREO) observations, however a resize parameter of 4 and seeding parameter `alpha` in the range \[0.8,0.9\] are recommended for STEREO data. 
- `acweRestoreScale.py`: Tools/functions for resizing a segmentation to match the spatial resolution of the input image.
  - Upscale a confidence map using `upscaleConMap`
  - Upscale a single segmentation using `upscale` 
  - Both functions take in the ACWE header and the segmentation or confidence map and return the same segmentation or confidence map, upscaled to match the resolution of the original EUV image.
- `acweSaveSeg_v5.py`: Tools/functions for saving and opening segmentations. 
  - The function `saveSeg` takes in the header of the original EUV image, the final segmentation(s), and the list of ACWE parameters. It generates an .npz file which saves the final segmentation with a header outlining the ACWE parameters and a copy of the header for the original EUV image. 
  - The function `openSeg` opens and returns the header of the original EUV image, as a dictionary, the header outlining the options used to generate the ACWE segmentation, organized as a dictionary, and the final ACWE segmentation(s).
  - Both functions work for both single segmentations and for confidence maps.

### Standard Segmentation
The script `runACWEdefault.py`, which generates the default implementation of ACWE on Solar EUV images generated by AIA is located in the folder `Standard`.

Note: 
- User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories and desired EUV wavelength (193 angstroms is the assumed default).
- The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately.

### Confidence Maps
The standard implementation of ACWE confidence maps is generated using the script `runACWEconfidenceLevelSet_Default.py`, within the `ConfidenceMapping` folder. 

Note: 
- User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories and desired EUV wavelength (193 angstroms is the assumed default).
- The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately.
- This script will generate all specified segmentations, regardless of whether or not a change of target will occur with the given parameters chosen in in the `Key Variables` cell. When change of target occurs, a valid confidence map can be extracted from the ensemble using the `smartConMap` function provided in `acweConfidenceMapTools_v3.py` (in the `ACWE_python_spring_2023` folder).

### Other Segmentations
- Segmentations generated at any spatial resolution other than 512x512 pixels should be performed using the script `runACWEscaledDefault.py` located within the `Scaled` folder.
  - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories and desired EUV wavelength (193 angstroms is the assumed default).
  - The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately.
  - Under the assumption that the input image is 4096x4096 pixels (the resolution of AIA), if the resize parameter variable `resize_param = 8`, this function will generate a standard segmentation.
- The script `runACWEconfidenceIndependent_Old.py` (inside of the `ConfidenceMapping`) can be used to generate confidence maps.
  - This script is not optimized, and operates by generating each map independently, starting from the initial seed each time. It is >2.9 times slower than the optimized script `runACWEconfidenceLevelSet_Default.py`.
  -  User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories and desired EUV wavelength (193 angstroms is the assumed default).
  - The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately.
  - This script will generate all specified segmentations, regardless of whether or not a change of target will occur with the given parameters chosen in in the `Key Variables` cell. When change of target occurs, a valid confidence map can be extracted from the ensemble using the `smartConMap` function provided in `acweConfidenceMapTools_v3.py` (in the `ACWE_python_spring_2023` folder).
  
## Analyzing ACWE Segmentations
Analysis of the stability and consistency of ACWE can be performed using the following tools.

### Spatial Resolution Effects
The tools in the folder `Scaled/Analysis/` can be used to to determine the effect that spatial decimation of the input image has on final segmentation.

- The script `analizeACWEscaledDefault.py` generates an `.npz` file which outlines the similarity of segmentation at the specified scales, compared to an ACWE segmentation generated from the EUV image at full scale. 
  - This script will return the intersection over union (Jacquard Index/IOU), structural similarity (SSIM), global consistency error (GCE), and local consistency error (LCE) for each segmentation at each spatial resolution as a function of interpolation method used to upscale the image to the original resolution.
  - Before running this script be sure to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
- The Jupyter Notebook `visulization_earlyData-BW.ipynb` generates a box and whisker plot for a block of CRs, from the data generated by `analizeACWEscaledDefault.py`.
- The Jupyter Notebook `visulization_singleCR-BW.ipynb` generates a box and whisker plot for the specific CR the user chooses in the `Key Variables` cell, from the data generated by `analizeACWEscaledDefault.py`
- The Jupyter Notebook `Scaling Samples.ipynb` will generate, display, and save a figure showing the effects of spatial decimation for the files specified in the second cell \(`In[2]`\)
  - The figures will be saved in a folder within the project space that the notebook creates
  - User will need to adjust the variables in the second cell (`In[2]`) to point to the correct directories.
- The script `TimeCheckScaled.py`, which is located in the `Scaled/` folder, can be used to determine the speed increase offered by reducing the spatial resolution of the input EUV images.
  - In order to run this script the user will need to un-comment, in `Scaled/runACWEscaledDefault.py` and `Standard/runACWEdefault.py`, all blocks with the header `Time` as well as the variables in the `Time ACWE` section of the `Key Variables` cell ***prior*** to running either script on the specified CR.
  - The variable `CR` in the `Key Variables` cell must match across all three aforementioned scripts.
  - The function assumes that the standard resolution segmentations were generated using `Standard/runACWEdefault.py` and will search the `Standard/` folder for data related to the time needed to generate the standard resolution segmentations.

### Effects Across Small Temporal Changes
The tools in the folder `TemporalEffects/Analysis/` can be used to determine the consistency of ACWE segmentation across short time spans where CH evolution is expected to be minimal. This can be used to determine how consistent ACWE is in identifying and segmenting CH regions.

- The script `analizeTempEffects.py` generates an `.npz` file which outlines the similarity of segmentation at the specified scale (the default 1/8th spatial resolution), compared to ACWE segmentation at the same scale generated from all EUV images in the same CR that are within +-12 hours. 
  - This script will return the IOU, SSIM, GCE, and LCE for each segmentation compared to the succeeding and preceding 12 hours of segmentations.
  - Before running this script be sure to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
- The Jupyter Notebook `visulization_earlyData.ipynb` generates a plot from the data created by `analizeTempEffects.py` for a block of CRs. 
- The Jupyter Notebook `visulization_earlyData.ipynb` generates a plot from the data created by `analizeTempEffects.py` for the specific CR the user chooses in the `Key Variables` cell.
- The Jupyter Notebook `visulization_IndividualExamples-SingleCR.ipynb` generates a plot showing the IOU into the future for each user-specified entry in the user-specified CR. Both variables are in the `Key Variables` cell. 

### Effects of Confidence Map Optimization 
The standard implementation of ACWE confidence maps optimizes ensemble generation by relying on the previous finalized segmentation as input for generating the next segmentation in the ensemble. In order to determine how this process effects final segmentation quality, the following tools are provided:

- The script `TimeCheck.py` can be used to determine the speed increase offered by the optimized confidence map process (`runACWEconfidenceLevelSet_Default.py`) over the non-optimized approach (`runACWEconfidenceIndependent_Old.py`).
  - In order to run this script the user will need to un-comment, in `runACWEconfidenceLevelSet_Default.py` and `runACWEconfidenceIndependent_Old.py`, all blocks with the header `Time` as well as the variables in the `Time ACWE` section of the `Key Variables` cell ***prior*** to running either script on the specified CR.
  - The variable `CR` in the `Key Variables` cell must match across all three aforementioned scripts.
- The script `analizeConmapSim.py`, in the folder `ConfidenceMapping/AnalysisInternal/`, performs three analyses:
  1. Determine, using weighted IOU, SSIM, GCE, and LCE how similar the confidence maps generated using the optimized approach are to confidence maps generated from an ensemble of independent segmentations.
  2. Determine, for the confidence maps generated from an ensemble of independent segmentations, the area of each segmentation that was not present in the previous segmentation, and the area that was not carried over to the succeeding segmentation.
  3. Determine, for the  optimized confidence map process, the area of each segmentation that was not present in the previous segmentation, and the area that was not carried over to the succeeding segmentation.
  - Before running this script be sure to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
- The jupyter notebook `ConMapStats_FullData.ipynb` (in the folder `ConfidenceMapping/AnalysisInternal/`), reports on all analyses performed in `analizeConmapSim.py`.
- The jupyter notebook `ConmapSamples.ipynb` (in the folder `ConfidenceMapping/AnalysisInternal/`) will generate, display, and save visualizations of the confidence maps that exist for the files specified in the second cell \(`In[2]`\).
  - The figures will be saved in a folder within the project space that the notebook creates.
  - User will need to adjust the variables in the second cell (`In[2]`) to point to the correct directories.

### Confidence Map Change of Target Analsis
The folder `ConfidenceMapping/AnalysisGrowthAndIntensity/` contains tools to evaluate the behaviour of the ensemble of segmentations used to generate confidence maps to help identify the characteristics that can be used to find and correct change of target cases.

- The script `analizeGrowthAndIntensity.py`compares segmentations to the input seed to allow the user to determine what differences exist between change of target cases and valid segmentations.
  - The script will report on the intensity of the input seed and each segmentation, providing min, mean, and max for each.
  - The script will report the area of the input seed, area of each segmentation, and the percentage of the initial seed that is retained in each segmentation.
  - Before running this script be sure to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
- The Jupyter Notebook `Growth Rate and mean intensity_SingleCR.ipynb` reports the results from `analizeGrowthAndIntensity.py` for the user-specified CR.
- The Jupyter Notebook `Change of Target Methology Check 5percent.ipynb` reports the list of change of target cases that were identified via the method implemented in the `smartConMap` function.
- The Jupyter Notebooks `Special Cases CR2099.ipynb` and `Special Cases CR2133.ipynb` provide figures for the specified cases, showing the effects of the `smartConMap` function. These examples include genuine cases of change of target, as well as cases without change of target to ensure the method operates properly.

### Validation of Confidence Maps 
The folder `ConfidenceMapping/AnalysisMagentogram/` contains tools used to determine the skew of the underlying magnetic field, as a function of confidence level, for the confidence maps generated via ACWE.

- The script `ConMapSkewness.py` calculates and returns the skew of the underlying magnetic field as a function of confidence for each CH identified in each segmentation.
  - Before running this script be sure to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
- The script `Skew Check SingleCR.py` converts the results from `ConMapSkewness.py` into a series of figures to aid the user in determining if CH regions correspond to regions of high unpopularity in the magnetic field.
  - Before running this script be sure to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
- the script `FigureMassCopy.py` takes the figures generated by `Skew Check SingleCR.py`, which are organized by record time, and groups them in a single folder, to help the user parse through all figures. 
  - Before running this script be sure to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.

### Intensity and Dynamic range Effects
The `Intensity/` folder and the `Intensity/Analysis/` subfolder contain tools used to determine the effects that compressing SDO-AIA observations by represetiing all intesnity using 265 discretized values (such as `.png`) have on the quality of ACWE segmentations.

- The script `runIntensityInvDefault.py` generates a series of segmentations. Each segmentation is generated from a copy of the EUV image that has been remapped to `uint8` (similar to popular image formats like `.png` and `jpeg`) via various remapping schemes.
  - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories and desired EUV wavelength (193 angstroms is the assumed default).
  - The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately.
  - The script relies on the functions in `pngScale2.py`
- The file `pngScale2.py` contains tools/functions to facilitate remapping of the EUV images to the 256 desecrate values.
- the script `Intensity Data.py` reports, for the specified CR, the range of intensities present in the original EUV images
  - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories and desired EUV wavelength (193 angstroms is the assumed default).
  - The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately.
- The script `analizeIntInvEffects.py`, located within the `Intensity/Analysis/` subfolder determines the similarity of segmentations generated from the remapped images compared to segmentations generated from the original EUV image.
  - This script will return the IOU, SSIM, GCE, and LCE for each comparison.
  - Before running this script be sure to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
- The Jupyter Notebook `visulization_earlyData-BW.ipynb` generates a box and whisker plot for a block of CRs, from the data generated by `analizeIntInvEffects.py`.
- The Jupyter Notebook `visulization_singleCR-BW.ipynb` generates a box and whisker plot for the specific CR the user chooses in the `Key Variables` cell, from the data generated by `analizeIntInvEffects.py`
- The Jupyter Notebook `Intensity Samples.ipynb` will generate, display, and save a figure showing the effects of intensity remapping and dynamic range decimation for the files specified in the second cell \(`In[2]`\)
  - The figures will be saved in a folder within the project space that the notebook creates.
  - User will need to adjust the variables in the second cell (`In[2]`) to point to the correct directories.
