# FLINO
Copyright (c) 2021 General Electric Company

The immuno**FL**uorescence **I**mage **NO**rmalization (FLINO) repository provides R-scripts and data to perform and evaluate image normalizations methods and workflows.
The results of the FLINO study have been described in a manuscript that is currently under review and will be cited here once published. This repository contains the data and R-scripts used to generate the results of that analysis. This research was supported by the National Cancer Institute of the National Institutes of Health under award number R01CA208179.

The importance of the FLINO work has relevance to many types of downstream analysis including performing [multi-omics-heterogeneity-analysis (MOHA)](https://github.com/thrive-itcr/multi-omics-heterogeneity-analysis). 

The FLINO repository contains approximately 203 MB of grid object and segmented cell object data that can be used to represent Virtual Slides and become a ground truth dataset to evaluate alternative normalization methods and workflows. The data consists of 14 rounds of DAPI staining and imaging of the same physical tissue microarray (TMA) slide with 85 samples to represent a ground truth. We abstracted the individual rounds of DAPI staining and imaging to represent virtual slides. Each virtual TMA slide is the exact same 85 physical samples that have undergone a set of experimental conditions that introduce both random variation and systematic offsets between the virtual slides.

## Getting Started

The R-script files provided are dependent upon a number of R libraries. The R packages that must be installed prior to running the FLINO R-scripts include:

#### Dependencies: 
1. [fCI](https://bioconductor.org/packages/release/bioc/html/fCI.html)
2. [NOISeq](https://bioconductor.org/packages/release/bioc/html/NOISeq.html)
3. [qsmooth](https://bioconductor.org/packages/release/bioc/html/qsmooth.html)
4. [stringr](https://cran.r-project.org/web/packages/stringr/index.html)
5. [plyr](https://www.rdocumentation.org/packages/plyr/versions/1.8.6)
6. [png](https://cran.r-project.org/web/packages/png/index.html)

Download the content of the FLINO repository to your workstation, start a R sesssion, install the required R library dependencies, and then within the R console set the working directory to the location of the FLINO-main directory.
> setwd(C:\Users\...\FLINO-main)

#### Example 1: Running the Run_FLINO_Evaluator.R script.

This first example may require one to two minutes of computational time to complete. First set the analysisRun parameter to the file name of the evaluation run input. The evaluation run input files are stored within the **\FLINO-main\eRuns** directory. After the analysisRun parameter is defined, the Run_FLINO_Evaluator.R script is run using the source command:

> analysisRun = "eRuns_Grid256_Q75NZ_14VS.txt"

> source("Rcode/Run_FLINO_Evaluator.R")

Upon completion of the R script, the evaluation run output will be saved as a tab delimited text file in the results directory as: **\FLINO-main\Results\results_eRuns_Grid256_Q75NZ_14VS.txt**. The first row contains the column names of all the input parameters as well as the output values. The subsequent rows contain the output values for each evaluation. For this example there is only one evaluation. 

Several key input parameters for this one evaluation that can be found in the output file include:

* PARM_NORM_METHOD	with the value of **Q75NZ** for 75% quantile normalization excluding zero intensity objects.
* PARM_SEG_OBJ_NAME	with the value of **Grid256** for using grid objects of size dimension of 256 pixels for normalization.
* PARM_EVAL_SEG_OBJ_NAME with the value of **NucleiSCA** for using the Nuclei Segmented Cell objects for evaluating the performance of the normalization.
* PARM_EVAL_DAPI_RNDS with the value of that lists the 14 virutal slides (i.e. 14 rounds of DAPI staining) over which the normalization is being evaluated.

The two key output values for the evaluation include:

* EvalSegObj_SegObjErrCV_RawLogSpace	with the value of 0.673 representing the slide-to-slide batch effect error for the uncorrected images. This error is quantified as the **M**ean of all individual **E**valuation **O**bject - **C**oefficient of **V**ariations (MEO-CV) across the virtual slides. 
* EvalSegObj_SegObjErrCV_NormLogSpace	with the value of 0.104 representing the error in the evaluation segemented cell objects (i.e. **NucleiSCA**) after normalizing the images using grid objects (**Grid256**) and the 75% quantile normalization (**Q75NZ**) method.

#### Example 2: Running the Run_FLINO_Evaluator.R script.
This second example may require five to seven minutes of computational time to complete. The difference between this  example and prior one is that it uses the Nuclei Segmented Cell objects for both normalization and evaluation. The commands are:

> analysisRun = "eRuns_NucleiSCA_Q50NZ_14VS.txt"

> source("Rcode/Run_FLINO_Evaluator.R")

The output from example 2 will be saved as **\FLINO-main\Results\results_eRuns_NucleiSCA_Q50NZ_14VS.txt**. The key output values for this example 2 evaluation which can be compared to the above example include:

* PARM_NORM_METHOD	with the value of **Q50NZ** for 50% quantile normalization excluding zero intensity objects.
* PARM_SEG_OBJ_NAME	with the value of **NucleiSCA** for using the Nuclei Segmented Cell objects for normalization.
* PARM_EVAL_SEG_OBJ_NAME with the value of **NucleiSCA** for using the Nuclei Segmented Cell objects for evaluating the performance of the normalization.
* EvalSegObj_SegObjErrCV_RawLogSpace	with the value of 0.673 representing the slide-to-slide batch effect error for the uncorrected images. This error is quantified as the **M**ean of all individual **E**valuation **O**bject - **C**oefficient of **V**ariations (MEO-CV) across the virtual slides. 
* EvalSegObj_SegObjErrCV_NormLogSpace	with the value of 0.0996 representing the error in the evaluation segemented cell objects (i.e. **NucleiSCA**) after normalizing the images using segemented cell objects (**NucleiSCA**) and the 50% quantile normalization (**Q50NZ**) method.

#### Example 3: Performing multiple evaluations.
This next example may require five to seven minutes of computational time to complete. The input file **\FLINO-main\eRuns\eRuns_Grid256_14VS.txt** is a tab delimited file with the first row representing the column names and input parameters for the evaluation run. Each subsequent row represents one evaluation run. There are five evaluation runs. The fourth column name is called **PARM_NORM_METHOD** and is the only change for the five evaluation runs. This input parameter is the normalization method that is being applied to grid objects which is then used to correct the images. The performance of the method is then quantified using the segmented cell evaluation objects. The five normalization methods listed in the input file are:
1. TMM - trimmed mean of the M-values (Robinson and Oshlack 2010, Tarazona 2011, 2015)
2. MRN - median ratio normalization (Maza 2013)
3. Q75NZ - 75% quantile normalization excluding zero intensity objects
4. Q50NZ - 50% quantile normalization excluding zero intensity objects
5. MEDIAN - Median normalization

The commands for example 3 are:

> analysisRun = "eRuns_Grid256_14VS.txt"

> source("Rcode/Run_FLINO_Evaluator.R")

The output from example 3 will be saved as **\FLINO-main\Results\results_eRuns_Grid256_14VS.txt**. The output file will include the performance of the five normalization methods at correcting the slide-to-slide batch effect across 14 virtual slides. The performance for each normalization method is quantified by the DAPI segmented nuclei objects computed error (MEO-CV) and can be found in the output file under the two respective columns **PARM_NORM_METHOD** and **EvalSegObj_SegObjErrCV_NormLogSpace**. 

PARM_NORM_METHOD | EvalSegObj_SegObjErrCV_NormLogSpace
---------------- | -----------------------------------
TMM | 0.1067
MRN | 0.1079
Q75NZ | 0.1041
Q50NZ | 0.1067
MEDIAN | 0.1586


It is possible to run all of the FLINO study evaluations.  This would however take a vary long time if run within a single core.

> analysisRun = " eRuns_FLINO.txt"

The output for all evaluation runs can be found in the merged output file **\FLINO-main\Data\results_eRuns_FLINO.txt**. This output file is used by the next example.


#### Running the Gen_FLINO_Figures.R script.

The Gen_FLINO_Figures.R script will generate plots and figures presenting the FLINO study results. This may require one to two minutes of computational time to complete.
> source("Rcode/Gen_FLINO_Figures.R")

The output of this R script will be saved to the **\FLINO-main\figures** directory. The data presented in **\FLINO-main\figures\Fig_4.tif** is the application of grid-object normalization to BAX staining of three physical TMA slides with 85 samples that include four cell lines.

