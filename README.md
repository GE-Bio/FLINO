# FLINO
Copyright (c) 2021 General Electric Company

The immuno**FL**uorescence **I**mage **NO**rmalization (FLINO) repository provides R-scripts and data to perform and evaluate image normalizations methods and workflows.
The results of the FLINO study have been described in a manuscript that is currently under review and will be cited here once published. This repository contains the data and R-scripts used to generate the results of that analysis. This research was supported by the National Cancer Institute of the National Institutes of Health under award number R01CA208179.

The importance of the FLINO work has relevance to many types of downstream analysis including performing [multi-omics-heterogeneity-analysis (MOHA)](https://github.com/thrive-itcr/multi-omics-heterogeneity-analysis). 

The FLINO repository contains approximately 203 MB of grid object and segmented cell object data that can be used to represent Virtual Slides and become a ground truth dataset to evaluate alternative normalization methods and workflows. The data consists of 14 rounds of DAPI staining and imaging of the same physical tissue microarray (TMA) slide with 85 samples to represent a ground truth. We abstracted the individual rounds of DAPI staining and imaging to represent virtual slides. Each virtual TMA slide is the exact same 85 physical samples that have undergone a set of experimental conditions that introduce both random variation and systematic offsets between the virtual slides.

## Getting Started

The R-script files provided are dependent upon a number of R packages. These R libraries must be installed prior to running the FLINO R-scripts.

#### Dependencies: 
1. [fCI](https://bioconductor.org/packages/release/bioc/html/fCI.html)
2. [NOISeq](https://bioconductor.org/packages/release/bioc/html/NOISeq.html)
3. [qsmooth](https://bioconductor.org/packages/release/bioc/html/qsmooth.html)
4. [stringr](https://cran.r-project.org/web/packages/stringr/index.html)
5. [plyr](https://www.rdocumentation.org/packages/plyr/versions/1.8.6)
6. [png](https://cran.r-project.org/web/packages/png/index.html)

Download the content of the FLINO repository to your workstation, start a R sesssion, install the required R library dependencies, and then within the R console set the working directory to the location of the FLINO-main directory.
> setwd(C:\Users\...\FLINO-main)

#### Running the Run_FLINO_Evaluator.R script.

This first example may require one to two minutes of computational time to complete. First set the analysisRun parameterby typing the following into the R console. This is providing the file name of the evaluation run input. The evaluation run input files are stored within the **\FLINO-main\eRuns** directory. After setting the analysisRun, then simply run the Run_FLINO_Evaluator.R script:

> analysisRun = "eRuns_Grid256_Q75NZ_14VS.txt"

> source("Rcode/Run_FLINO_Evaluator.R")

Upon completion of the script the output will be saved in a results directory as: **\FLINO-main\Results\results_eRuns_Grid256_Q75NZ_14VS.txt**

Describe the output

##### This next example may require five to seven minutes of computational time to complete.
> analysisRun = "eRuns_NucleiSCA_Q50NZ_14VS.txt"

> source("Rcode/Run_FLINO_Evaluator.R")

The output will be saved as
Results/results_eRuns_NucleiSCA_Q50NZ_14VS.txt

Describe the output

##### Performing multiple evaluations in a single run.
This next example may require five to seven minutes of computational time to complete. The input file **\FLINO-main\eRuns\eRuns_Grid256_14VS.txt** is a tab delimited file with the first row representing the column names and input parameters for the evaluation run. Each subsequent row represents one evaluation run. There are five evaluation runs. The fourth column name is called **PARM_NORM_METHOD** and is the only change for the five evaluation runs. This input parameter is the normalization method that is being applied to the evaluation data. The five normalization methods being applied are:
1. TMM - trimmed mean of the M-values (Robinson and Oshlack 2010, Tarazona 2011, 2015)
2. MRN - median ratio normalization (Maza 2013)
3. Q75NZ - 75% quantile normalization excluding zero intensity objects
4. Q50NZ - 50% quantile normalization excluding zero intensity objects
5. MEDIAN - Median normalization



> analysisRun = "eRuns_Grid256_14VS.txt"

> source("Rcode/Run_FLINO_Evaluator.R")

The output will be saved as
Results/results_eRuns_Grid256_14VS.txt

Describe the output that compares TMM, MRN, etc.


##### Running all of the FLINO study evaluation runs.  This will take a very long time.
> analysisRun = " eRuns_FLINO.txt"

Would not advise running them in series.

Data/results_eRuns_FLINO.txt

Describe the output


#### Running the Gen_FLINO_Figures.R script.

The Gen_FLINO_Figures.R script will generate plots and figures presenting the FLINO study results. This may require one to two minutes of computational time to complete.
> source("Rcode/Gen_FLINO_Figures.R")

The output of this R script will be saved to the **\FLINO-main\figures** directory. The data presented in **\FLINO-main\figures\Fig_4.tif** is the application of grid-object normalization to BAX staining of three physical TMA slides with 85 samples that include four cell lines.

