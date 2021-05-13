# FLINO
The immunoFLuorescence Image NOrmalization (FLINO) repository provides R-scripts and data to perform and evaluate image normalizations methods and workflows.

Reference the manuscript.

Reference other repositories (Thrive, and the MOHA method)

The FLINO repository contains approximately 203 MB of data that can be used to represent Virtual Slides and become a ground truth dataset to evaluate alternative normalization methods and workflows.

Dependencies: The R-script files use the following four library packages:
library("png")
library("stringr")
library("plyr")   https://www.rdocumentation.org/packages/plyr/versions/1.8.6
library(fCI)      https://bioconductor.org/packages/release/bioc/html/fCI.html
library(NOISeq)   https://bioconductor.org/packages/release/bioc/html/NOISeq.html
library(qsmooth)  https://bioconductor.org/packages/release/bioc/html/qsmooth.html

Start R session
Within the R console
> setwd(<path to FLINO-main directory>)

Running the FLINO Evaluator R-script.

This example may require one to two minutes of computational time to complete.
> analysisRun = "eRuns_Grid256_Q75NZ_14VS.txt"
> source("Rcode/Run_FLINO_Evaluator.R")
> 
The output will be saved as
Results/results_eRuns_Grid256_Q75NZ_14VS.txt

Describe the output

This next example may require five to seven minutes of computational time to complete.
> analysisRun = "eRuns_NucleiSCA_Q50NZ_14VS.txt"
> source("Rcode/Run_FLINO_Evaluator.R")

The output will be saved as
Results/results_eRuns_NucleiSCA_Q50NZ_14VS.txt

Describe the output

Performing multiple evaluation in a single run.
This next example may require five to seven minutes of computational time to complete.
> analysisRun = "eRuns_Grid256_14VS.txt"
> source("Rcode/Run_FLINO_Evaluator.R")

The output will be saved as
Results/results_eRuns_Grid256_14VS.txt

Describe the output that compares TMM, MRN, etc.


Running all of the FLINO study evaluation runs.  This will take a very long time.
> analysisRun = " eRuns_FLINO.txt"

Would not advise running them in series.

Data/results_eRuns_FLINO.txt

Describe the output


Generate plots and figures presenting FLINO results on a CRC dataset. This may require one to two minutes of computational time to complete.
> source("Rcode/Gen_FLINO_Figures.R")

The output will be saved to the “figures” subdirectory 

Describe the output
Fig_4.tif
