#!/bin/bash

#Run LCL saturation analysis

bsub -n 16 Rscript saturation_analysis_LCL.R 106 1
bsub -n 16 Rscript saturation_analysis_LCL.R 100 100
bsub -n 16 Rscript saturation_analysis_LCL.R 90 100
bsub -n 16 Rscript saturation_analysis_LCL.R 80 100
bsub -n 16 Rscript saturation_analysis_LCL.R 70 100
bsub -n 16 Rscript saturation_analysis_LCL.R 60 100
bsub -n 16 Rscript saturation_analysis_LCL.R 50 100
bsub -n 16 Rscript saturation_analysis_LCL.R 40 100
bsub -n 16 Rscript saturation_analysis_LCL.R 30 100
bsub -n 16 Rscript saturation_analysis_LCL.R 20 100
bsub -n 16 Rscript saturation_analysis_LCL.R 10 100
