#!/bin/bash

#Run Fib saturation analysis

bsub -n 16 Rscript saturation_analysis_fib.R 99 1
bsub -n 16 Rscript saturation_analysis_fib.R 90 100
bsub -n 16 Rscript saturation_analysis_fib.R 80 100
bsub -n 16 Rscript saturation_analysis_fib.R 70 100
bsub -n 16 Rscript saturation_analysis_fib.R 60 100
bsub -n 16 Rscript saturation_analysis_fib.R 50 100
bsub -n 16 Rscript saturation_analysis_fib.R 40 100
bsub -n 16 Rscript saturation_analysis_fib.R 30 100
bsub -n 16 Rscript saturation_analysis_fib.R 20 100
bsub -n 16 Rscript saturation_analysis_fib.R 10 100
