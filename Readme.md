# Genetic susceptibility to schizophrenia through neuroinflammatory pathways is associated with retinal thinning: Findings from the UK-Biobank

## Getting Started
This repository contains all the data and analysis code to reproduce the
manuscript Retina OCT Project. These instructions describe how to obtain a copy
of the project up and running on your local machine for reproducing the
analysis described in the manuscript. The repository contains a Makefile
which reflects the dependencies of the analysis; analysis, figures and
manuscript can be produced by simply typing 'make' from the Unix command
line.


### Prerequisites
All analyses were conducted with the R software R version 3.4.4
(2018-03-15). Mixed models were estimated using the lme4 library, Python
2.7.14 and pysurfer (0.8.0) were also used. The full session info under
R can be found at the end of this file


## Installing
Clone the repository or download the zip file.

## Running the analysis
Change to the retpsy directory and run 'make analysis'.

## Producing the figures
Change to the retpsy directory and run 'make figures'. The figures can then
be found in output/figures.

## Producing the manuscript
Change to the retpsy directory and run 'make manuscript'. The
manuscript will be in src/.
