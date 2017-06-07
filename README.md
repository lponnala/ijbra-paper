
#### Description

This repository contains supplementary material and software code to accompany the following publication:

**Lalit Ponnala**. *Detecting slow-translating regions in E.coli*. International Journal of Bioinformatics Research and Applications 2010; 6(5): 522-530

<hr>

#### Instructions

The code has been tested using MATLAB version 7.6.0 (R2008a) on a Windows XP machine. 

- mainscript.m : this is the main script that detects clusters in Ecoli genes
- kscanstat.m : this function performs the actual clustering procedure
- mergeclust.m : this function merges the detected clusters to find the distinct ones
- evalclust.m : this script evaluates various features of the distinct clusters


You need not run mainscript.m, since all of its outputs have been provided as data files. To verify the results presented in the paper, simply run the script evalclust.m
