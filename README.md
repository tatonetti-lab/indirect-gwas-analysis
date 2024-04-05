# Analysis for the Indirect GWAS project

Michael Zietz, Undina Gisladottir, Kathleen LaRow Brown, and Nicholas Tatonetti

Indirect GWAS is a method to produce GWAS summary statistics (coefficients, standard errors, p-values) for a linear combination of phenotypes using only other summary statistics.
This repository holds workflow files and notebooks used for the manuscript.
This analysis should be fully-reproducible and uses Snakemake.


## Outline

This project has the following general structure:

1. Validate that Indirect GWAS works and that it works for a variety of GWAS methods.
2. Demonstrate that it can be used to accelerate pan-biobank GWAS.
3. Demonstrate that it can be used to extrapolate to new phenotypes.
4. Evaluate the method's performance when phenotypes are not fit perfectly, and compare performance across categories of phenotypes.
