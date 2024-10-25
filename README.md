# Wasserstein-regression-with-empirical-measures

This repository contains codes necessary to replicate **Zhou and Müller (2024+)**: “Wasserstein regression with empirical measures and density estimation for sparse data.” The `grem` and `lrem` functions in the `code` folder, short for global and local Regression with Empirical Measures (REM), are designed for analyzing distributional data from heterogeneous populations.

## Folder Structure

The folder structure of this repo is as follows:

| Folder      | Detail                                                                                                    |
|:------------|:----------------------------------------------------------------------------------------------------------|
| code        | R scripts for the proposed approach                                                                       |
| data        | Data files                                                                                                |
| sim         | R scripts for simulations                                                                                 |

## code

| Data file | Detail                                                   |
|:----------|:---------------------------------------------------------|
| grem.R    | Global Regression with Empirical Measures (REM)          |
| lrem.R    | Local Regression with Empirical Measures (REM)           |
| lcm.R     | Least common multiple for a vector of integers           |
| bwCV.R    | Bandwidth selection for local REM using cross-validation |
| kerFctn.R | Kernel function                                          |

## data

| Data file                | Detail                                         |
|:-------------------------|:-----------------------------------------------|
| iseg.RData, isel.RData   | Simulation results for Settings I, II, III, IV |
| isegb.RData, iselb.RData | Simulation results for binomial distribution   |

## sim

R scripts to replicate simulation results in subsection 5.2 of the main text and Section S.2 of the Supplementary Material.

| Data file  | Detail                                                |
|:-----------|:------------------------------------------------------|
| simg.R     | Simulation for global REM with Gaussian distributions |
| siml.R     | Simulation for local REM with Gaussian distributions  |
| simgb.R    | Simulation for global REM with binomial distributions |
| simlb.R    | Simulation for local REM with binomial distributions  |
| simVisua.R | R script to replicate Figures 1, 2, and 3             |

## Report Errors

To report errors, please contact <ydzhou@ucdavis.edu>. Comments and suggestions are welcome.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-imbensxu" class="csl-entry">

Zhou, Y. and Müller, H.G., 2024. Wasserstein regression with empirical measures and density estimation for sparse data. Biometrics, (just-accepted), pp.1-12.

</div>
