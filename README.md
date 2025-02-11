# CBPT-GLME-Code
Last Updated by Seth Konig 2/11/2025

## Description

This repository contains example code and data that can be used to replicate a subset of results found in the follow publication (TBD):

"Flexible Multi-Step Hypothesis Testing of Human ECoG Data using Cluster-based Permutation Tests with GLMEs".


## Matlab Code

* Main Functions
  * **clusterBasedGlmeERD.m:** main cluster based permutation test (CBPT) code.
      * Uses both LMEs and GLMEs depending on input parameters.
 *  **clusterBasedGlmeERD2D.m:** cluster based permutation test (CBPT) code for 2D data (e.g. time-frequency) code.
     * Uses LMEs to anlayze normal data e.g. time-frequency power data.
* Supporting Functions
  * findgaps.m: function that finds gaps in vector indices.
* Example Code
  * exampleFaceHouseAnlaysisBB.m: example code for analyzing broadband (BB) power data.
    * Includes example code for group-level analysis.
  * exampleFaceHouseAnlaysisBurst.m: example code for analyzing high-gamma burst data.
  * exampleFaceHouseAnlaysisERP.m: example code for analyzing ERP data.


## 2D anlaysis Disclaimer
We provide the 2D cluster-based permutation LMEs code as a courtesy for those interested in using our proposed method for 2D (e.g. time-frequency) analysis. However, please note that this code has not been peer-reviewed though we have done our best to ensure its validity. 

Analyzing 2D data (e.g., time-frequency) with multiple fixed effects presents additional complexities, as clusters inherently span three dimensions (e.g., time × frequency × fixed effect). This makes it more challenging to properly control for family-wise error rate (FWER) without sacrificing statistical power. Additionally, defining clusters in 2D analyses is more complex than in 1D analyses. To address this, we offer multiple cluster-definition methods, but results may vary depending on the data's complexity. Therefore, we strongly recommend comparing cluster-based permutation results with false discovery rate (FDR) corrected results (e.g.  https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh).



## Preprocessed Data

Example preprocsesed data can be found here:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7703148.svg)](https://doi.org/10.5281/zenodo.7703148).

The data included in the **exampleGlmeFaceHouseAnalysisData.mat** file was derived from the following literature and data repository:

Miller, Kai Joshua. (2016). A library of human electrocorticographic data and analyses. Stanford Digital Repository. Available at: https://purl.stanford.edu/zk881ps0522

This data is also licensed under a Creative Commons Attribution Share Alike 4.0 International license (CC BY-SA).



## License

Shield: [![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg

&copy; Herman-Darrow Lab, University of Minnesota, 2023
