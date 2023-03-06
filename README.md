# CBPT-GLME-Code
Last Updated by Seth Konig 1/17/2023

## Description

This repository contains example code and data that can be used to replicate a subset of results found in the follow publication (TBD):

"Flexible Multi-Step Hypothesis Testing of Human ECoG Data using Cluster-based Permutation Tests with GLMEs".


## Matlab Code

* Main Functions
  * **clusterBasedGlmeERD.m:** main cluster based permutation test (CBPT) code.
      * Uses both LMEs and GLMEs depending on input parameters.
* Supporting Functions
  * findgaps.m: function that finds gaps in vector indices.
* Example Code
  * exampleFaceHouseAnlaysisBB.m: example code for analyzing broadband (BB) power data.
    * Includes example code for group-level analysis.
  * exampleFaceHouseAnlaysisBurst.m: example code for analyzing high-gamma burst data.
  * exampleFaceHouseAnlaysisERP.m: example code for analyzing ERP data.




## Preprocessed Data
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
