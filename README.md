# FIGS
FIGS (Featured-ions Guided Stoichiometry) is a universal method for peptide quantification of data-independent acquisition mass spectrometry data. FIGS consists of two modules, namely spectral deconvolution and peptide quantification.

#### Parameters
FIGS_Deconvolute.py
- mzmlPath: The path to the mzML file, including the file name.
- libPath: The path to the spectral library file, including the file name.
- tol: The instrument mass tolerance in *p.p.m* (parts per million). Defaults to 10.
- Top10First: A *bool* parameter indicating whether to apply the top10-first strategy. Defaults to *False*.
- distance: The minimum swap distance for generate a decoy library. Defaults to 100.

FIGS_Quant.R
- dirPath: The directory path to *Coeffs.csv* exported from *FIGS_Deconvolute.py*.
- FillCoeffsFlag: A *bool* parameter indicating whether to improve the coefficient profile. Defaults to *FALSE*.
