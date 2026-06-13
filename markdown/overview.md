# About AquaLNCRpred

- AquaLNCRpred is an XGB-based easy to use webserver for the prediction of LncRNA sequences. 
- AquaLNCRpred integrates a diverse set of sequence-based feature encodings, including k riboNucleotide composition (kNUComposition; (mono-, di-, and tri-nucleotide composition), codon Usage (CodonU), GC content (GC-c), maximum open reading frame length (mORFlen), Z curve (Zcurve), adaptive skip di-ribonucleotide composition (ASDC), di-riboNucleotide autocorrelation-autocovariance (AutoCorDiNUC), amphiphilic pseudo-k riboNucleotide composition-di(series) (ApkNUCdi), and their hybrid combination. 
- The AquaLNCRpred working module is based on optimal hybrid features that predicts whether a given sequence is a LncRNA or mRNA.


## Reference

- ** AquaLNCRpred: A sequence-based framework for lncRNA prediction in *Eriocheir sinensis* (Chinese Mitten Crab) and related aquaculture species (Manuscript Submitted).**

## AquaLNCRpred Algorithm

- Generation of training and independent datasets.
- Extraction of various feature encodings from the primary RNA sequences as described above. 
- 10-fold cross validation using different ML-based classifiers (e.g. KNN, NB, RF, SVM, &amp; XGB). 
- Optimization of hyperparameters. 
- The performance of the selected models are evaluated on the independent datasets separately. 
- Finally, the target sequence is predicted to be as a LncRNA or mRNA.
 
## Overview of AquaLNCRpred methodology

<img src="figures/Overview.jpg" alt="drawing" width="800" height="800"/>
