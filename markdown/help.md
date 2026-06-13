# 1. Instructions for using Web-server

## Input
- Input to the AquaLNCRpred  web-server consists of fasta formatted sequences. 
- These sequences can either be pasted directly on the form or uploaded as a file.
- Once you submit the button, please wait for a while before the results are available.

## Output
- As soon as the necessary calculations are completed, the prediction results are displayed on the page.
- If a user has submitted large number of sequences, the user may have to wait a bit before the results are available.
- When the predictions are done, the results can also be downloaded as a CSV file.
- The output consists of a dataframe of **FOUR** columns described as below:

### Sample Output

<img src="figures/sample_output.png" alt="drawing"/>

### Explaination of the output

```sh 
      Column 1: The accession number (ID) of the sequence that was predicted as a LncRNA or mRNA sequence.
      Columns 2-3: Probability scores.
      Column 4: The class of the predicted sequence, i.e. class with the highest probability i.e. LncRNA or mRNA
```

# 2. Datasets
- CV Dataset: [Download](CV_data.fasta)
- Independent Validation Set: [Download](IDS1.fasta)


# 4. Reference
  - ** AquaLNCRpred: A sequence-based framework for lncRNA prediction in *Eriocheir sinensis* (Chinese Mitten Crab) and related aquaculture species (Manuscript Submitted).**
