# SARS-CoV-2_analysis
Library Prep and Analysis files

This folder contains all of the files used for initial PCR through analysis of SARS-CoV-2 sequences.  

Files/folders include:  

## SARS-COV2_SOP_AmplificationandLibraryPrep.pdf -  
  - contains full protocol for amplification and library preparation of SARS-CoV-2 clinical samples. Primer sequences and necessary reagents are listed.

## pipeline - 
  - contains all necessary files and scripts for trimming, alignment and variant calling of SARS-CoV-2 sequences. 
  - config file is the only file that needs editing in order to run the pipeline
  - Usage is demonstrated at the top of the Snakefile (snakemake) file
  *(many consensus sequences are generated from this pipeline. gatk-timo_consensus files are used most frequently)

## timo -
  - Full read me file with proper usage contained within timo folder. Also available at https://github.com/GhedinLab/timo
    - In-house variant calling pipeline which can be used for finding both major and minor SNVs 
     - ConsensusCoverage.py - generate consensus sequences from timo output
     - Addaminogenes.py - adds Amino acid and reference information to timo output
