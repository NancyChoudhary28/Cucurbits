Link to Bonndata repository: https://doi.org/10.60507/FK2/ZNVSLA

This repository includes code used in the article:   
# [Out of the blue: Family-wide loss of anthocyanin biosynthesis in Cucurbitaceae](https://doi.org/10.1101/2025.10.06.680802)
https://doi.org/10.1101/2025.10.06.680802

### Getting datasets and identifying flavonoid and carotenoid biosynthesis genes
**1) Retrieving data and extracting coding sequences and longest isoforms**  
Post-downloading genome sequence (FASTA) and gene annotation (GFF3 format), **run_agat_wrapper.sh** was used to extract the coding sequences. Additionally, it also outputs single-isoform coding sequences.

   ```
   Input:
   GENOME_DIR         directory containing all genome fasta files
   GFF_DIR            directory containing all gff3 files

   Output:
   CDS_DIR            directory containing all coding sequence fasta files
   SINGLEISO_GFF_DIR  directory containing all gff3 files with the longest isoform kept
   SINGLEISO_CDS_DIR  directory containing all single-isoform coding sequence fasta files
   LOG_DIR            directory containing all log files
   ```

_*Note: The base names of the input genome fasta and gff3 files belonging to each other should match. The output files have the same base names in the output directories. Please change the input and output directory paths in the wrapper script before running it._

**2) Translating coding sequences into polypeptide sequence** 
**transeq.py** available [here](https://github.com/bpucker/PBBtools/tree/main/transeq) was used to convert all coding sequences to polypeptide sequences.

**3) Identifying contaminated datasets (both genomics and transcriptomics)** 
**run_kraken_braken_wrapper.sh** was used to classify coding sequences as contaminated or not. Datasets with more than 10% of coding sequences identified as Non-Cucurbitaceae were classified as contaminated and removed from subsequent analyses in the study.

   ```
   Input:
   input_dir          path to the directory containing all coding sequences (FASTA format)
   db_path            path to the database used for Kraken classification

   Output:
   output_base        path to output directory
   ```
   
_*Note: Please change the input and output directory paths in the wrapper script before running it._

**4) Identifying pigment biosynthesis genes in non-contaminated datasets** 
**KIPEs3.py** available [here](https://github.com/bpucker/KIPEs) was used to identify flavonoid and carotenoid biosynthesis pathway genes in all datasets' polypeptide sequences.

**5) Collecting candidate flavonoid biosynthesis genes for phylogenetic tree construction** 
**collect_KIPEs_sequence.py** was used to collect flavonoid biosynthesis candidate gene sequences from KIPEs results, and the sequences were used for phylogenetic tree construction.
  ```
   Usage
   python3 collect_KIPEs_sequence.py
   --directory        path to the KIPEs result folder
   --enzymes          list of genes of interest (comma-separated if multiple)
   --output           output file
   --threshold        minimum percentage of conserved residues (default:70.0)
   ``` 

**6) Preparing iTOL annotation file for visualizing ortholog clades** 
**get_conserved_residue_perc.py** was used to create an iTOL gradient annotation file to show the percent of conserved residues required for the said function. 
```
Usage
python3 get_conserved_residue_perc.py
-k               path to KIPEs result directory
-i               input file containing fasta headers (without '>')
-o               path to the output file
-g               candidate gene name to match in summary.txt (e.g. DFR)
```
The FASTA headers must be candidates of a single gene, and the name of this gene is specified using the `-g` flag. Once the output file is created using the above script, some manual edits are required to match the iTOL annotation file criteria. An example **_DFR_gradient.txt_** file is provided for reference. 


### Expression analyses
**1) Selecting reference genes** 
Once a set of reference genes and their orthologs in selected species were identified, the next step was to identify the top 10 most stable reference genes across all samples in all 16 selected species (outgroups and Cucurbitaceae). For this, the Python script **find_stable_reference.py** was used. The script first identifies the top 10 most stable reference genes from the selected references based on those with the least variable expression and prints those. This is followed by plotting the relative expression of these selected top 10 genes in all samples separately for each species. This figure can be found in the study as **_Fig. S19_**.  

**2) Normalizing count tables using reference genes** 
**normalize_CTs.py** was used to normalize all count tables based on the selected reference genes.
```
Usage
python3 normalize_CTs.py
```
It is an interactive-python script that asks for path to input folder containing count tables followed by path to output folder. 

**3) Plotting expression of Pigment biosynthesis as seen in Fig. 5 in the study**


# Data Availability
All data sets underlying this study are publicly available. The data supporting the results and conclusions are included in the article and its supporting information. All sequences and additional dataset files are accessible via bonndata https://doi.org/10.60507/FK2/ZNVSLA.
Customised Python scripts for the analyses in this study are available through GitHub https://github.com/NancyChoudhary28/Cucurbits.

# Citation
If you find the code useful in your research, please cite:

Out of the blue: Family-wide loss of anthocyanin biosynthesis in Cucurbitaceae. Nancy Choudhary, Marie Hagedorn, Boas Pucker. bioRxiv 2025.10.06.680802; doi: https://doi.org/10.1101/2025.10.06.680802

