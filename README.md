#### This repository contains the code used in the article:   
# [Out of the blue: Family-wide loss of anthocyanin biosynthesis in Cucurbitaceae](https://doi.org/10.1101/2025.10.06.680802)
https://doi.org/10.1101/2025.10.06.680802

---

**Software requirements**
- Python ≥ 3.8
- pandas
- seaborn
- matplotlib
- scipy
- numpy
- [KIPEs3](https://github.com/bpucker/KIPEs)
- [AGAT](https://github.com/NBISweden/AGAT)
- [Kraken](https://github.com/DerrickWood/kraken)
- [Bracken](https://github.com/jenniferlu717/Bracken)

---

## Getting datasets and identifying flavonoid and carotenoid biosynthesis genes 

### 1) Retrieving data and extracting coding sequences and longest isoforms

After downloading genome sequences (FASTA format) and gene annotations (GFF3 format), the script **run_agat_wrapper.sh** was used to extract coding sequences. In addition, the script outputs coding sequences corresponding to the longest-isoform per gene.

   ```
   Input:
   GENOME_DIR         directory containing all genome FASTA files
   GFF_DIR            directory containing all GFF3 files

   Output:
   CDS_DIR            directory containing all coding sequence FASTA files
   SINGLEISO_GFF_DIR  directory containing GFF3 files with only the longest isoform retained
   SINGLEISO_CDS_DIR  directory containing single-isoform coding sequence FASTA files
   LOG_DIR            directory containing all log files
   ```

_*Note: The base names of the genome FASTA and corresponding GFF3 files must match. Output files retain the same base names in their respective output directories. Please update the input and output directory paths in the wrapper script before running._

---

### 2) Translating coding sequences into polypeptide sequences

The script **transeq.py**, available [here](https://github.com/bpucker/PBBtools/tree/main/transeq), was used to translate all coding sequences into polypeptide sequences.

---

### 3) Identifying contaminated datasets (genomic and transcriptomic)

The script **run_kraken_braken_wrapper.sh** was used to classify coding sequences as contaminated or non-contaminated. Datasets in which more than 10% of coding sequences were classified as non-Cucurbitaceae were considered contaminated and excluded from downstream analyses.

   ```
   Input:
   input_dir          path to directory containing coding sequences (FASTA format)
   db_path            path to the Kraken classification database

   Output:
   output_base        path to the output directory
   ```
   
_*Note: Please update the input and output directory paths in the wrapper script before running._

---

### 4) Identifying pigment biosynthesis genes in non-contaminated datasets 

**KIPEs3.py**, available [here](https://github.com/bpucker/KIPEs), was used to identify flavonoid and carotenoid biosynthesis pathway genes from polypeptide sequences of all non-contaminated datasets.

---

### 5) Collecting candidate flavonoid biosynthesis genes for phylogenetic analysis

The script **collect_KIPEs_sequence.py** was used to extract candidate flavonoid biosynthesis gene sequences from KIPEs results for phylogenetic tree construction.
  ```
   Usage
   python3 collect_KIPEs_sequence.py
   --directory        path to the KIPEs result directory
   --enzymes          list of target genes (comma-separated if multiple)
   --output           output FASTA file
   --threshold        minimum percentage of conserved residues (default:70.0)
   ``` 

---

### 6) Preparing iTOL annotation files for visualizing ortholog clades

The script **get_conserved_residue_perc.py** was used to generate an iTOL gradient annotation file representing the percentage of conserved residues required for the said gene function. 
```
Usage
python3 get_conserved_residue_perc.py
-k               path to the KIPEs result directory
-i               input file containing fasta headers (without '>')
-o               path to the output file
-g               candidate gene name to match in summary.txt (e.g. DFR)
```
The input FASTA headers must correspond to candidates of a single gene, specified using the `-g` flag. After generating the annotation file, minor manual edits are required to conform to iTOL annotation file guidelines. An example file (**_DFR_gradient.txt_**) is provided for reference. 

---

## Expression analyses

### 1) Selecting reference genes

After identifying a set of reference genes and their orthologs across selected species, the top 10 most stable reference genes were identified across all samples from 16 selected species (including outgroups and Cucurbitaceae). This was done using **find_stable_reference.py**, which selects genes with the least variable expression and plots their relative expression across all samples for each species. The resulting figure is shown in the supplementary file of the study as **_Fig. S19_**.  

---

### 2) Normalizing count tables using reference genes

The script **normalize_CTs.py** was used to normalize all count tables based on the selected reference genes.
```
Usage
python3 normalize_CTs.py
```
This is an interactive Python script that prompts the user to provide the input directory containing count tables and the output directory for normalized count tables. 

---

### 3) Plotting pigment biosynthesis gene expression

For this, two approaches were used: 
(a) plotting the expression of each key pathway gene (flavonoid and carotenoid) per species, and 
(b) plotting the combined expression of each key gene per lineage (Cucurbitaceae vs outgroups).
For species-level plots, **expression_violin_plot_per_species.py** was used (separately for flavonoid and carotenoid genes). For lineage-level plots, **expression_violin_plot_per_lineage.py** was used. The four resulting plots were manually combined using Inkscape. 

These scripts are specific to the species and gene candidates used in the study. To reuse them, users should modify the species and gene candidates lists within the scripts accordingly. 

---

# Data Availability
All data sets underlying this study are publicly available. The data supporting the results and conclusions are included in the article and its supporting information. All sequences and additional dataset files are accessible via bonndata https://doi.org/10.60507/FK2/ZNVSLA.
Customised Python scripts for the analyses in this study are available through GitHub https://github.com/NancyChoudhary28/Cucurbits.

---

# Citation
If you find the code useful in your research, please cite:

Out of the blue: Family-wide loss of anthocyanin biosynthesis in Cucurbitaceae. Nancy Choudhary, Marie Hagedorn, Boas Pucker. bioRxiv 2025.10.06.680802; doi: https://doi.org/10.1101/2025.10.06.680802

