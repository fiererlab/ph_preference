# optimum_ph
A project to develop a method for predicting the optimal growth pH of microorganisms from 'omics data.

# outline
1. We will identify some datasets with 16S (or other marker gene data)
2. Identify the pH where each ASVs is at its highest relative abundance
3. Match the ASVs to a GTDB (or other) genome
4. Annotate the genome with functions
5. Build a ML model that maps functions -> optimum pH 

# Datasets
## Panama
ASV table: `/data/mhoffert/lab-data/datasets/Oliverio_2020a/processed/03_tabletax/seqtab_final.tsv`
ASV fasta: `/data/mhoffert/lab-data/datasets/Oliverio_2020a/processed/03_tabletax/repset.fasta`
metadata:  `/data/mhoffert/lab-data/datasets/Oliverio_2020a/metadata`
## Australia
ASV table: `/data/mhoffert/lab-data/datasets/Oliverio_2020b/processed/03_tabletax/seqtab_final.tsv`
ASV fasta: `/data/mhoffert/lab-data/datasets/Oliverio_2020b/processed/03_tabletax/repset.fasta`
metadata:  `/data/mhoffert/lab-data/datasets/Oliverio_2020b/metadata`

# Genomics databases
Most of the files from the GTDB database are here: `/data/mhoffert/genomes/GTDB/`

# Process
## Calculating the optimal pH for a strain
The pH where an ASV is most abundant could be calculated by fitting a smooth spline to the data and picking the maximum value. I've put a notebook demonstrating this [here](notebooks/Panama_ML_analysis.ipynb):

## Matching ASVs to GTDB
If you haven't (in a conda environment!) install vsearch:
```
conda install -c bioconda vsearch
```
The 16S sequences are here: `/data/mhoffert/genomes/GTDB/bac120_ssu_reps.fna` (or something like that)  
First, I made the 16S sequences for GTDB release 202 into a Usearch-formatted database:
```vsearch --makeudb_usearch bac120_blah.fna -o bac120_blah.udb```
Then aligned the ASVs to the database with vsearch:
 ```
 vsearch --usearch_global [asvs_file] --db [udb] --id 0.95 --strand both --notrunclabels --outfmt 6 --output blah.txt --maxhits 1 --maxrejects 100 --threads 16
 ```
 And parse the file with the library of your choice. The file is in Blast format 6. The second column with ids like `RS_GCF_002514152.1` is the genomes: this is the subset that need to be annotated as inputs to the ML method
 ## Annotate the genomes with functions
 It will be something like:
 ```hmmsearch pfam.hmm [genome aa fasta]```
 And then some parsing of the results to grab the complete set of functions
