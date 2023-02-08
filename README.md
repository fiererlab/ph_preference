# ph_preference
A project to develop a method for predicting the optimal growth pH of microorganisms from 'omics data.

# Repo layout
- R folder contains the code for the bootstrapping process used to calculate the optimal pH for a strain using biogeographical information (pH_preference_bootstrapping_Panama_example.R).
- Python folder contains the machine learning (ML) code for prediction of pH preferences using genomic information.

# Outline
1. Gather datasets with 16S rDNA (or other marker gene data) that cover broad ranges in pH (or any environmental factor of interest)
2. Identify the pH where each ASV achieves its highest relative abundance (statistical inference of pH preference)
3. Identify associations between genes and pH preference
4. Match the the 16S rDNA reads to reference genomes
5. Annotate the genomes with functions
6. Build a ML model that uses the presence/absence of genes to predict pH preference 

# Example dataset
-  /data/ramonedaj/Panama/mapfile_panama.txt
-  /data/ramonedaj/Panama/otu_table_panama.txt
-  /data/ramonedaj/Panama/refseq.txt

# Genome database
-  /data/ramonedaj/pH_optima/bac120_metadata_r207.tsv
-  /data/mhoffert/genomes/GTDB_r207/pfam/

# Process
## Calculating the pH preference for ASVs in a dataset
1. Generate 1000 randomized distributions of the relative abundance of each ASV across samples with replacement (i.e. where each relative abundance value can be sampled more than once).
2. Calculated the maximum value for each of these distributions. 
3. Obtain 95% confidence intervals of these relative abundance maxima using the boot package (v1.3.28) in R. 
4. Match the extremes of these intervals of relative abundance maxima to the pH of the samples where these ASVs achieved these relative abundance values obtaining the range of pH in which a given ASV consistently achieves maximal abundance across randomizations. 
5. Remove ASVs with inferred pH preferences with ranges greater than 0.5 pH units.
6. Take the midpoint of the pH range as the estimated pH preference for a given ASV. 

## Matching ASVs to GTDB
### Download genome/16S rDNA SSU database
wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/bac120_ssu_reps_r207.tar.gz
tar -xzvf bac120_ssu_reps_r207.tar.gz
### Transform database to a vsearch reference database
vsearch --makeudb_usearch bac120_ssu_reps_r207.fna -output bac120_ssu_reps_r207.udb
### Align the ASVs representative sequences to the reference database
vsearch --usearch_global refseq.fasta --db bac120_ssu_reps_r207.udb --strand both --notrunclabels --iddef 0 --id 0.99 --maxrejects 100 --maxaccepts 100 --blast6out aligned_ssu.tsv --threads 16

## Annotate the genomes with functions

