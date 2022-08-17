# Get data from GISAID and create distance matrix and phylogenetic tree
#
# Prerequisites: GISAID account, Anaconda or Miniconda, and RStudio
#
# 1. You will need a GISAID account to download the data. Do that first.
#    https://gisaid.org/register/
# 2. Install parsnp, harvesttools, and snp-dists with conda. For example:
#
#    $ conda create -y --name parsnp-env
#    $ conda activate parsnp-env
#    $ conda install -y -c bioconda parsnp
#    $ conda install -y -c bioconda harvesttools
#    $ conda install -y -c bioconda snp-dists
#
#    Then launch RStudio and open this script in RStudio and run it.
#
#    You may need to modify the parsnp_path variable below so R can find it.

# Load packages, installing as needed
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(rstudioapi, lubridate, dplyr, ape, fs, here)
pacman::p_load_gh("Wytamma/GISAIDR")
pacman::p_load_gh("deohs/folders")

# Define path to parsnp, harvesttools, and snp-dists
bin_path <- file.path(path_home_r(), 'miniconda3/envs/parsnp-env/bin')

# Setup folders
folders <- get_folders()
data_dir <- here(folders$data)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Get GISAID username and password
username <- rstudioapi::askForPassword("GISAID username")
password <- rstudioapi::askForPassword("GISAID password")
credentials <- login(username = username, password = password)

# Define a function for writing fasta sequences to files
write_fasta <- function (df, data_dir = '.') {
    sapply(1:nrow(df), function(x) {
        if (!is.na(df$accession_id[x]) & !is.na(df$sequence[x])) {
            id <- gsub('EPI_ISL_', '', df$accession_id[x])
            writeLines(c(c(paste0('>', id, '_', df$host_type[x])), df$sequence[x]), 
                       con = file.path(data_dir, paste0(id, '_', df$host_type[x], '.fasta')))
        }
    })
}

# Define a vector of 5 GISAID IDs to be downloaded
ids5 <- c('EPI_ISL_7845318', 'EPI_ISL_7845316', 'EPI_ISL_7845317', 
          'EPI_ISL_7845315', 'EPI_ISL_8897004')

# Get a random sample of 20 accession IDs from query
df_ids <- query(
    credentials = credentials, 
    location = "North America / USA / Washington / King County", 
    from = "2021-04-23", 
    to = "2021-08-18",
    complete = TRUE,
    fast = TRUE
)

# Get random sample of 20 IDs
set.seed(1)
ids <- sample(df_ids$accession_id, 20)

# Add additional 5 IDs
ids <- c(ids, ids5)

# Download sequences
df_seqs <- download(
    credentials = credentials, 
    list_of_accession_ids = ids, 
    get_sequence = TRUE
)

# Add host_type variable to use in sequence label and fasta filename
df_seqs <- df_seqs %>% 
    mutate(host_type = ifelse(host == 'Felis catus', 'Feline',
                              ifelse(host == 'Canis lupus familiaris', 'Canine', 
                                     gsub(' ', '_', host))))

# Save metadata as CSV
csv_file <- file.path(data_dir, 'GISAID_metadata.csv')
write.csv(df_seqs %>% select(-sequence), csv_file, row.names = FALSE)

# Split sequence data into a single fasta file per accession ID
fasta_dir <- file.path(data_dir, 'genomes_fasta')
dir.create(fasta_dir, showWarnings = FALSE, recursive = TRUE)
res <- write_fasta(df_seqs, fasta_dir)

# Get the reference genome and save sequence in fasta format
reference_genome_id <- "NC_045512.2"
reference_genome <- read.GenBank(reference_genome_id)
ref_file <- file.path(data_dir, paste0(reference_genome_id, "_Reference.fasta"))
write.FASTA(reference_genome, ref_file) 

# Create output folder for parsnp
results_dir <- here(folders$results)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Run parsnp to align the genomes and produce a phylogenetic tree file
parsnp_cmd <- 
    paste(file.path(bin_path, "parsnp"), "-p 4 -c -r", ref_file, 
          "-d", fasta_dir, "-o", results_dir)
res <- system(parsnp_cmd, intern = TRUE)
writeLines(res, file.path(results_dir, "parsnp_output.txt"))

# Relabel parsnip.tree by removing filename suffixes from genome labels
parsnp_tree_fn <- file.path(results_dir, "parsnp.tree")
parsnp_tree <- readLines(parsnp_tree_fn)
parsnp_tree <- gsub('\\.(?:fa(?:sta)?|gbk\\.fna|ref)', '', parsnp_tree)
writeLines(parsnp_tree, parsnp_tree_fn)

# Run harvesttools to create a FASTA file containing the aligned genomes
harvesttools_cmd <- 
    paste(file.path(bin_path, "harvesttools"), 
          "-i", file.path(results_dir, "parsnp.ggr"), 
          "-M", file.path(results_dir, "parsnp.aln"))
res <- system(harvesttools_cmd, intern = TRUE)

# Run snp-dists to create a distance matrix file with pairwise snp distances
snp_dists_cmd <- 
    paste(file.path(bin_path, "snp-dists"), 
          "-b", file.path(results_dir, "parsnp.aln"))
res <- system(snp_dists_cmd, intern = TRUE)
res <- gsub('\\.(?:fa(?:sta)?|gbk\\.fna|ref)', '', res)
writeLines(res, file.path(results_dir, "distances.tab"))

# Make distance matrix with normalized snp distances (from tree) and save
tree <- read.tree(parsnp_tree_fn)
PatristicDistMatrix <- cophenetic(tree)
write.csv(PatristicDistMatrix, file.path(results_dir, "distances.csv"))

# Create output folder for figures
figures_dir <- here(folders$figures)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Make phylogenetic tree and save
pdf(file = file.path(figures_dir, "parsnp_tree.pdf"), width = 8, height = 8)
plot(tree, type = "fan", cex = 0.6)
dev.off()
