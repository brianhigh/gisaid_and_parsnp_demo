# Get data from GISAID and create distance matrix and phylogenetic tree
#
# Prerequisites: GISAID account, Anaconda or Miniconda, parsnp, and RStudio
#
# 1. You will need a GISAID account to download the data. Do that first.
#    https://gisaid.org/register/
# 2. Make sure parsnp has been installed in a conda environment, for example:
#
#    $ conda create -y --name parsnp-env
#    $ conda activate parsnp-env
#    $ conda install -y -c bioconda parsnp
#
#    Then launch RStudio and open this script in RStudio and run it.
#
#    You may need to modify the parsnp_path variable below so R can find it.

# Load packages, installing as needed
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(rstudioapi, lubridate, dplyr, ape, fs, here)
pacman::p_load_gh("Wytamma/GISAIDR")
pacman::p_load_gh("deohs/folders")

# Define path to parsnp
#parsnp_path <- file.path(path_home_r(), 'miniconda3/envs/parsnp-env/bin/parsnp')
parsnp_path <- 'parsnp'

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
df <- query(
    credentials = credentials, 
    location = "North America / USA / Washington / King County", 
    from = "2021-04-23", 
    to = "2021-08-18",
    complete = TRUE,
    fast = TRUE
)

# Get random sample of 20 IDs
set.seed(1)
ids <- sample(df$accession_id, 20)

# Add additional 5 IDs
ids <- c(ids, ids5)

# Download sequences
full_df_with_seq <- download(
    credentials = credentials, 
    list_of_accession_ids = ids, 
    get_sequence = TRUE
)

# Add host_type variable to use in sequence label and fasta filename
full_df_with_seq <- full_df_with_seq %>% 
    mutate(host_type = ifelse(host == 'Felis catus','Feline',
                              ifelse(host == 'Canis lupus familiaris','Canine', 
                                     gsub(' ', '_', host))))

# Save metadata as CSV
csv_file <- file.path(data_dir, 'GISAID_metadata.csv')
write.csv(full_df_with_seq %>% select(-sequence), csv_file, row.names = FALSE)

# Split sequence data into a single fasta file per accession ID
fasta_dir <- file.path(data_dir, 'genomes')
dir.create(fasta_dir, showWarnings = FALSE, recursive = TRUE)
res <- write_fasta(full_df_with_seq, fasta_dir)

# Get the reference genome and save sequence in fasta format
reference_genome_id <- "NC_045512.2"
reference_genome <- read.GenBank(reference_genome_id)
ref_file <- file.path(data_dir, paste0(reference_genome_id, "_Reference.fasta")
write.FASTA(reference_genome, ref_file) 

# Create output folder for parsnp
output_dir <- here(folders$results)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Run parsnp 
parsnp_cmd <- 
    paste(parsnp_path, "-p 4 -c -r", ref_file, 
          "-d", fasta_dir, "-o", output_dir)
res <- system(parsnp_cmd, intern = TRUE)
writeLines(res, file.path(output_dir, "parsnp_output.txt"))

# Relabel parsnip.tree by removing filename suffixes from genome labels
parsnp_tree_fn <- file.path(output_dir, "parsnp.tree")
parsnp_tree <- readLines(parsnp_tree_fn)
parsnp_tree <- gsub('\\.(?:fa(?:sta)?|gbk\\.fna|ref)', '', parsnp_tree)
writeLines(parsnp_tree, parsnp_tree_fn)

# Make distances matrix and save as CSV
tree <- read.tree(parsnp_tree_fn)
PatristicDistMatrix <- cophenetic(tree)
write.csv(PatristicDistMatrix, file.path(output_dir, "distances.csv"))

# Create output folder for figures
figures_dir <- here(folders$figures)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Make phylogenetic tree and save as PDF
pdf(file = file.path(figures_dir, "parsnp_tree.pdf"), width = 8, height = 8)
plot(tree)
dev.off()
