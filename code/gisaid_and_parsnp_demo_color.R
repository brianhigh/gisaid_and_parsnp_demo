# Get data from GISAID and create distance matrix and phylogenetic tree
#
# Prerequisites: GISAID account, Anaconda or Miniconda, and RStudio
#
# 1. You will need a GISAID account to download the data. Do that first.
#    https://gisaid.org/register/
# 2. Install parsnp, harvesttools, and snp-dists with conda. For example:
#
#    $ conda create -y -c bioconda -n parsnp-env parsnp harvesttools snp-dists
#    $ conda activate parsnp-env
#
#    Then launch RStudio and open this script in RStudio and run it.
#
#    You may need to modify the bin_path variable below so R can find it.

# Load packages, installing as needed
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(rstudioapi, lubridate, dplyr, ape, fs, here, 
               ggtree, ggimage, RColorBrewer, emojifont, rsvg)
#pacman::p_install_gh("Wytamma/GISAIDR", force = TRUE)
pacman::p_load_gh("Wytamma/GISAIDR")
pacman::p_load_gh("deohs/folders")

# Define path to parsnp, harvesttools, and snp-dists
if (Sys.getenv("CONDA_PREFIX") == "") {
  bin_path <- '.'
  env_name <- 'parsnp-env'
  env_paths <- c('.conda/envs', 'miniconda3/envs')
  for (env_path in env_paths) {
    bin_dir <- file.path(path_home_r(), env_path, env_name, 'bin')
    if (dir.exists(bin_dir)) bin_path <- bin_dir
  }
} else {
  bin_path <- file.path(Sys.getenv("CONDA_PREFIX"), 'bin')
}

# Setup folders
conf <- here::here('folders.yml')
folders <- get_folders(conf)
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
                                   gsub(' ', '_', host))),
         tip.label = paste0(gsub('^EPI_ISL_', '', accession_id), '_', host_type))

# Save metadata as CSV
csv_file <- file.path(data_dir, 'GISAID_metadata.csv')
write.csv(df_seqs %>% select(-sequence), csv_file, row.names = FALSE)

# Split sequence data into a single fasta file per accession ID
fasta_dir <- file.path(data_dir, 'genomes_fasta')
dir.create(fasta_dir, showWarnings = FALSE, recursive = TRUE)
res <- write_fasta(df_seqs, fasta_dir)

# Get the lineage A reference genome and save as a FASTA file
# ref_genome_df <- query(credentials, lineage = 'A',
#                        location = "Asia / China / Hubei / Wuhan",
#                        from = "2020-01-05", to = "2020-01-05")
# ref_genome_df$virus_name
# ## "hCoV-19/Wuhan/WH04/2020"
# ref_genome_df$accession_id[1]
# ## "EPI_ISL_406801"
# ref_genome_id <- ref_genome_df$accession_id[1]
ref_genome_id <- "EPI_ISL_406801"
ref_genome <- download(credentials, ref_genome_id, get_sequence = TRUE)
id <- paste0(gsub('EPI_ISL_', '', ref_genome_id), '_Reference')
ref_file <- file.path(data_dir, paste0(id, ".fasta"))
writeLines(c(c(paste0('>', id)), ref_genome$sequence), con = ref_file)

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
tree <- root(tree, outgroup = grep('Reference', tree$tip.label))
tree$tip.label <- gsub("'", "", tree$tip.label)
PatristicDistMatrix <- cophenetic(tree)
write.csv(PatristicDistMatrix, file.path(results_dir, "distances.csv"))

# Create output folder for figures
figures_dir <- here(folders$figures)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Customize tree tip labels for plotting
palette <- "Dark2"
tree.labels <- df_seqs %>% 
  select(tip.label, 'Lineage' = pangolin_lineage)
tree.labels <- tibble(tip.label = tree$tip.label) %>% 
  left_join(tree.labels, by = 'tip.label') %>%
  mutate(Lineage = ifelse(tip.label == '406801_Reference', 'A', Lineage),
         Prefix = gsub('\\..*$', '', Lineage))
tree.prefixes <- tibble(Prefix = sort(unique(tree.labels$Prefix)))
tree.colors <- tibble(Prefix = tree.prefixes$Prefix, 
                      Color = brewer.pal(nrow(tree.prefixes), name = palette))
tree.labels <- tree.labels %>% left_join(tree.colors, by = 'Prefix') %>%
  mutate(tip.label = paste0(gsub('^\\d+_', '', tip.label), '_', Lineage))

# Assign image IDs for phylopic
# cat_uuid <- ggimage::phylopic_uid("Felis catus")$uid
# dog_uuid <- ggimage::phylopic_uid("Canis lupus")$uid
# human_uuid <- ggimage::phylopic_uid("Homo sapiens")$uid
cat_uuid <- 'f79384d6-2cee-47fb-bdae-25d491f82f9e'
dog_uuid <- '4d83a0cd-cf06-4a32-9a5a-0a6b644158c1'
human_uuid <- 'acf1cbec-f6ef-4a82-8ab5-be2f963b93f5'

tree.labels <- tree.labels %>% 
  mutate(uuid = ifelse(grepl('Canine', tip.label), dog_uuid, NA),
         uuid = ifelse(grepl('Feline', tip.label), cat_uuid, uuid),
         uuid = ifelse(grepl('Human|Reference', tip.label), human_uuid, uuid))

# Update tree tip labels in tree object
tree$tip.label <- tree.labels$tip.label

# Make phylogenetic tree plot and save
ggtree(tree) %<+% tree.labels +
  geom_tiplab(aes(colour = Prefix), offset = 0.0004) + xlim(NA, 0.02) + 
  geom_tiplab(aes(image = uuid, colour = Prefix), size = 0.015, 
              geom="phylopic") +
  theme(legend.position = "none") + scale_color_brewer(palette = palette)
ggsave(file.path(figures_dir, "parsnp_ggtree.pdf"), width = 11, height = 8.5)
