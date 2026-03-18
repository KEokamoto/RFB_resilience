# Assumes 00_setup.R has already been sourced

# Import sample metadata, feature table, taxonomy, and tree
sample_metadata <- read_excel("data/raw/sample-metadata.xlsx", col_names = TRUE)
otu <- read_excel("data/raw/feature-table.xlsx", skip = 1)
taxonomy <- read_excel("data/raw/taxonomy.xlsx", skip = 1)
TREE <- read_tree("data/raw/tree.nwk")

# Clean metadata
metadata <- sample_metadata %>%
  mutate(
    Description = if_else(
      str_detect(Description, "WholeWheatFresh|OatFresh"),
      Description,
      str_replace_all(Description, "F", "G")
    ),
    Treatment = str_extract(Description, "ControlA|ExpA|ControlB|ExpB"),
    Generation = str_extract(Description, "G[0-3]")
  )

metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1, drop = FALSE]

# Clean taxonomy
taxonomy <- as.data.frame(taxonomy)
rownames(taxonomy) <- taxonomy[, 1]
taxonomy <- taxonomy[, -1, drop = FALSE]

tax <- taxonomy %>%
  select(Taxon) %>%
  separate(
    Taxon,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = "; ",
    fill = "right",
    extra = "merge"
  )

tax.clean <- data.frame(
  Kingdom = str_replace(tax[[1]], "d__", ""),
  Phylum  = str_replace(tax[[2]], "p__", ""),
  Class   = str_replace(tax[[3]], "c__", ""),
  Order   = str_replace(tax[[4]], "o__", ""),
  Family  = str_replace(tax[[5]], "f__", ""),
  Genus   = str_replace(tax[[6]], "g__", ""),
  Species = str_replace(tax[[7]], "s__", ""),
  stringsAsFactors = FALSE
)

row.names(tax.clean) <- row.names(tax)
tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean == "__"] <- ""

for (i in seq_len(nrow(tax.clean))) {
  if (tax.clean[i, 7] != "") {
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i, 2] == "") {
    kingdom <- paste("Unclassified", tax.clean[i, 1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i, 3] == "") {
    phylum <- paste("Unclassified", tax.clean[i, 2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i, 4] == "") {
    class <- paste("Unclassified", tax.clean[i, 3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i, 5] == "") {
    order <- paste("Unclassified", tax.clean[i, 4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i, 6] == "") {
    family <- paste("Unclassified", tax.clean[i, 5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i, 7] == "") {
    tax.clean$Species[i] <- paste("Unclassified", tax.clean$Genus[i], sep = " ")
  }
}

# Process OTU table
otu_df <- as.data.frame(otu)
rownames(otu_df) <- otu_df[, 1]
otu_df <- otu_df[, -1, drop = FALSE]

# Create phyloseq object
OTU <- otu_table(as.matrix(otu_df), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(metadata)

ps <- phyloseq(OTU, TAX, SAMPLE, TREE)

# Subset to AdultGut
ps_AdultGut <- subset_samples(ps, Tissue == "AdultGut")

# Ensure Treatment and Generation exist in AdultGut metadata
sample_df <- data.frame(sample_data(ps_AdultGut))

if (!"Treatment" %in% colnames(sample_df) || any(is.na(sample_df$Treatment))) {
  sample_df$Treatment <- str_extract(sample_df$Description, "ControlA|ExpA|ControlB|ExpB")
}

if (!"Generation" %in% colnames(sample_df) || any(is.na(sample_df$Generation))) {
  sample_df$Generation <- str_extract(sample_df$Description, "G[0-3]")
}

sample_df$Treatment <- factor(sample_df$Treatment, levels = trt_levels)
sample_df$Generation <- factor(sample_df$Generation, levels = gen_levels)

sample_data(ps_AdultGut) <- sample_df

# G0-only object for baseline figure
ps_AdultGut_G0 <- subset_samples(ps_AdultGut, Generation == "G0")
ps_AdultGut_G0 <- prune_taxa(taxa_sums(ps_AdultGut_G0) > 0, ps_AdultGut_G0)

# Richness table for downstream analyses
richness_df <- estimate_richness(ps_AdultGut, measures = c("Observed", "Shannon"))
meta_df <- data.frame(sample_data(ps_AdultGut))

rich_combined <- cbind(meta_df, richness_df)
rich_combined$Generation <- factor(rich_combined$Generation, levels = gen_levels)
rich_combined$Treatment <- factor(rich_combined$Treatment, levels = trt_levels)
