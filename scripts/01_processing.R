trt_levels <- c("ControlA","ExpA","ControlB","ExpB")
trt_shapes <- c("ControlA"=16, "ExpA"=17, "ControlB"=1, "ExpB"=2)

#===========Data Import and Processing================
# Import sample metadata, OTU table, taxonomy, and tree
sample_metadata <- read_excel("sample-metadata.xlsx", col_names = TRUE)
otu <- read_excel("feature-table.xlsx", skip = 1)
taxonomy <- read_excel("taxonomy.xlsx", skip = 1)
TREE <- read_tree("tree.nwk")

# Replace all occurrences of "F" with "G" in the "Description" column of metadata
metadata <- sample_metadata %>%
  mutate(Description = if_else(
    str_detect(Description, "WholeWheatFresh|OatFresh"),
    Description,
    str_replace_all(Description, "F", "G")
  ))


metadata <- metadata %>%
  
  mutate(
    
    Treatment = str_extract(Description, "ControlA|ExpA|ControlB|ExpB"),
    
    Generation = str_extract(Description, "G[0-3]")
    
  )

# Convert metadata to a dataframe and set rownames
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

# Process taxonomy: convert to dataframe, set rownames, and clean taxonomic strings
taxonomy <- as.data.frame(taxonomy)
rownames(taxonomy) <- taxonomy[,1]
taxonomy <- taxonomy[,-1]

tax <- taxonomy %>%
  select(Taxon) %>%
  separate(Taxon, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = "; ", 
           fill = "right", 
           extra = "merge")

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
tax.clean[tax.clean=="__"] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified", tax.clean$Genus[i], sep = " ")
  }
}

# Process OTU table: convert to dataframe and set rownames
otu_df <- as.data.frame(otu)
rownames(otu_df) <- otu_df[,1]
otu_df <- otu_df[,-1]

# Create phyloseq object from OTU, Taxonomy, Sample metadata, and Tree
OTU <- otu_table(as.matrix(otu_df), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(metadata)
ps <- phyloseq(OTU, TAX, SAMPLE, TREE)

# Subset to AdultGut samples
ps_AdultGut <- subset_samples(ps, Tissue == "AdultGut")

# Estimate richness (Observed & Shannon) on AdultGut phyloseq object
richness_df   <- estimate_richness(ps_AdultGut, measures = c("Observed", "Shannon"))

# Pull out the sample metadata
meta_df       <- data.frame(sample_data(ps_AdultGut))

# Combine into one data frame
rich_combined <- cbind(meta_df, richness_df)

# Make sure Generation is a factor in the right order
rich_combined$Generation <- factor(rich_combined$Generation, levels = c("G0","G1","G2","G3"))

# (Re)ensure that Treatment and Generation columns are present in ps_AdultGut
sample_df <- data.frame(sample_data(ps_AdultGut))
if (!"Treatment" %in% colnames(sample_df) || any(is.na(sample_df$Treatment))) {
  sample_df$Treatment <- str_extract(sample_df$Description, "ControlA|ExpA|ControlB|ExpB")
}
if (!"Generation" %in% colnames(sample_df) || any(is.na(sample_df$Generation))) {
  sample_df$Generation <- str_extract(sample_df$Description, "G[0-3]")
}
sample_data(ps_AdultGut) <- sample_df
