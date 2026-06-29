# Assumes 00_setup.R and 01_processing.R have already been sourced
# REVISED: integrates ANCOMBC2, G1-G3 analysis, pairwise Wilcoxon,
#          SIMPER, taxon overlap, G1G2 transition analyses

library(ANCOMBC)
library(microbiome)

# -------------------------
# Figure 3 alpha diversity stats
# -------------------------
alpha_df <- estimate_richness(ps_AdultGut, measures = c("Observed", "Shannon")) %>%
  rownames_to_column("SampleID") %>%
  left_join(
    data.frame(sample_data(ps_AdultGut)) %>% rownames_to_column("SampleID"),
    by = "SampleID"
  ) %>%
  mutate(
    Treatment  = factor(Treatment, levels = trt_levels),
    Generation = factor(Generation, levels = gen_levels)
  )

# Linear model diagnostics for alpha diversity
shannon_lm <- lm(Shannon ~ Generation * Treatment, data = alpha_df)
observed_lm <- lm(Observed ~ Generation * Treatment, data = alpha_df)

shapiro_shannon <- shapiro.test(residuals(shannon_lm))
shapiro_observed <- shapiro.test(residuals(observed_lm))

alpha_residual_normality <- tibble(
  Metric = c("Shannon diversity", "Observed richness"),
  W = c(unname(shapiro_shannon$statistic), unname(shapiro_observed$statistic)),
  p_value = c(shapiro_shannon$p.value, shapiro_observed$p.value)
)

alpha_long <- alpha_df %>%
  pivot_longer(
    cols = c("Observed", "Shannon"),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Observed", "Shannon"),
      labels = c("Observed richness", "Shannon diversity")
    )
  )

kw_results_fig3 <- alpha_long %>%
  group_by(Treatment, Metric) %>%
  summarise(
    p = kruskal.test(Value ~ Generation)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = paste0("Kruskal-Wallis p = ", signif(p, 3))
  )

# -------------------------
# Pairwise Wilcoxon tests for alpha diversity (Reviewer 3)
# -------------------------
pairwise_wilcox <- alpha_long %>%
  group_by(Treatment, Metric) %>%
  summarise(
    pairwise = list(
      pairwise.wilcox.test(Value, Generation, p.adjust.method = "BH")$p.value
    ),
    .groups = "drop"
  )

# -------------------------
# Beta diversity and PERMANOVA — full dataset
# -------------------------
bray_dist <- phyloseq::distance(ps_AdultGut, method = "bray")
meta <- data.frame(sample_data(ps_AdultGut))

adonis_terms <- vegan::adonis2(
  bray_dist ~ Generation + Treatment + Generation:Treatment,
  data = meta,
  permutations = 999,
  by = "terms"
)

adonis_margin <- vegan::adonis2(
  bray_dist ~ Generation * Treatment,
  data = meta,
  permutations = 999,
  by = "margin"
)

disp_gen <- vegan::betadisper(bray_dist, meta$Generation)
disp_anova <- anova(disp_gen)
disp_permutest <- vegan::permutest(disp_gen, pairwise = TRUE)

disp_trt <- vegan::betadisper(bray_dist, meta$Treatment)
disp_trt_permutest <- vegan::permutest(disp_trt)

# -------------------------
# G1-G3 restricted PERMANOVA and PERMDISP (Reviewer 1)
# -------------------------
ps_AdultGut_G1G3 <- subset_samples(ps_AdultGut, Generation %in% c("G1", "G2", "G3"))
bray_dist_G1G3   <- phyloseq::distance(ps_AdultGut_G1G3, method = "bray")
meta_G1G3        <- data.frame(sample_data(ps_AdultGut_G1G3))

adonis_G1G3 <- vegan::adonis2(
  bray_dist_G1G3 ~ Generation + Treatment + Generation:Treatment,
  data = meta_G1G3, permutations = 999, by = "terms"
)

disp_G1G3_gen <- vegan::betadisper(bray_dist_G1G3, meta_G1G3$Generation)
disp_G1G3_trt <- vegan::betadisper(bray_dist_G1G3, meta_G1G3$Treatment)
permutest_G1G3_gen <- vegan::permutest(disp_G1G3_gen)
permutest_G1G3_trt <- vegan::permutest(disp_G1G3_trt)

# -------------------------
# SIMPER — order level (Reviewer 2)
# -------------------------
ps_order_simper <- tax_glom(ps_AdultGut, taxrank = "Order")
otu_order_simper <- as.data.frame(t(otu_table(ps_order_simper)))
meta_simper      <- data.frame(sample_data(ps_order_simper))

simper_gen     <- vegan::simper(otu_order_simper, meta_simper$Generation, permutations = 999)
simper_summary <- summary(simper_gen)

# -------------------------
# Taxon overlap: G0 vs downstream generations (Reviewer 3)
# -------------------------
ps_order_overlap <- tax_glom(ps_AdultGut, taxrank = "Order")

orders_G0 <- taxa_names(prune_taxa(
  taxa_sums(subset_samples(ps_order_overlap, Generation == "G0")) > 0,
  subset_samples(ps_order_overlap, Generation == "G0")
))

orders_downstream <- taxa_names(prune_taxa(
  taxa_sums(subset_samples(ps_order_overlap, Generation != "G0")) > 0,
  subset_samples(ps_order_overlap, Generation != "G0")
))

overlap_n   <- length(intersect(orders_G0, orders_downstream))
pct_overlap <- overlap_n / length(orders_downstream) * 100
cat("% downstream orders also in G0:", round(pct_overlap, 1), "\n")

overlap_summary <- tibble(
  n_G0_orders = length(orders_G0),
  n_downstream_orders = length(orders_downstream),
  n_downstream_also_in_G0 = overlap_n,
  proportion_downstream_also_in_G0 = overlap_n / length(orders_downstream),
  percent_downstream_also_in_G0 = pct_overlap
)

# -------------------------
# Order-level phyloseq object for DA analyses
# -------------------------
ps_order <- tax_glom(ps_AdultGut, taxrank = "Order")

sample_data(ps_order)$Treatment <- factor(
  sample_data(ps_order)$Treatment,
  levels = c("ControlA", "ExpA", "ControlB", "ExpB")
)

sample_data(ps_order)$Generation <- factor(
  sample_data(ps_order)$Generation,
  levels = c("G0", "G1", "G2", "G3")
)

tax_map <- as.data.frame(tax_table(ps_order)) %>%
  rownames_to_column("taxon") %>%
  select(taxon, OrderName = Order) %>%
  mutate(
    OrderName = gsub("_[0-9]+$", "", OrderName),
    OrderName = ifelse(OrderName == "" | is.na(OrderName), "Unclassified", OrderName)
  )

# -------------------------
# DESeq2 pairwise function
# -------------------------
run_pairwise_deseq <- function(ps_obj, gen_val, trt1, trt0, lfc_name, padj_name) {
  samdf <- data.frame(sample_data(ps_obj))
  samdf$SampleID <- rownames(samdf)

  keep_samples <- samdf %>%
    dplyr::filter(Generation == gen_val, Treatment %in% c(trt0, trt1)) %>%
    dplyr::pull(SampleID)

  ps_sub <- prune_samples(keep_samples, ps_obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)

  sample_data(ps_sub)$Treatment <- factor(
    sample_data(ps_sub)$Treatment,
    levels = c(trt0, trt1)
  )

  dds <- phyloseq_to_deseq2(ps_sub, ~ Treatment)
  dds <- DESeq(dds, quiet = TRUE)

  res <- results(dds, contrast = c("Treatment", trt1, trt0))
  res_df <- as.data.frame(res)
  res_df$taxon <- rownames(res_df)

  out <- res_df[, c("taxon", "log2FoldChange", "padj")]
  colnames(out) <- c("taxon", lfc_name, padj_name)
  out
}

# -------------------------
# DESeq2 diet-shift contrasts at G2 and G3
# -------------------------
oat_G2   <- run_pairwise_deseq(ps_order, "G2", "ExpA", "ControlA", "lfc_oat_G2",   "padj_oat_G2")
oat_G3   <- run_pairwise_deseq(ps_order, "G3", "ExpA", "ControlA", "lfc_oat_G3",   "padj_oat_G3")
wheat_G2 <- run_pairwise_deseq(ps_order, "G2", "ExpB", "ControlB", "lfc_wheat_G2", "padj_wheat_G2")
wheat_G3 <- run_pairwise_deseq(ps_order, "G3", "ExpB", "ControlB", "lfc_wheat_G3", "padj_wheat_G3")

merged_df <- full_join(oat_G2, oat_G3, by = "taxon") %>%
  full_join(wheat_G2, by = "taxon") %>%
  full_join(wheat_G3, by = "taxon") %>%
  left_join(tax_map, by = "taxon")

alpha <- 0.05

merged_df <- merged_df %>%
  rowwise() %>%
  mutate(
    oat_same_direction = !is.na(lfc_oat_G2) & !is.na(lfc_oat_G3) &
      sign(lfc_oat_G2) == sign(lfc_oat_G3) &
      sign(lfc_oat_G2) != 0,
    oat_any_sig = any(c(padj_oat_G2, padj_oat_G3) < alpha, na.rm = TRUE),
    wheat_same_direction = !is.na(lfc_wheat_G2) & !is.na(lfc_wheat_G3) &
      sign(lfc_wheat_G2) == sign(lfc_wheat_G3) &
      sign(lfc_wheat_G2) != 0,
    wheat_any_sig = any(c(padj_wheat_G2, padj_wheat_G3) < alpha, na.rm = TRUE)
  ) %>%
  ungroup()

# -------------------------
# ANCOMBC2 pairwise function
# -------------------------
run_pairwise_ancombc <- function(ps_obj, gen_val, trt1, trt0, lfc_name, padj_name) {
  samdf <- data.frame(sample_data(ps_obj))
  samdf$SampleID <- rownames(samdf)

  keep_samples <- samdf %>%
    dplyr::filter(Generation == gen_val, Treatment %in% c(trt0, trt1)) %>%
    dplyr::pull(SampleID)

  ps_sub <- prune_samples(keep_samples, ps_obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)

  sample_data(ps_sub)$Treatment <- factor(
    sample_data(ps_sub)$Treatment,
    levels = c(trt0, trt1)
  )

  ancom_out <- ancombc2(
    data         = ps_sub,
    fix_formula  = "Treatment",
    p_adj_method = "BH",
    pseudo_sens  = TRUE,
    prv_cut      = 0.10,
    verbose      = FALSE
  )

  res <- ancom_out$res %>%
    select(
      taxon,
      lfc  = starts_with("lfc_Treatment"),
      padj = starts_with("q_Treatment")
    )

  colnames(res) <- c("taxon", lfc_name, padj_name)
  res
}

# -------------------------
# ANCOMBC2 diet-shift contrasts at G2 and G3
# -------------------------
anbc_oat_G2   <- run_pairwise_ancombc(ps_order, "G2", "ExpA", "ControlA", "anbc_lfc_oat_G2",   "anbc_padj_oat_G2")
anbc_oat_G3   <- run_pairwise_ancombc(ps_order, "G3", "ExpA", "ControlA", "anbc_lfc_oat_G3",   "anbc_padj_oat_G3")
anbc_wheat_G2 <- run_pairwise_ancombc(ps_order, "G2", "ExpB", "ControlB", "anbc_lfc_wheat_G2", "anbc_padj_wheat_G2")
anbc_wheat_G3 <- run_pairwise_ancombc(ps_order, "G3", "ExpB", "ControlB", "anbc_lfc_wheat_G3", "anbc_padj_wheat_G3")

anbc_merged <- full_join(anbc_oat_G2, anbc_oat_G3, by = "taxon") %>%
  full_join(anbc_wheat_G2, by = "taxon") %>%
  full_join(anbc_wheat_G3, by = "taxon") %>%
  left_join(tax_map, by = "taxon")

# -------------------------
# Concordance: DESeq2 vs ANCOMBC2
# Must run BEFORE labels_row and heatmap
# -------------------------
anbc_concordance <- anbc_merged %>%
  select(taxon, OrderName,
         anbc_lfc_oat_G2, anbc_padj_oat_G2,
         anbc_lfc_oat_G3, anbc_padj_oat_G3,
         anbc_lfc_wheat_G2, anbc_padj_wheat_G2,
         anbc_lfc_wheat_G3, anbc_padj_wheat_G3) %>%
  rowwise() %>%
  mutate(
    anbc_oat_sig = any(c(anbc_padj_oat_G2, anbc_padj_oat_G3) < 0.05, na.rm = TRUE),
    anbc_oat_dir = !is.na(anbc_lfc_oat_G2) & !is.na(anbc_lfc_oat_G3) &
      sign(anbc_lfc_oat_G2) == sign(anbc_lfc_oat_G3),
    anbc_wheat_sig = any(c(anbc_padj_wheat_G2, anbc_padj_wheat_G3) < 0.05, na.rm = TRUE),
    anbc_wheat_dir = !is.na(anbc_lfc_wheat_G2) & !is.na(anbc_lfc_wheat_G3) &
      sign(anbc_lfc_wheat_G2) == sign(anbc_lfc_wheat_G3)
  ) %>%
  ungroup() %>%
  select(taxon, anbc_oat_sig, anbc_oat_dir, anbc_wheat_sig, anbc_wheat_dir)

# Drop any existing concordance columns to avoid join conflicts on re-runs
merged_df <- merged_df %>%
  select(-any_of(c("anbc_oat_sig", "anbc_oat_dir", "anbc_wheat_sig", "anbc_wheat_dir",
                   "oat_concordant", "wheat_concordant", "oat_arrow", "wheat_arrow")))

merged_df <- merged_df %>%
  left_join(anbc_concordance, by = "taxon") %>%
  mutate(
    oat_concordant   = oat_same_direction & oat_any_sig &
      !is.na(anbc_oat_sig) & anbc_oat_sig & anbc_oat_dir,
    wheat_concordant = wheat_same_direction & wheat_any_sig &
      !is.na(anbc_wheat_sig) & anbc_wheat_sig & anbc_wheat_dir,
    # Concordance-filtered arrows (all "" since ANCOMBC2 found nothing)
    oat_arrow = case_when(
      oat_concordant & mean(c(lfc_oat_G2, lfc_oat_G3), na.rm = TRUE) > 0 ~ "\u2191O",
      oat_concordant & mean(c(lfc_oat_G2, lfc_oat_G3), na.rm = TRUE) < 0 ~ "\u2193O",
      TRUE ~ ""
    ),
    wheat_arrow = case_when(
      wheat_concordant & mean(c(lfc_wheat_G2, lfc_wheat_G3), na.rm = TRUE) > 0 ~ "\u2191W",
      wheat_concordant & mean(c(lfc_wheat_G2, lfc_wheat_G3), na.rm = TRUE) < 0 ~ "\u2193W",
      TRUE ~ ""
    )
  )

# -------------------------
# Heatmap matrix
# -------------------------
ps_order_rel <- transform_sample_counts(ps_order, function(x) x / sum(x))
sample_data(ps_order_rel)$TrtGen <- paste0(
  sample_data(ps_order_rel)$Treatment, "_",
  sample_data(ps_order_rel)$Generation
)

heatmap_mat <- psmelt(ps_order_rel) %>%
  group_by(Order, TrtGen) %>%
  summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = TrtGen, values_from = MeanAbundance, values_fill = 0) %>%
  as.data.frame()

rownames(heatmap_mat) <- heatmap_mat$Order
heatmap_mat <- heatmap_mat[, -1, drop = FALSE]

desired_col_order <- c(
  "ControlA_G0", "ControlA_G1", "ControlA_G2", "ControlA_G3",
  "ExpA_G0",     "ExpA_G1",     "ExpA_G2",     "ExpA_G3",
  "ControlB_G0", "ControlB_G1", "ControlB_G2", "ControlB_G3",
  "ExpB_G0",     "ExpB_G1",     "ExpB_G2",     "ExpB_G3"
)
desired_col_order <- desired_col_order[desired_col_order %in% colnames(heatmap_mat)]
heatmap_mat       <- heatmap_mat[, desired_col_order, drop = FALSE]
heatmap_mat       <- heatmap_mat[apply(heatmap_mat, 1, var) > 0, , drop = FALSE]
heatmap_mat_filtered <- heatmap_mat[rowMeans(heatmap_mat) > 0.001, , drop = FALSE]

# labels_row uses concordance-filtered merged_df (computed above)
labels_row <- rownames(heatmap_mat_filtered) %>%
  sapply(function(ord) {
    row_i    <- merged_df %>% filter(OrderName == ord)
    oat_lab  <- if (nrow(row_i) > 0) unique(row_i$oat_arrow[row_i$oat_arrow != ""])  else character(0)
    wheat_lab <- if (nrow(row_i) > 0) unique(row_i$wheat_arrow[row_i$wheat_arrow != ""]) else character(0)
    arrows   <- c(oat_lab, wheat_lab)
    arrows   <- arrows[arrows != ""]
    if (length(arrows) > 0) paste(ord, paste(arrows, collapse = " ")) else ord
  })

col_labels <- gsub("_", " ", colnames(heatmap_mat_filtered))

# -------------------------
# G1 flour-type DESeq2 function
# -------------------------
run_flourtype_G1_deseq <- function(ps_obj) {
  samdf <- data.frame(sample_data(ps_obj))
  samdf$SampleID <- rownames(samdf)

  keep_samples <- samdf %>%
    dplyr::filter(Generation == "G1") %>%
    dplyr::pull(SampleID)

  ps_sub <- prune_samples(keep_samples, ps_obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)

  sample_data(ps_sub)$FlourType <- factor(
    ifelse(sample_data(ps_sub)$Treatment %in% c("ControlA", "ExpA"), "Wheat", "Oat"),
    levels = c("Wheat", "Oat")
  )

  dds <- phyloseq_to_deseq2(ps_sub, ~ FlourType)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("FlourType", "Oat", "Wheat"))
  res_df <- as.data.frame(res)
  res_df$taxon <- rownames(res_df)
  res_df[, c("taxon", "log2FoldChange", "padj")]
}

# -------------------------
# G1 flour-type ANCOMBC2 function
# -------------------------
run_flourtype_G1_ancombc <- function(ps_obj) {
  samdf <- data.frame(sample_data(ps_obj))
  samdf$SampleID <- rownames(samdf)

  keep_samples <- samdf %>%
    dplyr::filter(Generation == "G1") %>%
    dplyr::pull(SampleID)

  ps_sub <- prune_samples(keep_samples, ps_obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)

  sample_data(ps_sub)$FlourType <- factor(
    ifelse(sample_data(ps_sub)$Treatment %in% c("ControlA", "ExpA"), "Wheat", "Oat"),
    levels = c("Wheat", "Oat")
  )

  ancom_out <- ancombc2(
    data         = ps_sub,
    fix_formula  = "FlourType",
    p_adj_method = "BH",
    pseudo_sens  = TRUE,
    prv_cut      = 0.10,
    verbose      = FALSE
  )

  ancom_out$res %>%
    select(taxon,
           lfc  = starts_with("lfc_FlourType"),
           padj = starts_with("q_FlourType"))
}

# -------------------------
# G1->G2 transition DESeq2 function
# -------------------------
run_G1G2_transition_deseq <- function(ps_obj, treatment_val, lfc_name, padj_name) {
  samdf <- data.frame(sample_data(ps_obj))
  samdf$SampleID <- rownames(samdf)

  keep_samples <- samdf %>%
    dplyr::filter(Treatment == treatment_val, Generation %in% c("G1", "G2")) %>%
    dplyr::pull(SampleID)

  ps_sub <- prune_samples(keep_samples, ps_obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)

  sample_data(ps_sub)$Generation <- factor(
    sample_data(ps_sub)$Generation, levels = c("G1", "G2")
  )

  dds <- phyloseq_to_deseq2(ps_sub, ~ Generation)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("Generation", "G2", "G1"))
  res_df <- as.data.frame(res)
  res_df$taxon <- rownames(res_df)
  out <- res_df[, c("taxon", "log2FoldChange", "padj")]
  colnames(out) <- c("taxon", lfc_name, padj_name)
  out
}

# -------------------------
# G1->G2 transition ANCOMBC2 function
# -------------------------
run_G1G2_transition_ancombc <- function(ps_obj, treatment_val, lfc_name, padj_name) {
  samdf <- data.frame(sample_data(ps_obj))
  samdf$SampleID <- rownames(samdf)

  keep_samples <- samdf %>%
    dplyr::filter(Treatment == treatment_val, Generation %in% c("G1", "G2")) %>%
    dplyr::pull(SampleID)

  ps_sub <- prune_samples(keep_samples, ps_obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)

  sample_data(ps_sub)$Generation <- factor(
    sample_data(ps_sub)$Generation, levels = c("G1", "G2")
  )

  ancom_out <- ancombc2(
    data         = ps_sub,
    fix_formula  = "Generation",
    p_adj_method = "BH",
    pseudo_sens  = TRUE,
    prv_cut      = 0.10,
    verbose      = FALSE
  )

  res <- ancom_out$res %>%
    select(taxon,
           lfc  = starts_with("lfc_Generation"),
           padj = starts_with("q_Generation"))
  colnames(res) <- c("taxon", lfc_name, padj_name)
  res
}

# -------------------------
# Run G1 flour-type and G1->G2 transition contrasts
# -------------------------
g1_flour_deseq   <- run_flourtype_G1_deseq(ps_order)
g1_flour_ancombc <- run_flourtype_G1_ancombc(ps_order)

expA_G1G2_deseq   <- run_G1G2_transition_deseq(ps_order,   "ExpA", "lfc_expA_G1G2",      "padj_expA_G1G2")
expB_G1G2_deseq   <- run_G1G2_transition_deseq(ps_order,   "ExpB", "lfc_expB_G1G2",      "padj_expB_G1G2")
expA_G1G2_ancombc <- run_G1G2_transition_ancombc(ps_order, "ExpA", "anbc_lfc_expA_G1G2", "anbc_padj_expA_G1G2")
expB_G1G2_ancombc <- run_G1G2_transition_ancombc(ps_order, "ExpB", "anbc_lfc_expB_G1G2", "anbc_padj_expB_G1G2")

g1g2_merged <- full_join(expA_G1G2_deseq, expB_G1G2_deseq, by = "taxon") %>%
  left_join(tax_map, by = "taxon") %>%
  left_join(expA_G1G2_ancombc, by = "taxon") %>%
  left_join(expB_G1G2_ancombc, by = "taxon")

# -------------------------
# Save all numeric outputs
# -------------------------
dir.create("results", showWarnings = FALSE)

# Figure 3 — Alpha diversity
write.csv(
  kw_results_fig3,
  "results/fig3_alpha_diversity_kruskal.csv",
  row.names = FALSE
)

write.csv(
  alpha_residual_normality,
  "results/fig3_alpha_residual_normality.csv",
  row.names = FALSE
)

write.csv(
  bind_rows(lapply(seq_len(nrow(pairwise_wilcox)), function(i) {
    mat <- pairwise_wilcox$pairwise[[i]]
    as.data.frame(mat) %>%
      rownames_to_column("Gen1") %>%
      pivot_longer(-Gen1, names_to = "Gen2", values_to = "p_adj") %>%
      mutate(
        Treatment = pairwise_wilcox$Treatment[i],
        Metric = pairwise_wilcox$Metric[i]
      )
  })),
  "results/fig3_alpha_pairwise_wilcox.csv",
  row.names = FALSE
)

# Figure 4 — Beta diversity
write.csv(
  as.data.frame(adonis_terms),
  "results/fig4_G0G3_PERMANOVA_terms.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(adonis_margin),
  "results/fig4_G0G3_PERMANOVA_margin.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(adonis_G1G3),
  "results/fig4_G1G3_PERMANOVA_BrayCurtis.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(disp_anova),
  "results/fig4_G0G3_PERMDISP_generation.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(anova(disp_trt)),
  "results/fig4_G0G3_PERMDISP_treatment.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(anova(disp_G1G3_gen)),
  "results/fig4_G1G3_PERMDISP_generation.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(anova(disp_G1G3_trt)),
  "results/fig4_G1G3_PERMDISP_treatment.csv",
  row.names = TRUE
)

write.csv(
  do.call(rbind, lapply(names(simper_summary), function(x) {
    df <- simper_summary[[x]]
    df$comparison <- x
    df
  })),
  "results/SIMPER_by_generation.csv",
  row.names = FALSE
)

capture.output(
  disp_permutest,
  file = "results/fig4_G0G3_PERMDISP_generation_permutest.txt"
)

capture.output(
  disp_trt_permutest,
  file = "results/fig4_G0G3_PERMDISP_treatment_permutest.txt"
)

capture.output(
  permutest_G1G3_gen,
  file = "results/fig4_G1G3_PERMDISP_generation_permutest.txt"
)

capture.output(
  permutest_G1G3_trt,
  file = "results/fig4_G1G3_PERMDISP_treatment_permutest.txt"
)                                        
                                        
# Figure 5 — Differential abundance
write.csv(
  merged_df,
  "results/fig5_deseq_diet_shift.csv",
  row.names = FALSE
)

write.csv(
  anbc_merged,
  "results/fig5_ancombc_diet_shift.csv",
  row.names = FALSE
)

write.csv(
  merged_df %>%
    select(
      taxon,
      OrderName,
      oat_arrow,
      wheat_arrow,
      oat_concordant,
      wheat_concordant
    ),
  "results/fig5_concordance_summary.csv",
  row.names = FALSE
)
                                        
# Supplementary analyses
write.csv(
  overlap_summary,
  "results/G0_downstream_order_overlap_summary.csv",
  row.names = FALSE
)

write.csv(
  g1_flour_deseq,
  "results/da_G1_flourtype_deseq.csv",
  row.names = FALSE
)

write.csv(
  g1_flour_ancombc,
  "results/da_G1_flourtype_ancombc.csv",
  row.names = FALSE
)

write.csv(
  g1g2_merged,
  "results/da_G1G2_transition.csv",
  row.names = FALSE
)
