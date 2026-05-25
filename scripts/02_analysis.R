# Assumes 00_setup.R and 01_processing.R have already been sourced

# -------------------------
# Figure 2 baseline stats
# -------------------------
alpha_G0 <- estimate_richness(ps_AdultGut_G0, measures = c("Observed", "Shannon")) %>%
  rownames_to_column("SampleID") %>%
  left_join(
    data.frame(sample_data(ps_AdultGut_G0)) %>% rownames_to_column("SampleID"),
    by = "SampleID"
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("ControlA", "ControlB", "ExpA", "ExpB"))
  )

alpha_G0_long <- alpha_G0 %>%
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

kw_results_fig2 <- alpha_G0_long %>%
  group_by(Metric) %>%
  summarise(
    p = kruskal.test(Value ~ Treatment)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = paste0("Kruskal-Wallis p = ", signif(p, 3))
  )

# -------------------------
# Figure 4 alpha diversity stats
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

kw_results_fig4 <- alpha_long %>%
  group_by(Treatment, Metric) %>%
  summarise(
    p = kruskal.test(Value ~ Generation)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = paste0("Kruskal-Wallis p = ", signif(p, 3))
  )

# -------------------------
# Beta diversity and PERMANOVA
# -------------------------
bray_dist <- phyloseq::distance(ps_AdultGut, method = "bray")
meta <- data.frame(sample_data(ps_AdultGut))

adonis_terms <- adonis2(
  bray_dist ~ Generation + Treatment + Generation:Treatment,
  data = meta,
  permutations = 999,
  by = "terms"
)

adonis_margin <- adonis2(
  bray_dist ~ Generation * Treatment,
  data = meta,
  permutations = 999,
  by = "margin"
)

disp <- betadisper(bray_dist, meta$Generation)
disp_anova <- anova(disp)
disp_permutest <- permutest(disp, pairwise = TRUE)

# -------------------------
# Figure 6 DESeq2 diet-shift analysis
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

oat_G2 <- run_pairwise_deseq(ps_order, "G2", "ExpA", "ControlA", "lfc_oat_G2", "padj_oat_G2")
oat_G3 <- run_pairwise_deseq(ps_order, "G3", "ExpA", "ControlA", "lfc_oat_G3", "padj_oat_G3")
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
    oat_arrow = case_when(
      oat_same_direction & oat_any_sig & mean(c(lfc_oat_G2, lfc_oat_G3), na.rm = TRUE) > 0 ~ "↑O",
      oat_same_direction & oat_any_sig & mean(c(lfc_oat_G2, lfc_oat_G3), na.rm = TRUE) < 0 ~ "↓O",
      TRUE ~ ""
    ),
    wheat_same_direction = !is.na(lfc_wheat_G2) & !is.na(lfc_wheat_G3) &
      sign(lfc_wheat_G2) == sign(lfc_wheat_G3) &
      sign(lfc_wheat_G2) != 0,
    wheat_any_sig = any(c(padj_wheat_G2, padj_wheat_G3) < alpha, na.rm = TRUE),
    wheat_arrow = case_when(
      wheat_same_direction & wheat_any_sig & mean(c(lfc_wheat_G2, lfc_wheat_G3), na.rm = TRUE) > 0 ~ "↑W",
      wheat_same_direction & wheat_any_sig & mean(c(lfc_wheat_G2, lfc_wheat_G3), na.rm = TRUE) < 0 ~ "↓W",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

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
heatmap_mat <- heatmap_mat[, desired_col_order, drop = FALSE]

heatmap_mat <- heatmap_mat[apply(heatmap_mat, 1, var) > 0, , drop = FALSE]
heatmap_mat_filtered <- heatmap_mat[rowMeans(heatmap_mat) > 0.001, , drop = FALSE]

labels_row <- rownames(heatmap_mat_filtered) %>%
  sapply(function(ord) {
    row_i <- merged_df %>% filter(OrderName == ord)
    oat_lab <- if (nrow(row_i) > 0) unique(row_i$oat_arrow[row_i$oat_arrow != ""]) else character(0)
    wheat_lab <- if (nrow(row_i) > 0) unique(row_i$wheat_arrow[row_i$wheat_arrow != ""]) else character(0)
    arrows <- c(oat_lab, wheat_lab)
    arrows <- arrows[arrows != ""]
    if (length(arrows) > 0) {
      paste(ord, paste(arrows, collapse = " "))
    } else {
      ord
    }
  })

col_labels <- gsub("_", " ", colnames(heatmap_mat_filtered))

# -------------------------
# Save numeric outputs
# -------------------------
dir.create("results", showWarnings = FALSE)

write.csv(as.data.frame(adonis_terms), "results/permanova_terms.csv", row.names = TRUE)
write.csv(as.data.frame(adonis_margin), "results/permanova_margin.csv", row.names = TRUE)
write.csv(as.data.frame(disp_anova), "results/betadisper_anova.csv", row.names = TRUE)
write.csv(as.data.frame(kw_results_fig2), "results/fig2_kruskal.csv", row.names = FALSE)
write.csv(as.data.frame(kw_results_fig4), "results/fig4_kruskal.csv", row.names = FALSE)
write.csv(as.data.frame(merged_df), "results/fig6_deseq_diet_shift.csv", row.names = FALSE)

# -------------------------
# Reviewer analysis: Were downstream taxa present in G0?
# -------------------------

df_overlap <- psmelt(ps_order) %>%
  group_by(Order, Generation) %>%
  summarise(
    total_abundance = sum(Abundance),
    n_samples_present = sum(Abundance > 0),
    .groups = "drop"
  ) %>%
  mutate(present = total_abundance > 0)

g0_orders <- df_overlap %>%
  filter(Generation == "G0", present) %>%
  pull(Order)

downstream_orders <- df_overlap %>%
  filter(Generation %in% c("G1", "G2", "G3"), present) %>%
  pull(Order) %>%
  unique()

overlap_summary <- tibble(
  n_G0_orders = length(g0_orders),
  n_downstream_orders = length(downstream_orders),
  n_downstream_also_in_G0 = sum(downstream_orders %in% g0_orders),
  proportion_downstream_also_in_G0 = n_downstream_also_in_G0 / n_downstream_orders
)

write.csv(
  overlap_summary,
  "results/reviewer_G0_downstream_order_overlap_summary.csv",
  row.names = FALSE
)

# -------------------------
# Reviewer analysis: Key order-level relative abundance shifts
# -------------------------

order_shift_summary <- psmelt(ps_order) %>%
  group_by(Sample, Generation, Order) %>%
  summarise(Count = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(RelAbund = Count / sum(Count)) %>%
  ungroup() %>%
  filter(Order %in% c(
    "Enterobacterales_A_737866",
    "Lactobacillales",
    "Bacteroidales",
    "Burkholderiales_592522",
    "Pseudomonadales_650611",
    "Pseudomonadales_660879"
  )) %>%
  group_by(Generation, Order) %>%
  summarise(
    mean_rel_abund = mean(RelAbund, na.rm = TRUE),
    sd_rel_abund = sd(RelAbund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Order, Generation)

write.csv(
  order_shift_summary,
  "results/reviewer_key_order_relative_abundance_shifts.csv",
  row.names = FALSE
)

# -------------------------
# Reviewer analysis: Full beta diversity statistics for G0-G3
# -------------------------

bray_all <- phyloseq::distance(ps_AdultGut, method = "bray")
meta_all <- data.frame(sample_data(ps_AdultGut))

set.seed(123)
adonis_all <- vegan::adonis2(
  bray_all ~ Generation * Treatment,
  data = meta_all,
  permutations = 999
)

disp_gen_all <- betadisper(bray_all, meta_all$Generation)
disp_trt_all <- betadisper(bray_all, meta_all$Treatment)

write.csv(
  as.data.frame(adonis_all),
  "results/fig5_G0G3_PERMANOVA_BrayCurtis.csv",
  row.names = TRUE
)

capture.output(
  list(
    "PERMANOVA G0-G3" = adonis_all,
    "Generation dispersion ANOVA" = anova(disp_gen_all),
    "Generation dispersion permutation test" = permutest(disp_gen_all, permutations = 999),
    "Treatment dispersion ANOVA" = anova(disp_trt_all),
    "Treatment dispersion permutation test" = permutest(disp_trt_all, permutations = 999)
  ),
  file = "results/fig5_G0G3_beta_diversity_stats.txt"
)

# -------------------------
# Reviewer analysis: Recovery phase beta diversity only, G1-G3
# -------------------------

ps_AdultGut_G1G3 <- subset_samples(
  ps_AdultGut,
  !is.na(Generation) & Generation %in% c("G1", "G2", "G3")
)

ps_AdultGut_G1G3 <- prune_samples(sample_sums(ps_AdultGut_G1G3) > 0, ps_AdultGut_G1G3)
ps_AdultGut_G1G3 <- prune_taxa(taxa_sums(ps_AdultGut_G1G3) > 0, ps_AdultGut_G1G3)

sample_data(ps_AdultGut_G1G3)$Generation <- droplevels(
  factor(sample_data(ps_AdultGut_G1G3)$Generation, levels = c("G1", "G2", "G3"))
)

sample_data(ps_AdultGut_G1G3)$Treatment <- droplevels(
  factor(sample_data(ps_AdultGut_G1G3)$Treatment, levels = trt_levels)
)

bray_G1G3 <- phyloseq::distance(ps_AdultGut_G1G3, method = "bray")
meta_G1G3 <- data.frame(sample_data(ps_AdultGut_G1G3))
meta_G1G3 <- meta_G1G3[labels(bray_G1G3), ]

set.seed(123)
adonis_G1G3 <- vegan::adonis2(
  bray_G1G3 ~ Generation * Treatment,
  data = meta_G1G3,
  permutations = 999
)

disp_gen_G1G3 <- betadisper(bray_G1G3, meta_G1G3$Generation)
disp_trt_G1G3 <- betadisper(bray_G1G3, meta_G1G3$Treatment)

write.csv(
  as.data.frame(adonis_G1G3),
  "results/reviewer_G1G3_PERMANOVA_BrayCurtis.csv",
  row.names = TRUE
)

capture.output(
  list(
    "Sample counts G1-G3" = table(meta_G1G3$Generation, meta_G1G3$Treatment),
    "PERMANOVA G1-G3" = adonis_G1G3,
    "Generation dispersion ANOVA" = anova(disp_gen_G1G3),
    "Generation dispersion permutation test" = permutest(disp_gen_G1G3, permutations = 999),
    "Treatment dispersion ANOVA" = anova(disp_trt_G1G3),
    "Treatment dispersion permutation test" = permutest(disp_trt_G1G3, permutations = 999)
  ),
  file = "results/reviewer_G1G3_beta_diversity_stats.txt"
)

# -------------------------
# Reviewer analysis: SIMPER taxa contributing to Bray-Curtis separation
# -------------------------

extract_simper_top <- function(simper_obj, top_n = 10) {
  out <- lapply(names(simper_obj), function(comp) {

    x <- as.data.frame(simper_obj[[comp]])
    x$Order <- rownames(x)
    x$Comparison <- comp

    if ("cusum" %in% colnames(x) && !"cumsum" %in% colnames(x)) {
      x$cumsum <- x$cusum
    }

    keep_cols <- intersect(
      c("Comparison", "Order", "average", "sd", "ratio", "ava", "avb", "cumsum", "p"),
      colnames(x)
    )

    x %>%
      arrange(desc(average)) %>%
      select(all_of(keep_cols)) %>%
      slice_head(n = top_n)
  })

  bind_rows(out)
}

ps_order_rel_all <- transform_sample_counts(ps_order, function(x) x / sum(x))

order_mat_all <- t(as(otu_table(ps_order_rel_all), "matrix"))
meta_order_all <- data.frame(sample_data(ps_order_rel_all))

set.seed(123)
simper_generation <- vegan::simper(
  order_mat_all,
  group = meta_order_all$Generation,
  permutations = 999
)

ps_order_G1G3 <- subset_samples(ps_order_rel_all, Generation %in% c("G1", "G2", "G3"))
ps_order_G1G3 <- prune_taxa(taxa_sums(ps_order_G1G3) > 0, ps_order_G1G3)

order_mat_G1G3 <- t(as(otu_table(ps_order_G1G3), "matrix"))
meta_order_G1G3 <- data.frame(sample_data(ps_order_G1G3))

set.seed(123)
simper_treatment_G1G3 <- vegan::simper(
  order_mat_G1G3,
  group = meta_order_G1G3$Treatment,
  permutations = 999
)

simper_generation_top <- extract_simper_top(simper_generation, top_n = 10)
simper_treatment_G1G3_top <- extract_simper_top(simper_treatment_G1G3, top_n = 10)

write.csv(
  simper_generation_top,
  "results/reviewer_SIMPER_top_orders_by_generation.csv",
  row.names = FALSE
)

write.csv(
  simper_treatment_G1G3_top,
  "results/reviewer_SIMPER_top_orders_G1G3_by_treatment.csv",
  row.names = FALSE
)

# -------------------------
# Reviewer analysis: Alpha diversity recovery phase only, G1-G3
# -------------------------

alpha_G1G3 <- alpha_long %>%
  filter(Generation %in% c("G1", "G2", "G3"))

kw_alpha_G1G3 <- alpha_G1G3 %>%
  group_by(Treatment, Metric) %>%
  summarise(
    p = kruskal.test(Value ~ Generation)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p, method = "BH")
  )

write.csv(
  kw_alpha_G1G3,
  "results/reviewer_G1G3_alpha_Kruskal_results.csv",
  row.names = FALSE
)

pairwise_alpha_G1G3_clean <- data.frame()

for (trt in levels(alpha_G1G3$Treatment)) {
  for (met in levels(alpha_G1G3$Metric)) {

    dat <- alpha_G1G3 %>%
      filter(Treatment == trt, Metric == met)

    pw <- pairwise.wilcox.test(
      x = dat$Value,
      g = dat$Generation,
      p.adjust.method = "BH"
    )$p.value

    pw_df <- as.data.frame(as.table(pw), stringsAsFactors = FALSE)
    colnames(pw_df) <- c("Generation_1", "Generation_2", "p_adj_BH")

    pw_df <- pw_df %>%
      filter(!is.na(p_adj_BH)) %>%
      mutate(Treatment = trt, Metric = met) %>%
      select(Treatment, Metric, Generation_1, Generation_2, p_adj_BH)

    pairwise_alpha_G1G3_clean <- bind_rows(pairwise_alpha_G1G3_clean, pw_df)
  }
}

write.csv(
  pairwise_alpha_G1G3_clean,
  "results/reviewer_G1G3_alpha_pairwise_wilcox_clean.csv",
  row.names = FALSE
)




