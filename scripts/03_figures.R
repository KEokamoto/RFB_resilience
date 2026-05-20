# Assumes 00_setup.R, 01_processing.R, and 02_analysis.R have already been sourced

dir.create("figures", showWarnings = FALSE)

# -------------------------
# Figure 2
# -------------------------
p_fig2 <- ggplot(alpha_G0_long, aes(x = Treatment, y = Value)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.7, size = 1.4) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  geom_text(
    data = kw_results_fig2,
    aes(x = 1, y = Inf, label = p_label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.1,
    size = 3
  ) +
  labs(
    x = "Downstream treatment assignment",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 8)),
    panel.spacing = unit(1.2, "lines")
  ) +
  coord_cartesian(clip = "off")

ggsave(
  "figures/Fig2_G0_baseline_alpha_diversity.png",
  plot = p_fig2,
  width = 6.85,
  height = 3.8,
  dpi = 300,
  bg = "white"
)
#====FIG 3: Order-level composition in flour and adult gut====================
# Flour = Fig. 3a
# Adult gut = Fig. 3b
#---1. Prepare flour dataset--------------------------------------------------
flour_types <- c("WholeWheatFresh", "WholeWheatUsed", "OatFresh", "OatUsed")

meta_df <- data.frame(sample_data(ps))
keep_samps <- rownames(meta_df)[meta_df$Description %in% flour_types]
ps_flour <- prune_samples(keep_samps, ps)
ps_flour_order <- tax_glom(ps_flour, taxrank = "Order")
df_flour_raw <- psmelt(ps_flour_order)

#---2. Prepare adult gut dataset---------------------------------------------
ps_AdultGut_Order <- tax_glom(ps_AdultGut, taxrank = "Order")
df_gut_raw <- psmelt(ps_AdultGut_Order)

#---3. Build shared order ranking across BOTH datasets------------------------
#    using per-sample relative abundance
calc_sample_relative <- function(df) {
  df %>%
    group_by(Sample, Order) %>%
    summarise(Count = sum(Abundance), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(
      Total = sum(Count),
      Relative = 100 * Count / Total
    ) %>%
    ungroup()
}

df_flour_sample <- calc_sample_relative(df_flour_raw)
df_gut_sample   <- calc_sample_relative(df_gut_raw)

overall_abundance_all <- bind_rows(
  df_flour_sample %>% select(Order, Relative),
  df_gut_sample %>% select(Order, Relative)
) %>%
  group_by(Order) %>%
  summarise(mean_abundance = mean(Relative, na.rm = TRUE), .groups = "drop")

# Collapse low-abundance orders consistently across BOTH panels
low_abundance_orders <- overall_abundance_all %>%
  filter(mean_abundance < 0.5) %>%
  pull(Order)

df_flour_raw <- df_flour_raw %>%
  mutate(Order = if_else(Order %in% low_abundance_orders, "<0.5%", as.character(Order)))

df_gut_raw <- df_gut_raw %>%
  mutate(Order = if_else(Order %in% low_abundance_orders, "<0.5%", as.character(Order)))

# Determine shared Order levels based on pooled mean relative abundance
df_all_collapsed <- bind_rows(
  df_flour_raw %>% mutate(Source = "Flour"),
  df_gut_raw   %>% mutate(Source = "AdultGut")
)

df_all_sample <- df_all_collapsed %>%
  group_by(Source, Sample, Order) %>%
  summarise(Count = sum(Abundance), .groups = "drop") %>%
  group_by(Source, Sample) %>%
  mutate(
    Total = sum(Count),
    Relative = 100 * Count / Total
  ) %>%
  ungroup()

order_levels <- df_all_sample %>%
  group_by(Order) %>%
  summarise(mean_rel = mean(Relative, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_rel)) %>%
  pull(Order)

# Put <0.5% last, if present
if ("<0.5%" %in% order_levels) {
  order_levels <- c(setdiff(order_levels, "<0.5%"), "<0.5%")
}

# Shared color palette across both panels
n_colors <- length(order_levels)
shared_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
names(shared_colors) <- order_levels

shared_color_scale <- scale_fill_manual(values = shared_colors, drop = FALSE)

#---4. Flour panel data (Fig. 3a)---------------------------------------------
df_counts_flour <- df_flour_raw %>%
  group_by(Sample, Description, Order) %>%
  summarise(Count = sum(Abundance), .groups = "drop")

group_flour <- df_counts_flour %>%
  group_by(Description, Order) %>%
  summarise(total_count = sum(Count), .groups = "drop") %>%
  complete(
    Description = flour_types,
    Order = order_levels,
    fill = list(total_count = 0)
  )

group_flour$Description <- factor(
  group_flour$Description,
  levels = flour_types,
  labels = c("Whole wheat\nfresh", "Whole wheat\nused", "Oat\nfresh", "Oat\nused")
)

group_flour$Order <- factor(group_flour$Order, levels = order_levels)

#---5. Adult gut panel data (Fig. 3b)------------------------------------------
df_counts_gut <- df_gut_raw %>%
  group_by(Sample, Treatment, Generation, Order) %>%
  summarise(Count = sum(Abundance), .groups = "drop")

group_gut <- df_counts_gut %>%
  group_by(Treatment, Generation, Order) %>%
  summarise(total_count = sum(Count), .groups = "drop") %>%
  complete(
    Treatment = c("ControlA", "ExpA", "ControlB", "ExpB"),
    Generation = c("G0", "G1", "G2", "G3"),
    Order = order_levels,
    fill = list(total_count = 0)
  )

group_gut$Order <- factor(group_gut$Order, levels = order_levels)
group_gut$Treatment <- factor(
  group_gut$Treatment,
  levels = c("ControlA", "ExpA", "ControlB", "ExpB")
)
group_gut$Generation <- factor(
  group_gut$Generation,
  levels = c("G0", "G1", "G2", "G3")
)

#---6. Plot Fig. 3a (Flour)---------------------------------------------------
fig3a_flour <- ggplot(group_flour, aes(x = Description, y = total_count, fill = Order)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  shared_color_scale +
  labs(
    x = "Flour substrate",
    y = "Relative abundance (%)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "right"
  )

fig3a_flour

#---7. Plot Fig. 3b (Adult gut)---------------------------------------------
fig3b_gut <- ggplot(group_gut, aes(x = Treatment, y = total_count, fill = Order)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  facet_wrap(~Generation, nrow = 1) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  shared_color_scale +
  labs(
    x = "Treatment",
    y = "Relative abundance (%)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 8)),
    panel.spacing = unit(1, "lines"),
    legend.position = "right"
  )

fig3b_gut

#---8. Save separately-------------------------------------------------------
# Fig. 3a: flour panel
ggsave(
  "Fig3a_flour_order_composition.png",
  plot = fig3a_flour,
  width = 6.85,
  height = 3.4,
  dpi = 300,
  bg = "white"
)

# Fig. 3b: adult gut panel
ggsave(
  "Fig3b_adultgut_order_composition.png",
  plot = fig3b_gut,
  width = 6.85,
  height = 3.6,
  dpi = 300,
  bg = "white"
)

#=== Reviewer-requested overlap: Were downstream taxa present in G0? ===

# Use order-level adult gut data
df_overlap <- psmelt(ps_AdultGut_Order) %>%
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

print(overlap_summary)

write.csv(
  overlap_summary,
  "Reviewer_G0_downstream_order_overlap_summary.csv",
  row.names = FALSE
)

#=== Quantify order-level relative abundance shifts across generations ===

order_shift_summary <- psmelt(ps_AdultGut_Order) %>%
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

print(order_shift_summary)

write.csv(
  order_shift_summary,
  "Reviewer_key_order_relative_abundance_shifts.csv",
  row.names = FALSE
)

#===FIG 4: Alpha diversity dynamics across generations========================
# Observed richness (top) + Shannon diversity (bottom)

# Set factor levels
sample_data(ps_AdultGut)$Treatment <- factor(
  sample_data(ps_AdultGut)$Treatment,
  levels = trt_levels
)

sample_data(ps_AdultGut)$Generation <- factor(
  sample_data(ps_AdultGut)$Generation,
  levels = gen_levels
)

# Build alpha-diversity dataframe
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

# Long format
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

# Kruskal-Wallis p-value within each Treatment x Metric
kw_results <- alpha_long %>%
  group_by(Treatment, Metric) %>%
  summarise(
    p = kruskal.test(Value ~ Generation)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = paste0("Kruskal-Wallis p = ", signif(p, 3))
  )

# Shared theme
fig4_theme <- theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 8)),
    panel.spacing = unit(1, "lines"),
    legend.position = "right",
    plot.margin = margin(6, 6, 6, 10)
  )

# Observed richness panel
p_observed <- alpha_long %>%
  filter(Metric == "Observed richness") %>%
  ggplot(aes(x = Generation, y = Value)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_jitter(
    aes(color = Generation, shape = Treatment),
    width = 0.12, height = 0,
    alpha = 0.7, size = 1.6
  ) +
  facet_wrap(~Treatment, nrow = 1) +
  scale_color_manual(values = gen_cols, drop = FALSE, name = "Generation") +
  scale_shape_manual(values = trt_shapes, drop = FALSE, name = "Treatment") +
  geom_text(
    data = kw_results %>% filter(Metric == "Observed richness"),
    aes(x = 1, y = Inf, label = p_label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.1,
    size = 2.9
  ) +
  labs(
    x = "Generation",
    y = "Observed richness"
  ) +
  coord_cartesian(clip = "off") +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  fig4_theme

# Shannon diversity panel
p_shannon <- alpha_long %>%
  filter(Metric == "Shannon diversity") %>%
  ggplot(aes(x = Generation, y = Value)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_jitter(
    aes(color = Generation, shape = Treatment),
    width = 0.12, height = 0,
    alpha = 0.7, size = 1.6
  ) +
  facet_wrap(~Treatment, nrow = 1) +
  scale_color_manual(values = gen_cols, drop = FALSE, name = "Generation") +
  scale_shape_manual(values = trt_shapes, drop = FALSE, name = "Treatment") +
  geom_text(
    data = kw_results %>% filter(Metric == "Shannon diversity"),
    aes(x = 1, y = Inf, label = p_label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.1,
    size = 2.9
  ) +
  labs(
    x = "Generation",
    y = "Shannon diversity"
  ) +
  coord_cartesian(clip = "off") +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  fig4_theme

# Combine with one shared legend
combined_plot <- p_observed / p_shannon +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(combined_plot)

# Save 
ggsave(
  "Fig4_alpha_diversity_dynamics.png",
  plot = combined_plot,
  width = 6.85,
  height = 6.4,
  dpi = 300,
  bg = "white"
)

#=== Figure 5. Beta diversity of adult gut microbiomes (Bray-Curtis PCoA) ===
# Ensure factors/levels
sample_data(ps_AdultGut)$Treatment <- factor(
  sample_data(ps_AdultGut)$Treatment,
  levels = trt_levels
)

sample_data(ps_AdultGut)$Generation <- factor(
  sample_data(ps_AdultGut)$Generation,
  levels = gen_levels
)

# Ordination
pcoa_bray <- ordinate(ps_AdultGut, method = "PCoA", distance = "bray")

# Variance explained
pc1_var <- round(100 * pcoa_bray$values$Relative_eig[1], 1)
pc2_var <- round(100 * pcoa_bray$values$Relative_eig[2], 1)

# Plot
pcoa_plot1 <- plot_ordination(ps_AdultGut, pcoa_bray, color = "Generation") +
  geom_point(aes(shape = Treatment), size = 2.6, alpha = 0.9) +
  stat_ellipse(
    aes(group = Generation),
    type = "t",
    linetype = "dashed",
    linewidth = 0.4,
    show.legend = FALSE
  ) +
  scale_shape_manual(
    values = trt_shapes,
    drop = FALSE,
    name = "Treatment"
  ) +
  scale_color_manual(
    values = gen_cols,
    drop = FALSE,
    name = "Generation"
  ) +
  labs(
    x = paste0("PCoA 1 (", pc1_var, "%)"),
    y = paste0("PCoA 2 (", pc2_var, "%)")
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "plain"),
    axis.text = element_text(face = "plain"),
    legend.title = element_text(face = "plain"),
    legend.text = element_text(face = "plain"),
    plot.margin = margin(8, 8, 8, 8)
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

print(pcoa_plot1)

ggsave(
  "Fig5_PCoA_AdultGut.png",
  plot = pcoa_plot1,
  width = 6.85,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

#=== Full beta diversity statistics for Figure 5: G0-G3 ===

bray_all <- phyloseq::distance(ps_AdultGut, method = "bray")
meta_all <- data.frame(sample_data(ps_AdultGut))

set.seed(123)
adonis_all <- vegan::adonis2(
  bray_all ~ Generation * Treatment,
  data = meta_all,
  permutations = 999
)

print(adonis_all)

write.csv(
  as.data.frame(adonis_all),
  "Fig5_G0G3_PERMANOVA_BrayCurtis.csv"
)

disp_gen_all <- betadisper(bray_all, meta_all$Generation)
disp_trt_all <- betadisper(bray_all, meta_all$Treatment)

capture.output(
  list(
    "PERMANOVA G0-G3" = adonis_all,
    "Generation dispersion ANOVA" = anova(disp_gen_all),
    "Generation dispersion permutation test" = permutest(disp_gen_all, permutations = 999),
    "Treatment dispersion ANOVA" = anova(disp_trt_all),
    "Treatment dispersion permutation test" = permutest(disp_trt_all, permutations = 999)
  ),
  file = "Fig5_G0G3_beta_diversity_stats.txt"
)

permutest(disp_gen_all, pairwise = TRUE, permutations = 999)

#=== Reviewer-requested analysis: Beta diversity during recovery phase only ===
# Goal: separate the initial G0 -> G1 perturbation from post-perturbation dynamics G1-G3
ps_AdultGut_G1G3 <- subset_samples(ps_AdultGut, !is.na(Generation) & Generation %in% c("G1", "G2", "G3"))
ps_AdultGut_G1G3 <- prune_samples(sample_sums(ps_AdultGut_G1G3) > 0, ps_AdultGut_G1G3)
ps_AdultGut_G1G3 <- prune_taxa(taxa_sums(ps_AdultGut_G1G3) > 0, ps_AdultGut_G1G3)

sample_data(ps_AdultGut_G1G3)$Generation <- droplevels(
  factor(sample_data(ps_AdultGut_G1G3)$Generation, levels = c("G1", "G2", "G3"))
)

sample_data(ps_AdultGut_G1G3)$Treatment <- droplevels(
  factor(sample_data(ps_AdultGut_G1G3)$Treatment, levels = trt_levels)
)

# Check sample numbers before analysis
print(table(sample_data(ps_AdultGut_G1G3)$Generation,
            sample_data(ps_AdultGut_G1G3)$Treatment))

bray_G1G3 <- phyloseq::distance(ps_AdultGut_G1G3, method = "bray")
meta_G1G3 <- data.frame(sample_data(ps_AdultGut_G1G3))

# Make sure metadata rownames match distance labels
meta_G1G3 <- meta_G1G3[labels(bray_G1G3), ]

set.seed(123)
adonis_G1G3 <- vegan::adonis2(
  bray_G1G3 ~ Generation * Treatment,
  data = meta_G1G3,
  permutations = 999
)

print(adonis_G1G3)
write.csv(as.data.frame(adonis_G1G3), "Reviewer_G1G3_PERMANOVA_BrayCurtis.csv")

# Dispersion tests
disp_gen_G1G3 <- betadisper(bray_G1G3, meta_G1G3$Generation)
disp_trt_G1G3 <- betadisper(bray_G1G3, meta_G1G3$Treatment)

capture.output(
  list(
    "PERMANOVA G1-G3" = adonis_G1G3,
    "Generation dispersion ANOVA" = anova(disp_gen_G1G3),
    "Generation dispersion permutation test" = permutest(disp_gen_G1G3, permutations = 999),
    "Treatment dispersion ANOVA" = anova(disp_trt_G1G3),
    "Treatment dispersion permutation test" = permutest(disp_trt_G1G3, permutations = 999)
  ),
  file = "Reviewer_G1G3_beta_diversity_stats.txt"
)

# Ordination
pcoa_bray_G1G3 <- ordinate(ps_AdultGut_G1G3, method = "PCoA", distance = "bray")

pc1_var_G1G3 <- round(100 * pcoa_bray_G1G3$values$Relative_eig[1], 1)
pc2_var_G1G3 <- round(100 * pcoa_bray_G1G3$values$Relative_eig[2], 1)

# Safer plot: no ellipses, because some groups may have too few samples
pcoa_plot_G1G3 <- plot_ordination(ps_AdultGut_G1G3, pcoa_bray_G1G3, color = "Generation") +
  geom_point(aes(shape = Treatment), size = 2.6, alpha = 0.9) +
  scale_shape_manual(values = trt_shapes, drop = FALSE, name = "Treatment") +
  scale_color_manual(
    values = gen_cols[c("G1", "G2", "G3")],
    drop = FALSE,
    name = "Generation"
  ) +
  labs(
    x = paste0("PCoA 1 (", pc1_var_G1G3, "%)"),
    y = paste0("PCoA 2 (", pc2_var_G1G3, "%)")
  ) +
  theme_classic(base_size = 10)

print(pcoa_plot_G1G3)

ggsave(
  "Reviewer_Fig_G1G3_PCoA_AdultGut.png",
  plot = pcoa_plot_G1G3,
  width = 6.85,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

#=== Reviewer-requested taxa contributing to PCoA separation ===
# SIMPER identifies taxa contributing most to Bray-Curtis dissimilarity.
# This directly complements the Bray-Curtis PCoA.

# Use order-level data for interpretability
extract_simper_top <- function(simper_obj, top_n = 10) {
  out <- lapply(names(simper_obj), function(comp) {
    
    x <- as.data.frame(simper_obj[[comp]])
    x$Order <- rownames(x)
    x$Comparison <- comp
    
    # Some vegan versions use "cusum" instead of "cumsum"
    if ("cusum" %in% colnames(x) && !"cumsum" %in% colnames(x)) {
      x$cumsum <- x$cusum
    }
    
    # Keep only columns that exist
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

simper_generation_top <- extract_simper_top(simper_generation, top_n = 10)
simper_treatment_G1G3_top <- extract_simper_top(simper_treatment_G1G3, top_n = 10)

print(simper_generation_top)
print(simper_treatment_G1G3_top)

write.csv(
  simper_generation_top,
  "Reviewer_SIMPER_top_orders_by_generation.csv",
  row.names = FALSE
)

write.csv(
  simper_treatment_G1G3_top,
  "Reviewer_SIMPER_top_orders_G1G3_by_treatment.csv",
  row.names = FALSE
)

#=== Figure 6. Diet-shift differential abundance + order-level heatmap ===
#---Step A: Order-glommed object------------------------------------------
ps_order <- tax_glom(ps_AdultGut, taxrank = "Order")

sample_data(ps_order)$Treatment <- factor(
  sample_data(ps_order)$Treatment,
  levels = c("ControlA", "ExpA", "ControlB", "ExpB")
)

sample_data(ps_order)$Generation <- factor(
  sample_data(ps_order)$Generation,
  levels = c("G0", "G1", "G2", "G3")
)

# taxon -> Order name map
tax_map <- as.data.frame(tax_table(ps_order)) %>%
  rownames_to_column("taxon") %>%
  select(taxon, OrderName = Order) %>%
  mutate(
    OrderName = gsub("_[0-9]+$", "", OrderName),
    OrderName = ifelse(OrderName == "" | is.na(OrderName), "Unclassified", OrderName)
  )

#---Helper: run DESeq2 for one generation and one pair only------------------
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

#---Step B: run diet-shift contrasts within generation---------------------

# Oat-associated contrasts: ExpA vs ControlA
oat_G2 <- run_pairwise_deseq(
  ps_order, "G2", "ExpA", "ControlA",
  "lfc_oat_G2", "padj_oat_G2"
)

oat_G3 <- run_pairwise_deseq(
  ps_order, "G3", "ExpA", "ControlA",
  "lfc_oat_G3", "padj_oat_G3"
)

# Whole-wheat-associated contrasts: ExpB vs ControlB
wheat_G2 <- run_pairwise_deseq(
  ps_order, "G2", "ExpB", "ControlB",
  "lfc_wheat_G2", "padj_wheat_G2"
)

wheat_G3 <- run_pairwise_deseq(
  ps_order, "G3", "ExpB", "ControlB",
  "lfc_wheat_G3", "padj_wheat_G3"
)

#---Step C: merge results and assign arrows---------------------------------
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

#---Step D: build heatmap matrix------------------------------------------
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

# column order 
desired_col_order <- c(
  "ControlA_G0", "ControlA_G1", "ControlA_G2", "ControlA_G3",
  "ExpA_G0",     "ExpA_G1",     "ExpA_G2",     "ExpA_G3",
  "ControlB_G0", "ControlB_G1", "ControlB_G2", "ControlB_G3",
  "ExpB_G0",     "ExpB_G1",     "ExpB_G2",     "ExpB_G3"
)

desired_col_order <- desired_col_order[desired_col_order %in% colnames(heatmap_mat)]
heatmap_mat <- heatmap_mat[, desired_col_order, drop = FALSE]

# Filter rows
heatmap_mat <- heatmap_mat[apply(heatmap_mat, 1, var) > 0, , drop = FALSE]
heatmap_mat_filtered <- heatmap_mat[rowMeans(heatmap_mat) > 0.001, , drop = FALSE]

#---Step E: build row labels with diet-specific arrows---------------------
labels_row <- rownames(heatmap_mat_filtered) %>%
  sapply(function(ord) {
    
    row_i <- merged_df %>% filter(OrderName == ord)
    
    oat_lab <- if (nrow(row_i) > 0)
      unique(row_i$oat_arrow[row_i$oat_arrow != ""])
    else character(0)
    
    wheat_lab <- if (nrow(row_i) > 0)
      unique(row_i$wheat_arrow[row_i$wheat_arrow != ""])
    else character(0)
    
    arrows <- c(oat_lab, wheat_lab)
    arrows <- arrows[arrows != ""]
    
    if (length(arrows) > 0) {
      
      pretty_arrows <- arrows
      
      pretty_arrows <- gsub("↑O", "[OAT↑]", pretty_arrows)
      pretty_arrows <- gsub("↓O", "[OAT↓]", pretty_arrows)
      pretty_arrows <- gsub("↑W", "[WHEAT↑]", pretty_arrows)
      pretty_arrows <- gsub("↓W", "[WHEAT↓]", pretty_arrows)
      
      paste(ord, paste(pretty_arrows, collapse = " "))
      
    } else {
      ord
    }
  })
#---Step F: cleaner column labels------------------------------------
col_labels <- gsub("_", " ", colnames(heatmap_mat_filtered))

#---Step G: draw and save heatmap--------------------------------------------
heat_colors <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100)

png(
  filename = "Fig6_diet_shift_heatmap.png",
  width = 2055,   # 6.85 in * 300 dpi
  height = 2500,  # extra space for angled labels
  res = 300
)

pheatmap(
  mat = as.matrix(heatmap_mat_filtered),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,   # <-- dendrogram ON
  fontsize_row = 8,
  fontsize_col = 7,
  labels_row = labels_row,
  labels_col = col_labels,
  color = heat_colors,
  border_color = "grey80",
  angle_col = 90         # <-- rotated labels
)

dev.off()

#=== Figure 7. Shannon diversity trajectories across generations ============
# Define levels and palettes
gen_levels <- c("G0", "G1", "G2", "G3")
trt_levels <- c("ControlA", "ExpA", "ControlB", "ExpB")

# Keep treatment shapes consistent
trt_shapes <- c(
  ControlA = 16,
  ExpA     = 17,
  ControlB = 15,
  ExpB     = 18
)

# Treatment colors
trt_cols <- c(
  ControlA = "#840032",
  ExpA     = "#f05006",
  ControlB = "#25998f",
  ExpB     = "#f36e98"
)

# Disturbance shading
fill_vals <- c(
  "Pre-disturbance"  = "#d9d9d9",
  "Post-disturbance" = "#bdbdbd"
)

# Ensure factors
rich_combined <- rich_combined %>%
  mutate(
    Generation = factor(Generation, levels = gen_levels),
    Treatment  = factor(Treatment, levels = trt_levels),
    Gen_num    = as.numeric(Generation)
  )

# Summaries by Treatment x Generation
rich_summary <- rich_combined %>%
  group_by(Treatment, Generation, Gen_num) %>%
  summarise(
    Mean     = mean(Shannon, na.rm = TRUE),
    SD       = sd(Shannon, na.rm = TRUE),
    N        = n(),
    Lower_CI = Mean - 1.96 * SD / sqrt(N),
    Upper_CI = Mean + 1.96 * SD / sqrt(N),
    .groups  = "drop"
  )

# Global baseline (G0) mean Shannon across all commercial-stock samples
baseline_df <- rich_combined %>%
  filter(Generation == "G0") %>%
  summarise(
    y0 = mean(Shannon, na.rm = TRUE)
  )

# Disturbance periods
disturbance <- data.frame(
  Period = c("Pre-disturbance", "Post-disturbance"),
  xmin = c(0.5, 1.5),
  xmax = c(1.5, 4.5),
  ymin = -Inf,
  ymax = Inf
)

fig7 <- ggplot() +
  # Disturbance shading
  geom_rect(
    data = disturbance,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Period),
    alpha = 0.25,
    inherit.aes = FALSE
  ) +
  
  # CI ribbons by treatment (hidden from legend)
  geom_ribbon(
    data = rich_summary,
    aes(x = Gen_num, ymin = Lower_CI, ymax = Upper_CI, group = Treatment, fill = Treatment),
    alpha = 0.14,
    show.legend = FALSE
  ) +
  
  # Mean trajectories by treatment
  geom_line(
    data = rich_summary,
    aes(x = Gen_num, y = Mean, group = Treatment, color = Treatment),
    linewidth = 0.9
  ) +
  
  # Dashed global G0 baseline mean
  geom_hline(
    data = baseline_df,
    aes(yintercept = y0),
    linetype = "dashed",
    linewidth = 0.7,
    color = "black",
    show.legend = FALSE
  ) +
  
  # Individual sample points
  geom_point(
    data = rich_combined,
    aes(x = Gen_num, y = Shannon, color = Treatment, shape = Treatment),
    size = 2.0,
    alpha = 0.8
  ) +
  
  scale_x_continuous(
    breaks = 1:4,
    labels = gen_levels,
    expand = c(0.02, 0.02)
  ) +
  coord_cartesian(xlim = c(0.5, 4.5)) +
  
  # Treatment colors
  scale_color_manual(
    values = trt_cols,
    drop = FALSE,
    name = "Treatment"
  ) +
  
  # Treatment shapes
  scale_shape_manual(
    values = trt_shapes,
    drop = FALSE,
    name = "Treatment"
  ) +
  
  # Fill scale: disturbance only in legend
  scale_fill_manual(
    values = c(fill_vals, trt_cols),
    breaks = c("Pre-disturbance", "Post-disturbance"),
    name = NULL
  ) +
  
  labs(
    x = "Generation",
    y = "Shannon diversity"
  ) +
  theme_classic(base_size = 10) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 1),
    fill  = guide_legend(order = 2)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  )

print(fig7)

ggsave(
  "Fig7_shannon_diversity_trajectories.png",
  plot = fig7,
  width = 6.85,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

#=== Reviewer-requested alpha diversity analysis: recovery phase only G1-G3 ===

alpha_G1G3 <- alpha_long %>%
  filter(Generation %in% c("G1", "G2", "G3"))

# Kruskal-Wallis tests across G1-G3 within each treatment and metric
kw_alpha_G1G3 <- alpha_G1G3 %>%
  group_by(Treatment, Metric) %>%
  summarise(
    p = kruskal.test(Value ~ Generation)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p, method = "BH")
  )

print(kw_alpha_G1G3)

write.csv(
  kw_alpha_G1G3,
  "Reviewer_G1G3_alpha_Kruskal_results.csv",
  row.names = FALSE
)

# Pairwise Wilcoxon tests within each treatment and metric
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

print(pairwise_alpha_G1G3_clean)

write.csv(
  pairwise_alpha_G1G3_clean,
  "Reviewer_G1G3_alpha_pairwise_wilcox_clean.csv",
  row.names = FALSE
)
