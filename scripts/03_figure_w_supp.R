# ── figures_and_tables.R ─────────────────────────────────────────────────────
# Assumes 00_setup.R, 01_processing.R, and 02_analysis.R have already been sourced.
#
# Figure numbering (post-revision):
#   Fig 1  — experimental design diagram (not generated here)
#   Fig 2  — community composition bars: flour substrates + gut (all samples)
#   Fig 3  — alpha diversity dynamics: observed richness + Shannon boxplots
#             NOTE: Fig 3 code is in 02_analysis.R / a separate plotting script.
#             It is NOT present in this file and must be verified before submission.
#   Fig 4  — two-panel PCoA: all generations (a) + G1–G3 post-disturbance (b)
#   Fig 5  — heatmap of diet-associated order-level shifts
#   Fig 6  — Shannon diversity trajectories with global G0 baseline
#
# Supplementary Tables:
#   Supp Table 1 — SIMPER results (corrected: MD5 hash suffix fix applied)
#   Supp Table 2 — PERMANOVA and PERMDISP summary
# ─────────────────────────────────────────────────────────────────────────────

dir.create("figures", showWarnings = FALSE)

library(dplyr)
library(readr)
library(flextable)
library(officer)
library(cowplot)

# ─────────────────────────────────────────────────────────────────────────────
# ARCHIVED: G0 baseline alpha diversity figure
# Removed from main paper per Reviewer 3 (statistics reported in text,
# reference to Fig 3). File retained here for supplementary use if needed.
# ─────────────────────────────────────────────────────────────────────────────
# p_fig_baseline <- ggplot(alpha_G0_long, aes(x = Treatment, y = Value)) +
#   geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
#   geom_jitter(width = 0.12, height = 0, alpha = 0.7, size = 1.4) +
#   facet_wrap(~Metric, scales = "free_y", nrow = 1) +
#   geom_text(
#     data = kw_results_fig2,
#     aes(x = 1, y = Inf, label = p_label),
#     inherit.aes = FALSE,
#     hjust = 0, vjust = 1.1, size = 3
#   ) +
#   labs(x = "Downstream treatment assignment", y = NULL) +
#   theme_classic(base_size = 10) +
#   theme(
#     strip.background = element_blank(),
#     strip.text       = element_text(face = "plain"),
#     axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
#     axis.title.x     = element_text(margin = margin(t = 8)),
#     panel.spacing    = unit(1.2, "lines")
#   ) +
#   coord_cartesian(clip = "off")
# ggsave("figures/Fig_G0_baseline_alpha_ARCHIVED.png",
#        plot = p_fig_baseline, width = 6.85, height = 3.8, dpi = 300, bg = "white")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — community composition bars: flour substrates (a) + gut (b)
# ─────────────────────────────────────────────────────────────────────────────

flour_types <- c("WholeWheatFresh", "WholeWheatUsed", "OatFresh", "OatUsed")

meta_df_all <- data.frame(sample_data(ps))
keep_samps  <- rownames(meta_df_all)[meta_df_all$Description %in% flour_types]
ps_flour    <- prune_samples(keep_samps, ps)
ps_flour_order <- tax_glom(ps_flour, taxrank = "Order")
df_flour_raw   <- psmelt(ps_flour_order)

ps_AdultGut_Order <- tax_glom(ps_AdultGut, taxrank = "Order")
df_gut_raw        <- psmelt(ps_AdultGut_Order)

calc_sample_relative <- function(df) {
  df %>%
    group_by(Sample, Order) %>%
    summarise(Count = sum(Abundance), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(Total = sum(Count), Relative = 100 * Count / Total) %>%
    ungroup()
}

df_flour_sample <- calc_sample_relative(df_flour_raw)
df_gut_sample   <- calc_sample_relative(df_gut_raw)

overall_abundance_all <- bind_rows(
  df_flour_sample %>% select(Order, Relative),
  df_gut_sample   %>% select(Order, Relative)
) %>%
  group_by(Order) %>%
  summarise(mean_abundance = mean(Relative, na.rm = TRUE), .groups = "drop")

low_abundance_orders <- overall_abundance_all %>%
  filter(mean_abundance < 0.5) %>%
  pull(Order)

df_flour_raw <- df_flour_raw %>%
  mutate(Order = if_else(Order %in% low_abundance_orders, "<0.5%", as.character(Order)))

df_gut_raw <- df_gut_raw %>%
  mutate(Order = if_else(Order %in% low_abundance_orders, "<0.5%", as.character(Order)))

df_all_sample <- bind_rows(
  df_flour_raw %>% mutate(Source = "Flour"),
  df_gut_raw   %>% mutate(Source = "AdultGut")
) %>%
  group_by(Source, Sample, Order) %>%
  summarise(Count = sum(Abundance), .groups = "drop") %>%
  group_by(Source, Sample) %>%
  mutate(Total = sum(Count), Relative = 100 * Count / Total) %>%
  ungroup()

order_levels <- df_all_sample %>%
  group_by(Order) %>%
  summarise(mean_rel = mean(Relative, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_rel)) %>%
  pull(Order)

# ── Build plot-ready data frames for flour and gut panels ────────────────────

# Flour: compute sample-level relative abundances with Description metadata
df_flour_sample_plot <- df_flour_raw %>%
  group_by(Sample, Order) %>%
  summarise(Count = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Total = sum(Count), Relative = 100 * Count / Total) %>%
  ungroup() %>%
  left_join(
    data.frame(sample_data(ps_flour)) %>%
      rownames_to_column("Sample") %>%
      select(Sample, Description),
    by = "Sample"
  )

# Gut: compute sample-level relative abundances with Treatment and Generation metadata
df_gut_sample_plot <- df_gut_raw %>%
  group_by(Sample, Order) %>%
  summarise(Count = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Total = sum(Count), Relative = 100 * Count / Total) %>%
  ungroup() %>%
  left_join(
    data.frame(sample_data(ps_AdultGut)) %>%
      rownames_to_column("Sample") %>%
      select(Sample, Treatment, Generation),
    by = "Sample"
  )

if ("<0.5%" %in% order_levels) {
  order_levels <- c(setdiff(order_levels, "<0.5%"), "<0.5%")
}

shared_colors        <- colorRampPalette(brewer.pal(12, "Set3"))(length(order_levels))
names(shared_colors) <- order_levels
shared_color_scale   <- scale_fill_manual(values = shared_colors, drop = FALSE)

# Panel a: flour substrates — narrow bar for single-sample groups
sample_counts <- df_flour_sample_plot %>%
  distinct(Sample, Description) %>%
  group_by(Description) %>%
  summarise(n_samples = n(), .groups = "drop")

df_flour_sample_plot <- df_flour_sample_plot %>%
  left_join(sample_counts, by = "Description") %>%
  mutate(bar_width = ifelse(n_samples == 1, 0.27, 0.8))

fig2a_flour <- ggplot(df_flour_sample_plot,
                      aes(x = Sample, y = Relative, fill = Order)) +
  geom_bar(stat = "identity", position = "fill",
           aes(width = bar_width)) +
  facet_wrap(~Description, scales = "free_x", nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  shared_color_scale +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "plain", size = 9),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    legend.position  = "none"
  )

# Panel b: gut — individual sample bars, faceted by Treatment x Generation
sample_reps <- df_gut_sample_plot %>%
  distinct(Treatment, Generation, Sample) %>%
  group_by(Treatment, Generation) %>%
  mutate(rep_num = row_number()) %>%
  ungroup()

df_gut_sample_plot <- df_gut_sample_plot %>%
  left_join(sample_reps, by = c("Treatment", "Generation", "Sample"))

fig2b_gut <- ggplot(df_gut_sample_plot,
                    aes(x = factor(rep_num), y = Relative, fill = Order)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  facet_grid(rows = vars(Treatment), cols = vars(Generation)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  shared_color_scale +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "plain"),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.spacing    = unit(0.4, "lines"),
    legend.position  = "right",
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm"),
    legend.title     = element_text(size = 8)
  )

fig2_combined <- fig2a_flour / fig2b_gut +
  plot_layout(heights = c(1, 2.2)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix  = "(",
                  tag_suffix  = ")")

ggsave("figures/Fig2_composition.png",
       plot = fig2_combined,
       width = 9.0, height = 9.5, dpi = 300, bg = "white")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — alpha diversity dynamics (observed richness + Shannon boxplots)
# WARNING: this figure's plotting code is NOT present in this script.
# Verify that it exists in 02_analysis.R or a dedicated plotting script,
# and that the output is saved to figures/Fig3_alpha_diversity.png before
# submission. The manuscript cites this figure throughout the Results section.
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 — two-panel PCoA: all generations (a) + G1–G3 post-disturbance (b)
# ─────────────────────────────────────────────────────────────────────────────

# ── Build G1–G3 subset phyloseq object and ordination ────────────────────────

ps_G1G3 <- subset_samples(ps_AdultGut, Generation %in% c("G1", "G2", "G3"))

# Drop any ASVs that are now zero across all G1-G3 samples
ps_G1G3 <- prune_taxa(taxa_sums(ps_G1G3) > 0, ps_G1G3)

# Ordinate independently on the G1-G3 subset
pcoa_G1G3 <- ordinate(ps_G1G3, method = "PCoA", distance = "bray")
pc1_G1G3  <- round(100 * pcoa_G1G3$values$Relative_eig[1], 1)
pc2_G1G3  <- round(100 * pcoa_G1G3$values$Relative_eig[2], 1)

# Sanity check
cat("G1-G3 samples:", nsamples(ps_G1G3), "\n")
cat("PCoA1:", pc1_G1G3, "% PCoA2:", pc2_G1G3, "%\n")


pcoa_bray <- ordinate(ps_AdultGut, method = "PCoA", distance = "bray")
pc1_var   <- round(100 * pcoa_bray$values$Relative_eig[1], 1)
pc2_var   <- round(100 * pcoa_bray$values$Relative_eig[2], 1)

# Panel a: all generations — no legend (shared legend extracted below)
fig4a <- plot_ordination(ps_AdultGut, pcoa_bray, color = "Generation") +
  geom_point(aes(shape = Treatment), size = 2.6, alpha = 0.9) +
  stat_ellipse(aes(group = Generation), type = "t",
               linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  scale_shape_manual(values = trt_shapes, drop = FALSE, name = "Treatment") +
  scale_color_manual(values = gen_cols,   drop = FALSE, name = "Generation") +
  labs(x = paste0("PCoA 1 (", pc1_var, "%)"),
       y = paste0("PCoA 2 (", pc2_var, "%)"),
       title = "All generations (G0\u2013G3)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        plot.title       = element_text(size = 9, face = "plain"),
        plot.margin      = margin(8, 8, 8, 8)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# Panel b: G1–G3 only — legend attached here for extraction
fig4b_with_legend <- plot_ordination(ps_G1G3, pcoa_G1G3, color = "Generation") +
  geom_point(aes(shape = Treatment), size = 2.6, alpha = 0.9) +
  stat_ellipse(aes(group = Generation), type = "t",
               linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  scale_shape_manual(values = trt_shapes, drop = FALSE, name = "Treatment") +
  scale_color_manual(
    values = c("G0" = gen_cols["G0"], "G1" = gen_cols["G1"],
               "G2" = gen_cols["G2"], "G3" = gen_cols["G3"]),
    limits = c("G0", "G1", "G2", "G3"),
    drop   = FALSE,
    name   = "Generation"
  ) +
  labs(x = paste0("PCoA 1 (", pc1_G1G3, "%)"),
       y = paste0("PCoA 2 (", pc2_G1G3, "%)"),
       title = "Post-disturbance only (G1\u2013G3)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "right",
        plot.title       = element_text(size = 9, face = "plain"),
        plot.margin      = margin(8, 8, 8, 8)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

fig4b <- fig4b_with_legend + theme(legend.position = "none")

shared_legend <- get_legend(fig4b_with_legend)

fig4_combined <- plot_grid(
  plot_grid(fig4a, fig4b, nrow = 1, labels = c("(a)", "(b)"), label_size = 10),
  shared_legend,
  nrow       = 1,
  rel_widths = c(1, 0.15)
)

ggsave("figures/Fig4_PCoA_AdultGut.png",
       plot = fig4_combined, width = 13.0, height = 4.8, dpi = 300, bg = "white")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 5 — heatmap of diet-associated order-level shifts
# labels_row and heatmap_mat_filtered are computed in 02_analysis.R
# ─────────────────────────────────────────────────────────────────────────────

heat_colors <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100)

png(
  filename = "figures/Fig5_diet_shift_heatmap.png",
  width    = 3000,
  height   = 4000,
  res      = 300
)

pheatmap(
  mat          = as.matrix(heatmap_mat_filtered),
  scale        = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 7,
  fontsize_col = 7,
  labels_row   = labels_row,
  labels_col   = col_labels,
  color        = heat_colors,
  border_color = "grey85",
  angle_col    = 45,
  cellheight   = 10,
  cellwidth    = 18
)

dev.off()


# ─────────────────────────────────────────────────────────────────────────────
# Figure 6 — Shannon diversity trajectories with single global G0 baseline
# ─────────────────────────────────────────────────────────────────────────────

rich_combined <- rich_combined %>%
  mutate(
    Generation = factor(Generation, levels = gen_levels),
    Treatment  = factor(Treatment,  levels = trt_levels),
    Gen_num    = as.numeric(Generation)
  )

rich_summary <- rich_combined %>%
  group_by(Treatment, Generation, Gen_num) %>%
  summarise(
    Mean     = mean(Shannon, na.rm = TRUE),
    SD       = sd(Shannon,   na.rm = TRUE),
    N        = n(),
    Lower_CI = Mean - 1.96 * SD / sqrt(N),
    Upper_CI = Mean + 1.96 * SD / sqrt(N),
    .groups  = "drop"
  )

# Single global G0 mean — baseline did not differ among treatment assignments
# (Kruskal-Wallis p=0.444), so one line is appropriate
global_baseline <- rich_combined %>%
  filter(Generation == "G0") %>%
  summarise(y0 = mean(Shannon, na.rm = TRUE)) %>%
  pull(y0)

disturbance <- data.frame(
  Period = c("Pre-disturbance", "Post-disturbance"),
  xmin   = c(0.5, 1.5),
  xmax   = c(1.5, 4.5),
  ymin   = -Inf,
  ymax   = Inf
)

fig6 <- ggplot() +
  geom_rect(
    data = disturbance,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Period),
    alpha = 0.25, inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = rich_summary,
    aes(x = Gen_num, ymin = Lower_CI, ymax = Upper_CI,
        group = Treatment, fill = Treatment),
    alpha = 0.14, show.legend = FALSE
  ) +
  geom_line(
    data = rich_summary,
    aes(x = Gen_num, y = Mean, group = Treatment, color = Treatment),
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = global_baseline,
    linetype = "dashed", linewidth = 0.7, color = "grey40"
  ) +
  geom_point(
    data = rich_combined,
    aes(x = Gen_num, y = Shannon, color = Treatment, shape = Treatment),
    size = 2.0, alpha = 0.8
  ) +
  scale_x_continuous(breaks = 1:4, labels = gen_levels, expand = c(0.02, 0.02)) +
  coord_cartesian(xlim = c(0.5, 4.5)) +
  scale_color_manual(values = trt_cols,   drop = FALSE, name = "Treatment") +
  scale_shape_manual(values = trt_shapes, drop = FALSE, name = "Treatment") +
  scale_fill_manual(
    values = c(fill_vals, trt_cols),
    breaks = c("Pre-disturbance", "Post-disturbance"),
    name   = NULL
  ) +
  labs(x = "Generation", y = "Shannon diversity") +
  theme_classic(base_size = 10) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 1),
    fill  = guide_legend(order = 2)
  ) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x     = element_text(margin = margin(t = 8)),
    legend.position  = "right",
    plot.margin      = margin(8, 8, 8, 8)
  )

ggsave("figures/Fig6_shannon_diversity_trajectories.png",
       plot = fig6, width = 6.85, height = 4.8, dpi = 300, bg = "white")


# ─────────────────────────────────────────────────────────────────────────────
# Supplementary Table 1 — SIMPER results
# Fix applied: MD5 hash row names in simper_by_generation.csv have a trailing
# comparison-index digit appended by R for all comparisons after the first.
# substr(taxon_raw, 1, 32) strips the suffix before joining to tax_map.
# ─────────────────────────────────────────────────────────────────────────────

simper_raw <- read.csv("results/simper_by_generation.csv", row.names = 1)

simper_joined <- simper_raw %>%
  rownames_to_column("taxon_raw") %>%
  mutate(taxon = substr(taxon_raw, 1, 32)) %>%
  left_join(tax_map, by = "taxon") %>%
  mutate(
    OrderName = ifelse(is.na(OrderName) | OrderName == "", "Unclassified", OrderName),
    average   = round(as.numeric(average), 4),
    sd        = round(as.numeric(sd),      4),
    ratio     = round(as.numeric(ratio),   3),
    cum_sum   = round(as.numeric(.data$cumsum), 3),  # .data$ avoids base R cumsum() conflict
    p_num     = as.numeric(as.character(p)),
    p         = ifelse(p_num < 0.001, "<0.001", as.character(round(p_num, 3)))
  )

# Rename all columns using base R
names(simper_joined)[names(simper_joined) == "comparison"]  <- "Comparison"
names(simper_joined)[names(simper_joined) == "OrderName"]   <- "Order"
names(simper_joined)[names(simper_joined) == "average"]     <- "Mean contribution"
names(simper_joined)[names(simper_joined) == "sd"]          <- "SD"
names(simper_joined)[names(simper_joined) == "ratio"]       <- "Mean/SD"
names(simper_joined)[names(simper_joined) == "cum_sum"]     <- "Cumulative sum"
names(simper_joined)[names(simper_joined) == "p"]           <- "p-value"

simper_named <- simper_joined[, c("Comparison", "Order", "Mean_contribution",
                                  "SD", "Mean_SD", "Cumulative_sum", "p_value")]

simper_named <- simper_named[order(simper_named$Comparison,
                                   -simper_named$Mean_contribution), ]

# Now rename to display-friendly names for the Word table
colnames(simper_named) <- c("Comparison", "Order", "Mean contribution",
                            "SD", "Mean/SD", "Cumulative sum", "p-value")

head(simper_named)
doc_s1 <- read_docx() %>%
  body_add_par(
    "Supplementary Table 1. SIMPER analysis of bacterial order contributions to Bray-Curtis dissimilarities between generation pairs. Mean contribution is the average contribution of each order to the overall dissimilarity between groups. Cumulative sum indicates the running total of contributions in descending order. P-values are from permutation tests (99 permutations).",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal")

for (comp in unique(simper_named$Comparison)) {
  sub <- simper_named %>%
    filter(Comparison == comp, `Cumulative sum` <= 0.80) %>%
    select(-Comparison)
  
  ft <- flextable(sub) %>%
    bold(part = "header") %>%
    fontsize(size = 9, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    autofit() %>%
    theme_vanilla()
  
  doc_s1 <- doc_s1 %>%
    body_add_par(paste0("Comparison: ", comp), style = "Normal") %>%
    body_add_flextable(ft) %>%
    body_add_par("", style = "Normal")
}

print(doc_s1, target = "results/Supplementary_Table1_SIMPER.docx")
cat("Supplementary Table 1 written\n")


# ─────────────────────────────────────────────────────────────────────────────
# Supplementary Table 2 — PERMANOVA and PERMDISP summary
# ─────────────────────────────────────────────────────────────────────────────

perm_full <- data.frame(
  Analysis = "Full dataset (G0-G3)",
  Term     = c("Generation", "Treatment", "Generation x Treatment", "Residual", "Total"),
  Df       = c(3, 3, 9, 44, 59),
  R2       = c(0.238, 0.056, 0.141, 0.565, 1.000),
  F        = c(6.17, 1.46, 1.22, NA, NA),
  p        = c("0.001", "0.030", "0.046", NA, NA)
)

perm_g1g3 <- data.frame(
  Analysis = "Post-disturbance only (G1-G3)",
  Term     = c("Generation", "Treatment", "Generation x Treatment", "Residual", "Total"),
  Df       = c(2, 3, 6, 33, 44),
  R2       = c(0.064, 0.094, 0.147, 0.695, 1.000),
  F        = c(1.52, 1.49, 1.17, NA, NA),
  p        = c("0.001", "0.001", "0.009", NA, NA)
)

perm_combined <- bind_rows(perm_full, perm_g1g3) %>%
  mutate(
    R2 = ifelse(is.na(R2), "", sprintf("%.3f", R2)),
    F  = ifelse(is.na(F),  "", sprintf("%.2f", F)),
    p  = ifelse(is.na(p),  "", p)
  )

ft_perm <- flextable(perm_combined) %>%
  set_header_labels(Analysis = "Analysis", Term = "Term", Df = "df",
                    R2 = "R\u00b2", F = "F", p = "p-value") %>%
  merge_v(j = "Analysis") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  autofit() %>%
  theme_vanilla()

permdisp <- data.frame(
  Analysis   = c("Full dataset (G0-G3)", "Full dataset (G0-G3)",
                 "Post-disturbance only (G1-G3)", "Post-disturbance only (G1-G3)"),
  Grouping   = c("Generation", "Treatment", "Generation", "Treatment"),
  F          = c(7.49, 1.23, 1.74, 6.32),
  p          = c("<0.001", "0.316", "0.186", "0.002"),
  Conclusion = c("Significant — G0 more dispersed than G1-G3",
                 "Not significant",
                 "Not significant",
                 "Significant — treatment-level dispersion in recovery phase")
)

ft_disp <- flextable(permdisp) %>%
  set_header_labels(Analysis = "Analysis", Grouping = "Grouping factor",
                    F = "F", p = "p-value", Conclusion = "Interpretation") %>%
  merge_v(j = "Analysis") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  autofit() %>%
  theme_vanilla()

doc_s2 <- read_docx() %>%
  body_add_par(
    "Supplementary Table 2. PERMANOVA and permutational multivariate dispersion (PERMDISP) results for Bray-Curtis dissimilarities among adult gut microbiome samples. PERMANOVA was run with 999 permutations using adonis2 (by terms). PERMDISP was evaluated using betadisper and permutest. Results are shown for the full dataset (G0-G3) and the post-disturbance restricted analysis (G1-G3).",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("PERMANOVA results", style = "Normal") %>%
  body_add_flextable(ft_perm) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("PERMDISP results", style = "Normal") %>%
  body_add_flextable(ft_disp) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par(
    "Pairwise PERMDISP comparisons (generation grouping, full dataset): G0 vs G1 p=0.005; G0 vs G2 p=0.008; G0 vs G3 p=0.003; G2 vs G3 p=0.020. All other pairwise comparisons non-significant. Permuted p-values from permutest (999 permutations).",
    style = "Normal"
  )

print(doc_s2, target = "results/Supplementary_Table2_PERMANOVA_PERMDISP.docx")
cat("Supplementary Table 2 written\n")


