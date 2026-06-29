# Assumes 00_setup.R, 01_processing.R, and 02_analysis_integrated.R have already been sourced.
# REVISED: Fig 2 sample-level composition bars, Fig 3 alpha diversity,
#          Fig 4 two-panel PCoA, Fig 5 heatmap, Fig 6 single global baseline
#          Uses concordance-filtered labels_row from 02_analysis_integrated.R

dir.create("figures", showWarnings = FALSE)

# -------------------------
# ARCHIVED: G0 baseline alpha diversity
# Removed from main paper per Reviewer 3.
# Statistics are reported in text.
# -------------------------
# p_fig_G0_baseline <- ggplot(alpha_G0_long, aes(x = Treatment, y = Value)) +
#   geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
#   geom_jitter(width = 0.12, height = 0, alpha = 0.7, size = 1.4) +
#   facet_wrap(~Metric, scales = "free_y", nrow = 1) +
#   geom_text(
#     data = kw_results_baseline_G0,
#     aes(x = 1, y = Inf, label = p_label),
#     inherit.aes = FALSE,
#     hjust = 0, vjust = 1.1, size = 3
#   ) +
#   labs(x = "Downstream treatment assignment", y = NULL) +
#   theme_classic(base_size = 10) +
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_text(face = "plain"),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#     axis.title.x = element_text(margin = margin(t = 8)),
#     panel.spacing = unit(1.2, "lines")
#   ) +
#   coord_cartesian(clip = "off")
#
# ggsave("figures/Fig_G0_baseline_alpha_ARCHIVED.png",
#        plot = p_fig_G0_baseline,
#        width = 6.85, height = 3.8, dpi = 300, bg = "white")

# -------------------------
# Figure 2 — shared flour/gut colour palette setup
# (required by both Fig 2a and 2b)
# -------------------------
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

if ("<0.5%" %in% order_levels) {
  order_levels <- c(setdiff(order_levels, "<0.5%"), "<0.5%")
}

shared_colors     <- colorRampPalette(brewer.pal(12, "Set3"))(length(order_levels))
names(shared_colors) <- order_levels
shared_color_scale   <- scale_fill_manual(values = shared_colors, drop = FALSE)

# -------------------------
# Figure 2a — flour substrates: sample-level bars
# -------------------------
df_flour_sample_plot <- df_flour_raw %>%
  group_by(Sample, Description, Order) %>%
  summarise(Count = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Total = sum(Count), Relative = Count / Total) %>%
  ungroup() %>%
  mutate(
    Order = factor(Order, levels = order_levels),
    Description = factor(Description,
                         levels = flour_types,
                         labels = c("Whole wheat\nfresh", "Whole wheat\nused",
                                    "Oat\nfresh", "Oat\nused"))
  )

fig2a_flour <- ggplot(df_flour_sample_plot,
                      aes(x = Sample, y = Relative, fill = Order)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  facet_wrap(~Description, scales = "free_x", nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  shared_color_scale +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "right"
  )

ggsave("figures/Fig2a_flour_order_composition.png",
       plot = fig2a_flour, width = 6.85, height = 3.4, dpi = 300, bg = "white")

# -------------------------
# Figure 2b — gut microbiome: sample-level bars
# -------------------------
df_gut_sample_plot <- df_gut_raw %>%
  group_by(Sample, Treatment, Generation, Order) %>%
  summarise(Count = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Total = sum(Count), Relative = Count / Total) %>%
  ungroup() %>%
  mutate(
    Order      = factor(Order, levels = order_levels),
    Treatment  = factor(Treatment,  levels = c("ControlA", "ExpA", "ControlB", "ExpB")),
    Generation = factor(Generation, levels = c("G0", "G1", "G2", "G3"))
  )

fig2b_gut <- ggplot(df_gut_sample_plot,
                    aes(x = Sample, y = Relative, fill = Order)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9) +
  facet_grid(rows = vars(Treatment), cols = vars(Generation),
             scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  shared_color_scale +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.4, "lines"),
    legend.position = "right"
  )

ggsave("figures/Fig2b_adultgut_order_composition.png",
       plot = fig2b_gut, width = 8.5, height = 6.0, dpi = 300, bg = "white")

# -------------------------
# Figure 3 — alpha diversity dynamics
# -------------------------
fig3_theme <- theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 8)),
    panel.spacing = unit(1, "lines"),
    legend.position = "right",
    plot.margin = margin(6, 6, 6, 10)
  )

p_observed <- alpha_long %>%
  filter(Metric == "Observed richness") %>%
  ggplot(aes(x = Generation, y = Value)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_jitter(aes(color = Generation, shape = Treatment),
              width = 0.12, height = 0, alpha = 0.7, size = 1.6) +
  facet_wrap(~Treatment, nrow = 1) +
  scale_color_manual(values = gen_cols, drop = FALSE, name = "Generation") +
  scale_shape_manual(values = trt_shapes, drop = FALSE, name = "Treatment") +
  geom_text(
    data = kw_results_fig3 %>% filter(Metric == "Observed richness"),
    aes(x = 1, y = Inf, label = p_label),
    inherit.aes = FALSE, hjust = 0, vjust = 1.1, size = 2.9
  ) +
  labs(x = "Generation", y = "Observed richness") +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  fig3_theme

p_shannon <- alpha_long %>%
  filter(Metric == "Shannon diversity") %>%
  ggplot(aes(x = Generation, y = Value)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_jitter(aes(color = Generation, shape = Treatment),
              width = 0.12, height = 0, alpha = 0.7, size = 1.6) +
  facet_wrap(~Treatment, nrow = 1) +
  scale_color_manual(values = gen_cols, drop = FALSE, name = "Generation") +
  scale_shape_manual(values = trt_shapes, drop = FALSE, name = "Treatment") +
  geom_text(
    data = kw_results_fig3 %>% filter(Metric == "Shannon diversity"),
    aes(x = 1, y = Inf, label = p_label),
    inherit.aes = FALSE, hjust = 0, vjust = 1.1, size = 2.9
  ) +
  labs(x = "Generation", y = "Shannon diversity") +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  fig3_theme

fig3 <- p_observed / p_shannon +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("figures/Fig3_alpha_diversity_dynamics.png",
       plot = fig3, width = 6.85, height = 6.4, dpi = 300, bg = "white")

# -------------------------
# Figure 4 — two-panel PCoA: all generations + G1-G3 only
# -------------------------
pcoa_bray <- ordinate(ps_AdultGut, method = "PCoA", distance = "bray")
pc1_var   <- round(100 * pcoa_bray$values$Relative_eig[1], 1)
pc2_var   <- round(100 * pcoa_bray$values$Relative_eig[2], 1)

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
  theme(legend.position = "right",
        plot.title = element_text(size = 9, face = "plain"),
        plot.margin = margin(8, 8, 8, 8)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

ps_G1G3    <- subset_samples(ps_AdultGut, Generation %in% c("G1", "G2", "G3"))
ps_G1G3    <- prune_taxa(taxa_sums(ps_G1G3) > 0, ps_G1G3)
pcoa_G1G3  <- ordinate(ps_G1G3, method = "PCoA", distance = "bray")
pc1_G1G3   <- round(100 * pcoa_G1G3$values$Relative_eig[1], 1)
pc2_G1G3   <- round(100 * pcoa_G1G3$values$Relative_eig[2], 1)
gen_cols_G1G3 <- gen_cols[c("G1", "G2", "G3")]

fig4b <- plot_ordination(ps_G1G3, pcoa_G1G3, color = "Generation") +
  geom_point(aes(shape = Treatment), size = 2.6, alpha = 0.9) +
  stat_ellipse(aes(group = Generation), type = "t",
               linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  scale_shape_manual(values = trt_shapes,   drop = FALSE, name = "Treatment") +
  scale_color_manual(values = gen_cols_G1G3, drop = FALSE, name = "Generation") +
  labs(x = paste0("PCoA 1 (", pc1_G1G3, "%)"),
       y = paste0("PCoA 2 (", pc2_G1G3, "%)"),
       title = "Post-disturbance only (G1\u2013G3)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "right",
        plot.title = element_text(size = 9, face = "plain"),
        plot.margin = margin(8, 8, 8, 8)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

fig4_combined <- fig4a + fig4b +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("figures/Fig4_PCoA_AdultGut.png",
       plot = fig4_combined, width = 13.0, height = 4.8, dpi = 300, bg = "white")

# -------------------------
# Figure 5 — heatmap
# labels_row and heatmap_mat_filtered are computed in 02_analysis_integrated.R
# using concordance-filtered merged_df — no concordance code needed here
# -------------------------
heat_colors <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100)

png(filename = "figures/Fig5_diet_shift_heatmap.png",
    width = 2055, height = 1650, res = 300)

pheatmap(
  mat          = as.matrix(heatmap_mat_filtered),
  scale        = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 7,
  labels_row   = labels_row,
  labels_col   = col_labels,
  color        = heat_colors,
  border_color = "grey85",
  angle_col    = 45
)

dev.off()

# -------------------------
# Figure 6 — Shannon diversity trajectories
# UPDATED: single global G0 baseline (Reviewer 1)
# -------------------------
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

global_baseline <- rich_combined %>%
  filter(Generation == "G0") %>%
  summarise(y0 = mean(Shannon, na.rm = TRUE)) %>%
  pull(y0)

disturbance <- data.frame(
  Period = c("Pre-disturbance", "Post-disturbance"),
  xmin = c(0.5, 1.5), xmax = c(1.5, 4.5),
  ymin = -Inf, ymax = Inf
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
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "right",
    plot.margin  = margin(8, 8, 8, 8)
  )

ggsave("figures/Fig6_shannon_diversity_trajectories.png",
       plot = fig6, width = 6.85, height = 4.8, dpi = 300, bg = "white")
