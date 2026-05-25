# Assumes 00_setup.R, 01_processing.R, and 02_analysis.R have already been sourced

dir.create("figures", showWarnings = FALSE)

#===Figure 2. Distinct Microbial Profile & Alpha Diversity (Observed and Shannon) by Treatment===
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
    x = "Treatment assignment",
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
  height = 4.5,
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
ggsave("figures/Fig3a_flour_order_composition.png",
  plot = fig3a_flour,
  width = 6.85,
  height = 3.4,
  dpi = 300,
  bg = "white"
)

# Fig. 3b: adult gut panel
ggsave("figures/Fig3b_adultgut_order_composition.png",
  plot = fig3b_gut,
  width = 6.85,
  height = 3.6,
  dpi = 300,
  bg = "white"
)

#===FIG 4: Alpha diversity dynamics across generations========================
# Observed richness (top) + Shannon diversity (bottom)
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
  ggplot(aes(x = Generation, y = Value, color = Generation)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_jitter(
    aes(shape = Treatment),
    width = 0.12, height = 0,
    alpha = 0.7, size = 1.6
  ) +
  facet_wrap(~Treatment, nrow = 1) +
  scale_color_manual(values = gen_cols, drop = FALSE, name = "Generation") +
  scale_shape_manual(values = trt_shapes, drop = FALSE, name = "Treatment") +
  geom_text(
    data = kw_results_fig4 %>% filter(Metric == "Observed richness"),
    aes(x = 1, y = Inf, label = p_label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.1,
    size = 2.9,
    color = "black"
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
    data = kw_results_fig4 %>% filter(Metric == "Shannon diversity"),
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
ggsave("figures/Fig4_alpha_diversity_dynamics.png",
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

ggsave("figures/Fig5_PCoA_AdultGut.png",
  plot = pcoa_plot1,
  width = 6.85,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

#=== Reviewer-requested Figure: Beta diversity during recovery phase only, G1-G3 ===
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

ggsave("figures/Reviewer_Fig_G1G3_PCoA_AdultGut.png",
  plot = pcoa_plot_G1G3,
  width = 6.85,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

#=== Figure 6. Diet-shift differential abundance + order-level heatmap ===

labels_row_pretty <- labels_row
labels_row_pretty <- gsub("↑O", "[OAT↑]", labels_row_pretty)
labels_row_pretty <- gsub("↓O", "[OAT↓]", labels_row_pretty)
labels_row_pretty <- gsub("↑W", "[WHEAT↑]", labels_row_pretty)
labels_row_pretty <- gsub("↓W", "[WHEAT↓]", labels_row_pretty)

heat_colors <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100)

png(
  filename = "figures/Fig6_diet_shift_heatmap.png",
  width = 2055,
  height = 2500,
  res = 300
)

pheatmap(
  mat = as.matrix(heatmap_mat_filtered),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 7,
  labels_row = labels_row_pretty,
  labels_col = col_labels,
  color = heat_colors,
  border_color = "grey80",
  angle_col = 90
)

dev.off()

#=== Figure 7. Shannon diversity trajectories across generations ============
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

ggsave("figures/Fig7_shannon_diversity_trajectories.png",
  plot = fig7,
  width = 6.85,
  height = 4.8,
  dpi = 300,
  bg = "white"
)


