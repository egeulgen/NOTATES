#!/usr/bin/env Rscript
# ============================================================================
# Chromosome-arm Somatic Copy Number (SCN) and SCN Alteration (SCNA) tables
# Input : ExomeCNV segment calls (*_cnv.txt) + hg38 chromosome-arm BED
# Output: per-arm tidy table + arm x {p,q} grids (copy number & alteration)
# ============================================================================

# ---- Parameters ------------------------------------------------------------
arguments <- commandArgs(trailingOnly = T)

cnv_file  <- arguments[1]
arm_file  <- arguments[2]
out_dir   <- "chromosome_arm_scna"

baseline_cn   <- 2      # expected diploid copy number
frac_thresh   <- 0.30   # fraction of the CALLED arm that must be altered to call Gain/Loss
min_covered   <- 0.02   # minimum fraction of the arm that must be called, else "Not called"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data -------------------------------------------------------------
seg <- read.delim(cnv_file, header = TRUE,
                  na.strings = c("NA", "", "NaN"))
# ExomeCNV columns: chr, probe_start, probe_end, ..., copy.number, logR, ratio, ...
seg <- data.frame(
  chr   = seg$chr,
  start = as.numeric(seg$probe_start),
  end   = as.numeric(seg$probe_end),
  cn    = suppressWarnings(as.numeric(seg$copy.number)),
  logR  = suppressWarnings(as.numeric(seg$logR))
)

arms <- read.delim(arm_file, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("chr", "start", "end", "arm"))
arms$start <- as.numeric(arms$start)
arms$end   <- as.numeric(arms$end)

# ---- Overlap segments with each arm (split at centromere) ------------------
# For every arm, clamp overlapping segments to the arm boundaries and weight by
# genomic overlap length (bp). This correctly handles segments that span the
# centromere by splitting their bp between the p and q arm.
compute_arm <- function(a) {
  s <- seg[seg$chr == a$chr, , drop = FALSE]
  # overlap of each segment with this arm
  ov_start <- pmax(s$start, a$start)
  ov_end   <- pmin(s$end,   a$end)
  ov_len   <- ov_end - ov_start
  keep     <- ov_len > 0
  s        <- s[keep, , drop = FALSE]
  ov_len   <- ov_len[keep]
  
  arm_len  <- a$end - a$start
  covered  <- sum(ov_len)
  
  # segments with a called (non-NA) copy number
  called   <- !is.na(s$cn)
  called_bp <- sum(ov_len[called])
  
  # length-weighted mean copy number & logR over called segments
  if (called_bp > 0) {
    w        <- ov_len[called]
    wcn      <- sum(s$cn[called]   * w) / sum(w)
    # dominant integer copy-number state (largest covered length)
    cn_by_state <- tapply(w, s$cn[called], sum)
    dom_cn   <- as.numeric(names(cn_by_state)[which.max(cn_by_state)])
    frac_gain <- sum(w[s$cn[called] >  baseline_cn]) / called_bp
    frac_loss <- sum(w[s$cn[called] <  baseline_cn]) / called_bp
    frac_neut <- sum(w[s$cn[called] == baseline_cn]) / called_bp
  } else {
    wcn <- NA_real_; dom_cn <- NA_real_
    frac_gain <- NA_real_; frac_loss <- NA_real_; frac_neut <- NA_real_
  }
  
  logR_ok  <- !is.na(s$logR)
  wlogR    <- if (sum(ov_len[logR_ok]) > 0)
    sum(s$logR[logR_ok] * ov_len[logR_ok]) / sum(ov_len[logR_ok]) else NA_real_
  
  data.frame(
    arm        = a$arm,
    chrom      = sub("chr", "", a$chr),
    p_or_q     = sub(".*chr[0-9XY]+", "", a$arm),  # trailing p / q
    arm_length = arm_len,
    frac_covered = covered / arm_len,
    n_segments = nrow(s),
    weighted_CN = wcn,
    dominant_CN = dom_cn,
    weighted_logR = wlogR,
    frac_gain  = frac_gain,
    frac_loss  = frac_loss,
    frac_neutral = frac_neut
  )
}

res <- do.call(rbind, lapply(seq_len(nrow(arms)), function(i) compute_arm(arms[i, ])))

# ---- Classify each arm: Gain / Loss / Neutral / Not called -----------------
classify <- function(fg, fl, fc) {
  if (is.na(fg) || fc < min_covered)                return("Not called")
  if (fg >= frac_thresh && fg > fl)                 return("Gain")
  if (fl >= frac_thresh && fl > fg)                 return("Loss")
  "Neutral"
}
res$SCNA <- mapply(classify, res$frac_gain, res$frac_loss, res$frac_covered)

# order chromosomes naturally 1..22, X, Y and p before q
chr_levels <- c(as.character(1:22), "X", "Y")
res$chrom  <- factor(res$chrom, levels = chr_levels)
res$p_or_q <- factor(res$p_or_q, levels = c("p", "q"))
res <- res[order(res$chrom, res$p_or_q), ]

# round for readability
res$frac_covered  <- round(res$frac_covered, 3)
res$weighted_CN   <- round(res$weighted_CN, 3)
res$weighted_logR <- round(res$weighted_logR, 4)
res$frac_gain     <- round(res$frac_gain, 3)
res$frac_loss     <- round(res$frac_loss, 3)
res$frac_neutral  <- round(res$frac_neutral, 3)

# ---- Build p/q grids -------------------------------------------------------
to_grid <- function(df, value_col) {
  m <- reshape(df[, c("chrom", "p_or_q", value_col)],
               idvar = "chrom", timevar = "p_or_q", direction = "wide")
  names(m) <- sub(paste0(value_col, "\\."), "", names(m))
  m <- m[order(factor(m$chrom, levels = chr_levels)), ]
  names(m)[names(m) == "chrom"] <- "chromosome"
  if (!"p" %in% names(m)) m$p <- NA
  if (!"q" %in% names(m)) m$q <- NA
  m[, c("chromosome", "p", "q")]
}

grid_cn    <- to_grid(res, "weighted_CN")   # length-weighted arm copy number
grid_domcn <- to_grid(res, "dominant_CN")   # dominant integer copy-number state
grid_scna  <- to_grid(res, "SCNA")          # Gain / Loss / Neutral

# ---- Write outputs ---------------------------------------------------------
write.csv(res,      file.path(out_dir, "arm_scna_table.csv"),        row.names = FALSE)
write.csv(grid_cn,  file.path(out_dir, "arm_copynumber_grid.csv"),   row.names = FALSE)
write.csv(grid_domcn, file.path(out_dir, "arm_dominant_cn_grid.csv"), row.names = FALSE)
write.csv(grid_scna, file.path(out_dir, "arm_scna_grid.csv"),        row.names = FALSE)


# ============================================================================
# Visualise chromosome-arm copy-number & SCNA grids as heatmaps (ggplot2)
# Reads the per-arm table produced by arm_scna.R
# ============================================================================
suppressMessages(library(ggplot2))

# ordering: chromosomes 1..22, X, Y along x ; arm p (top) / q (bottom) along y
chr_levels <- c(as.character(1:22), "X", "Y")
res$chrom  <- factor(as.character(res$chrom), levels = chr_levels)
res$p_or_q <- factor(res$p_or_q, levels = c("q", "p"))   # q first => p ends up on top

# ---------------------------------------------------------------------------
# Plot 1: length-weighted arm copy number (diverging scale centred on CN=2)
# ---------------------------------------------------------------------------
# clip fill to a tight window around baseline so subtle arm-level shifts show,
# while the printed number gives the exact value.
lim <- 0.30
res$fill_cn <- pmin(pmax(res$weighted_CN, baseline_cn - lim), baseline_cn + lim)

p_cn <- ggplot(res, aes(chrom, p_or_q, fill = fill_cn)) +
  geom_tile(color = "white", height = 0.6) +
  geom_text(aes(label = ifelse(is.na(weighted_CN), "NC",
                               sprintf("%.2f", weighted_CN))),
            size = 2.7, color = "grey15") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = baseline_cn,
                       limits = c(baseline_cn - lim, baseline_cn + lim),
                       na.value = "grey80",
                       name = "Weighted\ncopy number",
                       breaks = c(1.7, 2.0, 2.3),
                       labels = c("\u22641.7", "2.0", "\u22652.3")) +
  labs(title = "Chromosome-arm somatic copy number",
       subtitle = "Length-weighted mean CN per arm (NC = not covered by exome). Baseline diploid = 2.",
       x = "Chromosome", y = "Arm") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "grey20"),
        legend.position = "right")

ggsave(file.path(out_dir, "arm_copynumber_heatmap.png"), p_cn,
       width = 11, height = 2.8, dpi = 200, bg = "white")

# ---------------------------------------------------------------------------
# Plot 2: categorical SCNA call (Gain / Loss / Neutral / Not called)
# ---------------------------------------------------------------------------
res$SCNA <- factor(res$SCNA, levels = c("Loss", "Neutral", "Gain", "Not called"))
scna_cols <- c("Loss" = "#2166AC", "Neutral" = "#F0F0F0",
               "Gain" = "#B2182B", "Not called" = "grey65")

p_scna <- ggplot(res, aes(chrom, p_or_q, fill = SCNA)) +
  geom_tile(color = "white", height = 0.6) +
  geom_text(aes(label = substr(as.character(SCNA), 1, 1)),
            size = 2.7, color = "grey20") +
  scale_fill_manual(values = scna_cols, name = "SCNA call", drop = FALSE) +
  labs(title = "Chromosome-arm somatic copy number alterations",
       subtitle = paste0("Arm called Gain/Loss when \u226530% of covered length is altered vs CN=2 ",
                         "(G/L/N/N = Gain/Loss/Neutral/Not called)."),
       x = "Chromosome", y = "Arm") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "grey20"),
        legend.position = "right")

ggsave(file.path(out_dir, "arm_scna_heatmap.png"), p_scna,
       width = 11, height = 2.8, dpi = 200, bg = "white")

# ---------------------------------------------------------------------------
# Plot 3: fraction of each arm gained / lost (surfaces subclonal/partial events
# that the near-2 weighted mean hides, e.g. the partial chr7 gain)
# ---------------------------------------------------------------------------
res$signed_frac <- ifelse(is.na(res$frac_gain), NA, res$frac_gain - res$frac_loss)

p_frac <- ggplot(res, aes(chrom, p_or_q, fill = signed_frac)) +
  geom_tile(color = "white", height = 0.6) +
  geom_text(aes(label = ifelse(is.na(signed_frac), "NC",
                               ifelse(abs(signed_frac) < 0.005, "",
                                      sprintf("%+.0f%%", 100 * signed_frac)))),
            size = 2.6, color = "grey15") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), na.value = "grey80",
                       name = "Arm fraction\n(gain \u2212 loss)",
                       breaks = c(-1, 0, 1),
                       labels = c("100% loss", "0", "100% gain")) +
  labs(title = "Fraction of each arm gained vs lost",
       subtitle = "Signed fraction of covered length altered; reveals partial / subclonal events.",
       x = "Chromosome", y = "Arm") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "grey20"),
        legend.position = "right")

ggsave(file.path(out_dir, "arm_gain_loss_fraction_heatmap.png"), p_frac,
       width = 11, height = 2.8, dpi = 200, bg = "white")
