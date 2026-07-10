#!/usr/bin/env Rscript
# ============================================================================
# Chromosome-arm Somatic Copy Number (SCN) and SCN Alteration (SCNA) tables
# Input : ExomeCNV segment calls (*_cnv.txt) + hg38 chromosome-arm BED
# Output: per-arm tidy table + arm x {p,q} grids + heatmaps
# Drives arm calls off the continuous, recentred log2(tumor/normal) ratio
# (robust to low tumour purity, which collapses ExomeCNV's integer copy.number
# toward 2 and caps focal amplifications such as CDK4).
# ============================================================================
suppressMessages(library(ggplot2))

# ---- arguments -------------------------------------------------------------
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2)
  stop("Usage: Rscript arm_scna_pipeline.R <cnv_file> <arm_bed> [out_dir]")
cnv_file <- arguments[1]
arm_file <- arguments[2]
out_dir  <- ifelse(length(arguments) >= 3, arguments[3], "chromosome_arm_scna")

# ---- tunable parameters ----------------------------------------------------
seg_lr_cut   <- 0.25   # per-segment log2 cutoff (matches the report's SCNA def)
arm_lr_cut   <- 0.15   # |arm-mean recentred log2| beyond this => Gain/Loss call
mad_k        <- 2      # arm also flagged if beyond median +/- k*MAD of arm means
focal_lr_cut <- 0.80   # |seg log2| >= this => candidate focal event
focal_min_kb <- 100    # ...and at least this wide (excludes single-probe noise spikes)
focal_max_mb <- 20     # ...and at most this wide (else it's a broad, not focal, event)
min_covered  <- 0.05   # min fraction of arm covered, else "Not called"
lim          <- 0.30   # heatmap colour-scale window around 0 (log2 units)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load ------------------------------------------------------------------
raw <- read.delim(cnv_file, header = TRUE, na.strings = c("NA", "", "NaN"))
seg <- data.frame(
  chr   = raw$chr,
  start = as.numeric(raw$probe_start),
  end   = as.numeric(raw$probe_end),
  cn    = suppressWarnings(as.numeric(raw$copy.number)),
  logR  = suppressWarnings(as.numeric(raw$logR)),
  ratio = suppressWarnings(as.numeric(raw$ratio))
)
seg$len <- seg$end - seg$start

arms <- read.delim(arm_file, header = FALSE,
                   col.names = c("chr", "start", "end", "arm"))
arms$start <- as.numeric(arms$start); arms$end <- as.numeric(arms$end)

# ---- robust baseline recentring -------------------------------------------
# Subtract the bp-weighted median segment log2 (autosomes) so a balanced
# genome sits at 0, guarding against any global normalisation offset.
au <- seg[grepl("^chr[0-9]+$", seg$chr) & !is.na(seg$logR), ]
o  <- order(au$logR); cs <- cumsum(as.numeric(au$len[o]))
baseline <- au$logR[o][which(cs >= sum(as.numeric(au$len)) / 2)[1]]
seg$logR_c <- seg$logR - baseline
cat(sprintf("Global baseline (bp-weighted median log2) = %.4f  [subtracted]\n", baseline))

# ---- per-arm aggregation on the continuous signal --------------------------
compute_arm <- function(a) {
  s <- seg[seg$chr == a$chr, , drop = FALSE]
  ovl <- pmin(s$end, a$end) - pmax(s$start, a$start)
  keep <- ovl > 0; s <- s[keep, ]; ovl <- ovl[keep]
  arm_len <- a$end - a$start
  
  ok  <- !is.na(s$logR_c)
  w   <- ovl[ok]; lr <- s$logR_c[ok]
  covered   <- sum(ovl)
  called_bp <- sum(w)
  
  if (called_bp > 0) {
    wmean_lr  <- sum(lr * w) / called_bp
    rel_cn    <- 2 * 2^wmean_lr                 # relative CN implied by ratio (NOT purity-corrected)
    frac_gain <- sum(w[lr >=  seg_lr_cut]) / called_bp
    frac_loss <- sum(w[lr <= -seg_lr_cut]) / called_bp
  } else {
    wmean_lr <- NA_real_; rel_cn <- NA_real_; frac_gain <- NA_real_; frac_loss <- NA_real_
  }
  data.frame(arm = a$arm,
             chrom = sub("chr", "", a$chr),
             p_or_q = sub(".*chr[0-9XY]+", "", a$arm),
             frac_covered = covered / arm_len,
             n_seg = nrow(s),
             wmean_log2 = wmean_lr,
             rel_CN = rel_cn,
             frac_gain = frac_gain,
             frac_loss = frac_loss)
}
res <- do.call(rbind, lapply(seq_len(nrow(arms)), function(i) compute_arm(arms[i, ])))

# ---- arm-level Gain/Loss call ---------------------------------------------
# adaptive floor: combine a fixed cutoff with a MAD-based cutoff of arm means
arm_mad <- mad(res$wmean_log2, na.rm = TRUE)
cut_eff <- max(arm_lr_cut, mad_k * arm_mad)
cat(sprintf("Arm-mean MAD = %.4f ; effective |log2| call threshold = %.3f\n",
            arm_mad, cut_eff))

res$SCNA <- with(res, ifelse(is.na(wmean_log2) | frac_covered < min_covered, "Not called",
                             ifelse(wmean_log2 >=  cut_eff, "Gain",
                                    ifelse(wmean_log2 <= -cut_eff, "Loss", "Neutral"))))

# ---- focal high-amplitude events (rescued from arm averaging) --------------
foc <- seg[is.finite(seg$logR_c) & is.finite(seg$ratio) & seg$ratio > 0 &
             abs(seg$logR_c) >= focal_lr_cut &
             seg$len >= focal_min_kb * 1e3 &
             seg$len <= focal_max_mb * 1e6, ]
foc <- foc[order(-abs(foc$logR_c)), ]
focal <- data.frame(
  region = sprintf("%s:%d-%d", foc$chr, foc$start, foc$end),
  width_kb = round(foc$len / 1e3, 1),
  log2_recentred = round(foc$logR_c, 3),
  ratio = round(foc$ratio, 3),
  rel_CN = round(2 * 2^foc$logR_c, 2),
  exomeCNV_int_CN = foc$cn,
  type = ifelse(foc$logR_c > 0, "amplification", "deletion"))

# ---- ordering & rounding ---------------------------------------------------
chr_lv <- c(as.character(1:22), "X", "Y")
res$chrom  <- factor(res$chrom, levels = chr_lv)
res$p_or_q <- factor(res$p_or_q, levels = c("p", "q"))
res <- res[order(res$chrom, res$p_or_q), ]
for (col in c("frac_covered","wmean_log2","rel_CN","frac_gain","frac_loss"))
  res[[col]] <- round(res[[col]], 3)

# ---- grids -----------------------------------------------------------------
to_grid <- function(df, v) {
  m <- reshape(df[, c("chrom","p_or_q",v)], idvar="chrom", timevar="p_or_q", direction="wide")
  names(m) <- sub(paste0(v,"\\."), "", names(m))
  m <- m[order(factor(m$chrom, levels=chr_lv)), ]
  names(m)[1] <- "chromosome"; for (x in c("p","q")) if (!x %in% names(m)) m[[x]] <- NA
  m[, c("chromosome","p","q")]
}
grid_lr   <- to_grid(res, "wmean_log2")
grid_scna <- to_grid(res, "SCNA")

# ---- write -----------------------------------------------------------------
write.csv(res,       file.path(out_dir, "arm_scna_table_robust.csv"),  row.names = FALSE)
write.csv(grid_lr,   file.path(out_dir, "arm_log2ratio_grid.csv"),     row.names = FALSE)
write.csv(grid_scna, file.path(out_dir, "arm_scna_grid_robust.csv"),   row.names = FALSE)
write.csv(focal,     file.path(out_dir, "focal_events.csv"),           row.names = FALSE)

# ============================================================================
# Visualise the arm grids as heatmaps (ggplot2)
# ============================================================================
# plotting order: chromosomes 1..22, X, Y along x ; p on top, q below
res$chrom  <- factor(as.character(res$chrom), levels = chr_lv)
res$p_or_q <- factor(as.character(res$p_or_q), levels = c("q", "p"))

# flag arms that contain a focal event, so CDK4 / chr1q etc. are visible
res$focal <- FALSE
if (nrow(foc) > 0) {
  arm_of <- function(chr, pos) {
    a <- arms[arms$chr == chr & arms$start <= pos & arms$end > pos, "arm"]
    if (length(a)) a[1] else NA_character_
  }
  foc_arm <- mapply(arm_of, foc$chr, foc$start)
  res$focal <- res$arm %in% na.omit(foc_arm)
}

# ---------------------------------------------------------------------------
# Plot 1: arm relative copy number = length-weighted mean log2(tumor/normal)
# ---------------------------------------------------------------------------
res$fill_lr <- pmin(pmax(res$wmean_log2, -lim), lim)

p_lr <- ggplot(res, aes(chrom, p_or_q, fill = fill_lr)) +
  geom_tile(color = "white", height = 0.6) +
  geom_tile(data = subset(res, focal), color = "black", height = 0.9, fill = NA) +
  geom_text(aes(label = ifelse(is.na(wmean_log2), "NC", sprintf("%+.2f", wmean_log2))),
            size = 2.6, color = "grey15") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim), na.value = "grey80",
                       name = "mean log2\n(tumor/normal)",
                       breaks = c(-lim, 0, lim),
                       labels = c(sprintf("\u2264-%.1f", lim), "0", sprintf("\u2265+%.1f", lim))) +
  labs(title = "Chromosome-arm relative copy number",
       subtitle = "Length-weighted mean log2(tumor/normal), recentred. Black outline = arm contains a focal event.",
       x = "Chromosome", y = "Arm") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "grey20"), legend.position = "right")

ggsave(file.path(out_dir, "arm_log2ratio_heatmap.png"), p_lr,
       width = 11, height = 2.9, dpi = 200, bg = "white")

# ---------------------------------------------------------------------------
# Plot 2: categorical SCNA call (Gain / Loss / Neutral / Not called)
# ---------------------------------------------------------------------------
res$SCNA <- factor(res$SCNA, levels = c("Loss", "Neutral", "Gain", "Not called"))
scna_cols <- c("Loss" = "#2166AC", "Neutral" = "#F0F0F0",
               "Gain" = "#B2182B", "Not called" = "grey65")

p_scna <- ggplot(res, aes(chrom, p_or_q, fill = SCNA)) +
  geom_tile(color = "white", height = 0.6) +
  geom_tile(data = subset(res, focal), color = "black", height = 0.9, fill = NA) +
  geom_text(aes(label = substr(as.character(SCNA), 1, 1)), size = 2.7, color = "grey20") +
  scale_fill_manual(values = scna_cols, name = "SCNA call", drop = FALSE) +
  labs(title = "Chromosome-arm somatic copy number alterations (purity-robust)",
       subtitle = sprintf("Arm called Gain/Loss when |arm-mean log2(tumor/normal)| \u2265 %.2f. Black outline = arm harbours a focal event.", cut_eff),
       x = "Chromosome", y = "Arm") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "grey20"), legend.position = "right")

ggsave(file.path(out_dir, "arm_scna_heatmap.png"), p_scna,
       width = 11, height = 2.9, dpi = 200, bg = "white")

# ---------------------------------------------------------------------------
# Plot 3: signed fraction of each arm gained vs lost (surfaces partial events)
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
       subtitle = sprintf("Signed fraction of covered length with |segment log2| \u2265 %.2f; reveals partial / subclonal events.", seg_lr_cut),
       x = "Chromosome", y = "Arm") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "grey20"), legend.position = "right")

ggsave(file.path(out_dir, "arm_gain_loss_fraction_heatmap.png"), p_frac,
       width = 11, height = 2.8, dpi = 200, bg = "white")

cat("Done. Wrote tables and 3 heatmaps to: ", out_dir, "\n", sep = "")