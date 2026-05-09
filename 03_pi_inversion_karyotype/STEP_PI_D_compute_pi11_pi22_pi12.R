#!/usr/bin/env Rscript
# =============================================================================
# compute_pi11_pi22_pi12.R
# =============================================================================
# Compute sitewise and windowed nucleotide diversity / divergence between
# the two homokaryotype bands, given ANGSD .mafs.gz outputs that share the
# same major/minor allele coding (via -doMajorMinor 3 -sites <4-col file>).
#
# Per biallelic site, with the *same* allele defined as "minor" in both groups:
#
#   pi11 = 2 * p1 * (1 - p1)             # Homo_1 within-group diversity
#   pi22 = 2 * p2 * (1 - p2)             # Homo_2 within-group diversity
#   pi12 = p1 * (1 - p2) + p2 * (1 - p1) # pairwise divergence Homo_1 vs Homo_2
#
# pi12 is NOT the heterozygote group; it is computed entirely from p1, p2.
#
# Site basis warning
#   If the shared sites file contains only variant sites (the default when
#   building it from -doMaf 1 with -SNP_pval), all three pi values are
#   "per-variant-site" and inflated relative to true per-callable-site pi.
#   The --site-basis flag is recorded in outputs so downstream code can
#   rescale. To get unbiased per-callable-site pi, the caller must instead
#   provide a sites file that includes monomorphic callable sites (which
#   contribute 0 to all three pi).
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--homo1",        type = "character"),
  make_option("--homo2",        type = "character"),
  make_option("--inv",          type = "character"),
  make_option("--chr",          type = "character"),
  make_option("--start",        type = "integer"),
  make_option("--end",          type = "integer"),
  make_option("--win",          type = "integer", default = 50000L),
  make_option("--step",         type = "integer", default = 10000L),
  make_option("--site-basis",   type = "character", default = "variant_only",
              help = "'variant_only' or 'all_callable' [default %default]"),
  make_option("--out-sitewise", type = "character"),
  make_option("--out-window",   type = "character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(
  !is.null(opt$homo1), !is.null(opt$homo2), !is.null(opt$inv),
  !is.null(opt$chr),   !is.null(opt$start), !is.null(opt$end),
  !is.null(opt$`out-sitewise`), !is.null(opt$`out-window`)
)
if (!opt$`site-basis` %in% c("variant_only", "all_callable")) {
  stop("--site-basis must be 'variant_only' or 'all_callable'")
}

# ── Read mafs.gz ────────────────────────────────────────────────────────────
read_mafs <- function(path) {
  if (!file.exists(path)) stop(sprintf("mafs file not found: %s", path))
  m <- fread(path, sep = "\t", header = TRUE)
  # ANGSD writes one of: knownEM, unknownEM, phat (depending on flags).
  freq_col <- intersect(c("knownEM", "unknownEM", "phat"), colnames(m))[1]
  if (is.na(freq_col)) {
    stop(sprintf("No freq column (knownEM/unknownEM/phat) in %s; columns: %s",
                 path, paste(colnames(m), collapse = ",")))
  }
  m[, .(chromo, position, major, minor, freq = get(freq_col))]
}

m1 <- read_mafs(opt$homo1)
m2 <- read_mafs(opt$homo2)
setnames(m1, "freq", "p1")
setnames(m2, "freq", "p2")

# ── Merge on (chromo, position, major, minor) ──────────────────────────────
# This guarantees we are comparing frequencies of the SAME allele (minor)
# in both groups. If the shared sites pipeline worked, every row in m1 and
# m2 should match; we use inner join and report any drop.
key_cols <- c("chromo", "position", "major", "minor")
setkeyv(m1, key_cols); setkeyv(m2, key_cols)
mm <- merge(m1, m2, by = key_cols, all = FALSE)

n1 <- nrow(m1); n2 <- nrow(m2); nm <- nrow(mm)
if (nm < min(n1, n2)) {
  message(sprintf(
    "[compute_pi] dropped sites during merge: Homo_1=%d Homo_2=%d intersect=%d",
    n1, n2, nm))
}
if (nm == 0L) stop("No shared sites after merge; check shared sites file.")

# Restrict to interval just in case ANGSD wrote a few sites outside -r
mm <- mm[chromo == opt$chr & position >= opt$start & position <= opt$end]

# ── Sitewise pi ─────────────────────────────────────────────────────────────
# p1, p2 are minor-allele frequencies of the SAME minor allele in both groups.
mm[, pi11 := 2 * p1 * (1 - p1)]
mm[, pi22 := 2 * p2 * (1 - p2)]
mm[, pi12 := p1 * (1 - p2) + p2 * (1 - p1)]

mm[, inversion_id := opt$inv]
mm[, site_basis   := opt$`site-basis`]

setcolorder(mm, c("inversion_id", "chromo", "position", "major", "minor",
                  "p1", "p2", "pi11", "pi22", "pi12", "site_basis"))

# Sitewise output (gzipped TSV)
fwrite(mm, opt$`out-sitewise`, sep = "\t",
       compress = "gzip", na = "NA", quote = FALSE)
message(sprintf("[compute_pi] sitewise: %d sites -> %s",
                nrow(mm), opt$`out-sitewise`))

# ── Windowed pi ────────────────────────────────────────────────────────────
win  <- as.integer(opt$win)
step <- as.integer(opt$step)
stopifnot(win > 0L, step > 0L)

# Build sliding windows covering [start, end]; window = [w_start, w_start+win)
w_starts <- seq.int(from = as.integer(opt$start),
                    to   = as.integer(opt$end) - 1L,
                    by   = step)
windows <- data.table(
  inversion_id = opt$inv,
  chromo       = opt$chr,
  win_start    = w_starts,
  win_end      = w_starts + win - 1L,
  win_center   = w_starts + (win %/% 2L)
)

# Aggregate sites into windows
sites <- mm[, .(position, pi11, pi22, pi12)]
agg <- vector("list", nrow(windows))
for (i in seq_len(nrow(windows))) {
  ws <- windows$win_start[i]; we <- windows$win_end[i]
  s  <- sites[position >= ws & position <= we]
  if (nrow(s) == 0L) {
    agg[[i]] <- data.table(
      n_sites   = 0L,
      mean_pi11 = NA_real_, mean_pi22 = NA_real_, mean_pi12 = NA_real_,
      sum_pi11  = NA_real_, sum_pi22  = NA_real_, sum_pi12  = NA_real_
    )
  } else {
    agg[[i]] <- data.table(
      n_sites   = nrow(s),
      mean_pi11 = mean(s$pi11, na.rm = TRUE),
      mean_pi22 = mean(s$pi22, na.rm = TRUE),
      mean_pi12 = mean(s$pi12, na.rm = TRUE),
      sum_pi11  = sum(s$pi11,  na.rm = TRUE),
      sum_pi22  = sum(s$pi22,  na.rm = TRUE),
      sum_pi12  = sum(s$pi12,  na.rm = TRUE)
    )
  }
}
win_dt <- cbind(windows, rbindlist(agg))
win_dt[, site_basis := opt$`site-basis`]

fwrite(win_dt, opt$`out-window`, sep = "\t", na = "NA", quote = FALSE)
message(sprintf("[compute_pi] windowed: %d windows (win=%d step=%d) -> %s",
                nrow(win_dt), win, step, opt$`out-window`))
