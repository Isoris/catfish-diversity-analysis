#!/usr/bin/env Rscript
# =============================================================================
# average_sample_pestPG_by_karyotype.R
# =============================================================================
# Average per-sample ANGSD pestPG tP tracks by karyotype group, restricted
# to one inversion interval. NOTE: this is *not* a formal group pi; it is
# the mean of within-individual diversity estimates across samples in a
# band. Label outputs accordingly: mean_per_sample_tP / mean_within_individual_diversity.
#
# Inputs
#   --karyotype   TSV: sample, inversion_id, group
#   --pestpg-dir  Directory containing per-sample windowed pestPG files,
#                 one file per sample. Files are produced by STEP_A02 as:
#                   <SAMPLE>.win<WIN>.step<STEP>.pestPG
#                 We accept any file matching <SAMPLE>.*.pestPG and pick
#                 the first match.
#   --inversion   inversion_id
#   --chr, --start, --end  interval to restrict to
#   --out         output TSV
#
# Output columns
#   inversion_id, chromo, WinStart, WinStop, WinCenter, group,
#   n_samples, mean_tP, median_tP, sd_tP, se_tP
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--karyotype",  type = "character"),
  make_option("--pestpg-dir", type = "character"),
  make_option("--inversion",  type = "character"),
  make_option("--chr",        type = "character"),
  make_option("--start",      type = "integer"),
  make_option("--end",        type = "integer"),
  make_option("--out",        type = "character"),
  make_option("--win-pattern", type = "character", default = NULL,
              help = "Optional substring filter for pestPG filename, e.g. 'win50000.step10000'")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(
  !is.null(opt$karyotype), !is.null(opt$`pestpg-dir`),
  !is.null(opt$inversion), !is.null(opt$chr),
  !is.null(opt$start),     !is.null(opt$end), !is.null(opt$out)
)

kary <- fread(opt$karyotype, sep = "\t", header = TRUE)
need_kary <- c("sample", "inversion_id", "group")
miss_k <- setdiff(need_kary, colnames(kary))
if (length(miss_k) > 0L) {
  stop(sprintf("karyotype table missing columns: %s",
               paste(miss_k, collapse = ", ")))
}
kary <- kary[inversion_id == opt$inversion]
if (nrow(kary) == 0L) stop("No karyotype calls for this inversion")

# Locate per-sample pestPG file: <sample>.*.pestPG
find_pestpg <- function(sample) {
  files <- list.files(opt$`pestpg-dir`,
                      pattern = paste0("^", sample, "\\..*pestPG$"),
                      full.names = TRUE)
  if (!is.null(opt$`win-pattern`) && length(files) > 1L) {
    f2 <- files[grepl(opt$`win-pattern`, files, fixed = TRUE)]
    if (length(f2) >= 1L) files <- f2
  }
  if (length(files) == 0L) return(NA_character_)
  files[1]
}

# pestPG header (per ANGSD doc):
#   #(indexStart,indexStop)(posStart,posStop)(regStart,regStop) Chr WinCenter ...
# Columns we need: Chr, WinCenter, tP, nSites
read_pestpg <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  # First column has a leading "#" but pestPG files include a header line
  # like "#(...) Chr WinCenter tW tP Tajima ... nSites".
  hdr <- readLines(path, n = 1L)
  hdr <- sub("^#", "", hdr)
  cols <- strsplit(trimws(hdr), "\\s+")[[1]]
  d <- fread(path, sep = "\t", header = FALSE, skip = 1L)
  if (ncol(d) != length(cols)) {
    # Some ANGSD versions emit different separators; try whitespace.
    d <- fread(path, header = FALSE, skip = 1L, sep = " ")
  }
  if (ncol(d) != length(cols)) {
    warning(sprintf("pestPG column-count mismatch in %s (got %d, header %d)",
                    path, ncol(d), length(cols)))
    return(NULL)
  }
  setnames(d, cols)
  needed <- c("Chr", "WinCenter", "tP", "nSites")
  if (any(!needed %in% colnames(d))) {
    warning(sprintf("pestPG missing required cols in %s: have %s",
                    path, paste(colnames(d), collapse = ",")))
    return(NULL)
  }
  d[, .(Chr, WinCenter, tP, nSites)]
}

# Determine window size from first non-empty pestPG (for WinStart/WinStop)
all_samples <- unique(kary$sample)
sample_files <- sapply(all_samples, find_pestpg)
have_files <- !is.na(sample_files)
if (sum(have_files) == 0L) {
  stop("No pestPG files found for any sample in karyotype table")
}
message(sprintf("[mean_tP] pestPG files found for %d/%d samples",
                sum(have_files), length(all_samples)))

# Read all pestPG files, restrict to interval, attach group
rows <- list()
for (i in seq_along(all_samples)) {
  s <- all_samples[i]
  if (!have_files[i]) next
  d <- read_pestpg(sample_files[i])
  if (is.null(d)) next
  d <- d[Chr == opt$chr & WinCenter >= opt$start & WinCenter <= opt$end]
  if (nrow(d) == 0L) next
  d[, sample := s]
  d[, group  := kary[sample == s, group][1]]
  rows[[s]] <- d
}
if (length(rows) == 0L) stop("No pestPG rows in interval for any sample")
all <- rbindlist(rows, use.names = TRUE)

# pestPG tP is per-window SUM, not per-site density. We report tP / nSites
# (per-site tP) when nSites > 0, otherwise NA.
all[, tP_persite := ifelse(nSites > 0, tP / nSites, NA_real_)]

# Aggregate by (group, WinCenter)
sem <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x) / sqrt(length(x))
}
out <- all[, .(
  n_samples = sum(!is.na(tP_persite)),
  mean_tP   = mean(tP_persite,   na.rm = TRUE),
  median_tP = median(tP_persite, na.rm = TRUE),
  sd_tP     = sd(tP_persite,     na.rm = TRUE),
  se_tP     = sem(tP_persite)
), by = .(Chr, WinCenter, group)]

# Reconstruct WinStart/WinStop from WinCenter spacing if possible
centers <- sort(unique(out$WinCenter))
if (length(centers) >= 2L) {
  step <- centers[2L] - centers[1L]
  half <- step %/% 2L
  out[, WinStart := WinCenter - half]
  out[, WinStop  := WinCenter + half - 1L]
} else {
  out[, WinStart := NA_integer_]
  out[, WinStop  := NA_integer_]
}

out[, inversion_id := opt$inversion]
setnames(out, "Chr", "chromo")
setcolorder(out, c("inversion_id", "chromo", "WinStart", "WinStop", "WinCenter",
                   "group", "n_samples",
                   "mean_tP", "median_tP", "sd_tP", "se_tP"))

# Stable group ordering for consumers
out[, group := factor(group, levels = c("Homo_1", "Het", "Homo_2"))]
setorder(out, WinCenter, group)

fwrite(out, opt$out, sep = "\t", na = "NA", quote = FALSE)
message(sprintf("[mean_tP] wrote %d rows -> %s", nrow(out), opt$out))
