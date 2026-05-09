#!/usr/bin/env Rscript
# =============================================================================
# make_karyotype_bamlists.R
# =============================================================================
# Build Homo_1 / Homo_2 / Het bamlists for one inversion, by joining the
# karyotype calls table to the sample manifest.
#
# Inputs
#   --karyotype  TSV with columns: sample, inversion_id, group
#                  group must be one of: Homo_1, Het, Homo_2
#   --manifest   TSV with at least: Sample, <BAM cols...>
#                  Column 3 (1-based) is the filtered BAM, matching
#                  the existing convention (see STEP_A02).
#   --inversion  inversion_id to extract
#
# Outputs (one BAM path per line)
#   --out-homo1, --out-homo2, --out-het
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
})

opt_list <- list(
  make_option("--karyotype",  type = "character"),
  make_option("--manifest",   type = "character"),
  make_option("--inversion",  type = "character"),
  make_option("--out-homo1",  type = "character"),
  make_option("--out-homo2",  type = "character"),
  make_option("--out-het",    type = "character"),
  make_option("--bam-col",    type = "integer", default = 3L,
              help = "1-based column index of the filtered BAM in manifest [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(
  !is.null(opt$karyotype), !is.null(opt$manifest), !is.null(opt$inversion),
  !is.null(opt$`out-homo1`), !is.null(opt$`out-homo2`), !is.null(opt$`out-het`)
)

kary <- read.table(opt$karyotype, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
need_kary <- c("sample", "inversion_id", "group")
miss_k <- setdiff(need_kary, colnames(kary))
if (length(miss_k) > 0L) {
  stop(sprintf("karyotype table missing required columns: %s",
               paste(miss_k, collapse = ", ")))
}
kary <- kary[kary$inversion_id == opt$inversion, , drop = FALSE]
if (nrow(kary) == 0L) {
  stop(sprintf("No karyotype calls for inversion '%s' in %s",
               opt$inversion, opt$karyotype))
}

# Manifest: tab-separated, header row, sample in col 1, BAM in col --bam-col
man <- read.table(opt$manifest, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE,
                  comment.char = "")
sample_col <- colnames(man)[1L]
bam_col    <- colnames(man)[opt$`bam-col`]
if (is.na(bam_col)) {
  stop(sprintf("manifest has %d columns; --bam-col=%d is out of range",
               ncol(man), opt$`bam-col`))
}

# Build sample -> BAM lookup; drop empty BAM cells
man <- man[!is.na(man[[bam_col]]) & nzchar(man[[bam_col]]), , drop = FALSE]
sample2bam <- setNames(man[[bam_col]], man[[sample_col]])

write_band <- function(group_label, outfile) {
  s <- kary$sample[kary$group == group_label]
  bams <- sample2bam[s]
  missing_bam <- s[is.na(bams)]
  if (length(missing_bam) > 0L) {
    message(sprintf("[make_karyotype_bamlists] %s: %d sample(s) missing in manifest: %s",
                    group_label, length(missing_bam),
                    paste(head(missing_bam, 5), collapse = ", ")))
  }
  bams <- bams[!is.na(bams)]
  writeLines(bams, outfile)
  message(sprintf("[make_karyotype_bamlists] %s: wrote %d BAMs -> %s",
                  group_label, length(bams), outfile))
}

write_band("Homo_1", opt$`out-homo1`)
write_band("Homo_2", opt$`out-homo2`)
write_band("Het",    opt$`out-het`)
