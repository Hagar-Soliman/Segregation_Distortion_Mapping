# ---- Parameters (change these when adapting for a different map) ------------
map        <- "8C"
n_indiv    <- 481      # total individuals in this map
NN_thresh           <- 0.90   # exclude windows where >90% of individuals are NN
min_snps            <- 3      # exclude windows with median NumSNPs < 3 across individuals
recurrent_class     <- "AA"   # "BB" for maps BC to IM62; "AA" for maps BC to CCC9
NN_indiv_thresh     <- 0.80   # exclude individuals with >=80% NN of total windows
recur_hom_thresh    <- 0.70   # exclude individuals with >=70% recurrent-parent hom of called windows
het_thresh          <- 0.70   # exclude individuals with >=70% het of called windows
nonrecur_hom_thresh <- 0.10   # exclude individuals with >=10% non-recurrent-parent hom of called windows

# Derived threshold (number of NN individuals corresponding to NN_thresh)
NN_cutoff  <- floor(NN_thresh * n_indiv)

# ---- Paths ------------------------------------------------------------------
path         <- "/home/hks25/palmer_scratch/googaV3/refAlt/8C_refAlt/50kb_windows/"
output.path  <- "/home/hks25/palmer_scratch/googaV3/refAlt/8C_refAlt/50kb_windows/8C_g_files/"
dir.create(output.path, showWarnings = FALSE, recursive = TRUE)

setwd(path)

# ---- Load per-individual genotype summary (Genostats) ----------------------
files <- dir(path, pattern = "Genostats\\.")
data  <- data.frame(matrix(ncol = 5, nrow = length(files)))
for (k in 1:length(files)) {
  data[k, 1:5] <- read.delim(files[k], header = TRUE)
}
names(data) <- c("indiv", "AA", "AB", "BB", "NN")

# ---- Load per-window genotype counts (site.counts from step9A) -------------
site.counts <- read.delim(
  paste0(path, map, "_site.counts_hetdev=0.2+WS50K.txt"),
  header = TRUE
)
# Columns: scaffold, pos, AA, AB, BB, NN
n_windows <- nrow(site.counts)

# ---- Diagnostic counts (informational only, not used downstream) -----------

# How many windows have zero coverage across all individuals?
count <- 0
for (i in 1:length(site.counts[, 1])) {
  if (site.counts[i, 6] == n_indiv) count <- count + 1
}
cat(sprintf("Diagnostic: %d windows with no coverage in any individual\n", count))

# How many individuals have no coverage across all windows?
test <- read.delim(dir(path, pattern = "Genotypes\\.")[1], header = TRUE)
count <- 0
for (i in 1:length(data[, 1])) {
  if (data[i, 5] == nrow(test)) count <- count + 1
}
cat(sprintf("Diagnostic: %d individuals with no coverage across all windows\n", count))

# ============================================================================
# FILTER 1: High-NN windows
# Remove windows where more than NN_thresh of individuals have NN.
# ============================================================================
site.list <- data.frame(scaffold = character(), pos = integer(),
                        stringsAsFactors = FALSE)

for (i in 1:n_windows) {
  if (site.counts[i, 6] > NN_cutoff) {
    site.list <- rbind(site.list, site.counts[i, 1:2])
  }
}

cat(sprintf("Filter 1 (>%.0f%% NN): %d windows flagged\n",
            NN_thresh * 100, nrow(site.list)))

# ============================================================================
# FILTER 2: Low marker count windows
# Compute median NumSNPs per window across all individuals using window.*.txt
# files (produced by step8_window.py). Each file has columns:
#   Scaffold, WindowStart, NumSNPs, REF_Reads, ALT_Reads
# Exclude windows where the median NumSNPs across individuals < min_snps.
# ============================================================================
window.files <- dir(path, pattern = "^window\\.")

if (length(window.files) == 0) {
  warning("No window.*.txt files found in path. Skipping filter 2 (min marker count).")
} else {
  wtest      <- read.delim(window.files[1], header = TRUE)
  snp.matrix <- matrix(NA, nrow = nrow(wtest), ncol = length(window.files))

  for (k in 1:length(window.files)) {
    w <- read.delim(window.files[k], header = TRUE)
    snp.matrix[, k] <- w$NumSNPs
  }

  snp.median   <- apply(snp.matrix, 1, median, na.rm = TRUE)
  window.coords <- wtest[, 1:2]
  names(window.coords) <- c("scaffold", "pos")

  low.snp.windows <- window.coords[snp.median < min_snps, ]

  cat(sprintf("Filter 2 (<%.0f median SNPs): %d windows flagged\n",
              min_snps, nrow(low.snp.windows)))

  combined.keys    <- paste(site.list$scaffold,       site.list$pos,       sep = "_")
  lowsnp.keys      <- paste(low.snp.windows$scaffold, low.snp.windows$pos, sep = "_")
  new.low.snp      <- low.snp.windows[!(lowsnp.keys %in% combined.keys), ]
  site.list        <- rbind(site.list, new.low.snp)
}

cat(sprintf("Total windows in exclusion list (both filters): %d\n", nrow(site.list)))

site.list.file <- paste0(path, map, "_site_list_NN", NN_thresh*100, "pct_minSNP", min_snps, ".txt")
write.table(site.list, file = site.list.file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Site exclusion list written to: %s\n", site.list.file))

# ============================================================================
# FILTER 3: Individual quality filters
# ============================================================================
indiv.list    <- c()
indiv.reasons <- c()

for (i in 1:nrow(data)) {
  aa     <- data[i, 2]
  ab     <- data[i, 3]
  bb     <- data[i, 4]
  nn     <- data[i, 5]
  total  <- aa + ab + bb + nn
  called <- aa + ab + bb

  if (recurrent_class == "AA") {
    recur_hom    <- aa
    nonrecur_hom <- bb
  } else {
    recur_hom    <- bb
    nonrecur_hom <- aa
  }

  reason <- NA
  if      (called == 0)                                               reason <- "no_called_windows"
  else if (total  > 0 && nn / total            >= NN_indiv_thresh)   reason <- sprintf(">=80pct_NN (%.1f%%)",             100 * nn / total)
  else if (called > 0 && recur_hom / called    >= recur_hom_thresh)  reason <- sprintf(">=70pct_recurrent_hom (%.1f%%)", 100 * recur_hom / called)
  else if (called > 0 && ab / called           >= het_thresh)        reason <- sprintf(">=70pct_het (%.1f%%)",            100 * ab / called)
  else if (called > 0 && nonrecur_hom / called >= nonrecur_hom_thresh) reason <- sprintf(">=10pct_nonrecurrent_hom (%.1f%%)", 100 * nonrecur_hom / called)

  if (!is.na(reason)) {
    indiv.list    <- append(indiv.list,    data[i, 1])
    indiv.reasons <- append(indiv.reasons, reason)
  }
}

cat(sprintf("Filter 3 (individual quality): %d individuals flagged\n", length(indiv.list)))
if (length(indiv.list) > 0) {
  for (j in 1:length(indiv.list)) {
    cat(sprintf("  %s: %s\n", indiv.list[j], indiv.reasons[j]))
  }
}

indiv.list.file <- paste0(path, map, "_indiv_list_filtered.txt")
write.table(indiv.list, file = indiv.list.file, sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# Generate g. files
# ============================================================================
genotype.files <- dir(path, pattern = "Genotypes\\.")

indiv.list.loaded <- scan(indiv.list.file, what = "", sep = "\t")
indiv.list.loaded <- paste0("Genotypes.", indiv.list.loaded, ".txt")

site.list.loaded  <- read.delim(site.list.file, header = TRUE)
remove_keys       <- paste(site.list.loaded$scaffold, site.list.loaded$pos, sep = "_")

n_written <- 0
for (i in 1:length(genotype.files)) {

  if (genotype.files[i] %in% indiv.list.loaded) next

  temp      <- read.delim(genotype.files[i], header = TRUE)
  temp_keys <- paste(temp$Scaffold, temp$WindowStart, sep = "_")
  keep_rows <- !(temp_keys %in% remove_keys)
  genotypes <- temp[keep_rows, ]

  write.table(genotypes,
              file      = paste0(output.path, "g.", genotype.files[i]),
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  n_written <- n_written + 1
}

cat(sprintf("g. files written for %d individuals to: %s\n", n_written, output.path))
