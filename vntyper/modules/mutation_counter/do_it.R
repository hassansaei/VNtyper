#!/usr/bin/env Rscript
version <- 0.1
## Author: Martin Figeac / Adapted for our CSV file
## 
## Usage:
##   Rscript do_it.R count.csv ALL MUT1 {MUT2 ...}
##
##   where:
##     - count.csv is the CSV file produced by count_variants.sh (tab-delimited)
##     - The second argument is the control column (should be "All")
##     - The remaining arguments are the mutation columns to test.
##
## This script computes pairwise Fisher’s exact tests on the counts (CTL - MUT)
## and prints a matrix of p-values along with a summary of samples with significant p-values.
##
## Note: At least one count (either control or mutation) for each sample must be >= MIN_COUNTS.
MIN_COUNTS <- 100

cat("## script do_it.R running in version", version, "\n")
cat("## date is", date(), "\n")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  cat("Erreur usage, need at least 3 arguments: Rscript do_it.R count.csv CTL MUT1 {MUT2 ...}\n")
  quit(status = 1)
}

cat("Opening file", args[1], "...\n")
# Read the CSV file; adjust sep if needed (our count file is tab-delimited)
csv <- read.csv(args[1], sep="\t", header=TRUE, as.is=TRUE, comment.char="#")
cat("done\n")

# Check that the control column (args[2]) exists
if(!(args[2] %in% colnames(csv))) {
  cat("Erreur argument for control:", args[2], "not found in count file", args[1], "\n")
  cat("colnames =", colnames(csv), "\n")
  quit(status = 1)
}
ctl <- which(colnames(csv) == args[2])

# Loop over each mutation argument (from the third argument onward)
for(mutationsArg in 3:length(args)) {
  if(!(args[mutationsArg] %in% colnames(csv))) {
    cat("Erreur argument for mutation:", args[mutationsArg], "not found in count file", args[1], "\n")
    cat("colnames =", colnames(csv), "\n")
    quit(status = 1)
  }
  mut <- which(colnames(csv) == args[mutationsArg])
  cat("\n##\n")
  cat("## Testing mutation", args[mutationsArg], "against control", args[2], "\n")
  
  # Select rows (samples) that have at least MIN_COUNTS in one of the two columns (CTL or MUT)
  rows <- c()
  for (sample in 1:nrow(csv)) {
    if (any(csv[sample, c(mut, ctl)] >= MIN_COUNTS)) {
      rows <- c(rows, sample)
    } else {
      cat("WARNING: Sample", csv[sample, 1], "doesn't have enough counts (", MIN_COUNTS, ") and has been excluded\n", sep="")
    }
  }
  
  if(length(rows) < 2) {
    cat("Not enough samples with sufficient counts to compare.\n")
    next
  }
  
  cat("\nFisher Matrix:\n")
  
  # Create an empty matrix to hold the p-values for pairwise comparisons.
  pMatrix <- matrix(0, nrow=length(rows), ncol=length(rows))
  sampleNames <- substr(csv[rows, 1], 1, 20)
  colnames(pMatrix) <- sampleNames
  rownames(pMatrix) <- sampleNames
  
  i <- 1
  for (sample1 in rows) {
    cat(substr(csv[sample1, 1], 1, 20), "\t")
    j <- 1
    for (sample2 in rows) {
      # Build the 2x2 contingency table:
      #   [ (CTL - MUT) for sample2, (CTL - MUT) for sample1 ]
      #   [ MUT for sample2,         MUT for sample1 ]
      m <- matrix(c(
        csv[sample2, ctl] - csv[sample2, mut],
        csv[sample1, ctl] - csv[sample1, mut],
        csv[sample2, mut],
        csv[sample1, mut]
      ), nrow = 2)
      
      # Perform Fisher's exact test with alternative "greater"
      p <- fisher.test(m, alternative = "greater")$p.value
      pMatrix[i, j] <- p
      cat(p, "\t")
      j <- j + 1
    }
    cat("\n")
    i <- i + 1
  }
  
  cat("\n")
  # Set the diagonal to NA
  for (sample in 1:nrow(pMatrix)) {
    pMatrix[sample, sample] <- NA
  }
  
  # Report samples with significant p-values at different thresholds.
  for (pos in 1:min(6, nrow(pMatrix)-1)) {
    if (pos == 1) {
      cat("Who is significant with only", pos, "positive sample in this library (p<=0.001):")
    } else {
      cat("Who is significant with only", pos, "positives samples in this library (p<=0.001):")
    }
    pass <- 0
    for (i in 1:nrow(pMatrix)) {
      id <- order(pMatrix[i, ], decreasing = TRUE)
      if (pMatrix[i, id[pos]] <= 0.001) {
        cat("\n", rownames(pMatrix)[i], "with a p-value of", pMatrix[i, id[pos]])
        pass <- pass + 1
      }
    }
    if (pass == 0) {
      cat(" none!\n")
    } else {
      cat("\n")
    }
    cat("\n")
  }
}
