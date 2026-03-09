# =============================================================================
# Example: Processing Omniscope TCR sequencing data
# Author: Rémy Pétremand
#
# This script demonstrates the full processing pipeline from raw Omniscope
# TCR output files to a clonotype-wise summary table.
#
# Replace the file paths below with your own Omniscope output files.
# =============================================================================

# Source the pipeline functions (run from the repo root)
source("R/omniscope_tcr_r.R")

# -----------------------------------------------------------------------------
# 1. Set file paths
#    Inputs are CSVs exported from an Omniscope TCR run.
#    Each argument also accepts a directory path — CSV files will be
#    auto-discovered recursively.
# -----------------------------------------------------------------------------

file.name.contigs <- "path/to/contigs.csv"   # cell-level contig table
file.name.counts  <- "path/to/counts.csv"    # contig-level UMI counts
file.name.pairs   <- "path/to/pairs.csv"     # alpha-beta chain pairing

# -----------------------------------------------------------------------------
# 2. Load and QC-filter TCR data
#
#    Returns a list with three tables:
#      $cellwise   — one row per cell × chain
#      $contigwise — one row per unique contig_id (with counts and pairing)
#      $summary    — QC counts (raw, filtered, merged)
# -----------------------------------------------------------------------------

TCR.table.list <- load_omniscope_tcr(
  input.contigs           = file.name.contigs,
  input.counts            = file.name.counts,
  input.pairs             = file.name.pairs,
  filter_nonproductive    = TRUE,
  filter_stop_codon       = TRUE,
  filter_out_of_frame     = TRUE,
  require_complete_vdj    = FALSE,
  remove_multi_chain      = FALSE,
  resolve_ambiguous_genes = "none",
  verbose                 = TRUE
)

# Retrieve the three output tables
df.cellwise   <- TCR.table.list$cellwise    # one row per cell × chain
df.contigwise <- TCR.table.list$contigwise  # one row per unique contig
df.qc         <- TCR.table.list$summary     # QC summary

# Quick look
head(df.cellwise)
head(df.contigwise)
print(df.qc)

# -----------------------------------------------------------------------------
# 3. Assign clonotypes
#
#    Adds clonotype_id and clone_size columns to the contig-wise table.
#    Available methods: "cdr3_aa" (default), "cdr3_nt", "gene+cdr3_aa"
# -----------------------------------------------------------------------------

df.contigwise <- assign_clonotypes(tcr_df = df.contigwise, method = "cdr3_aa")

# -----------------------------------------------------------------------------
# 4. Build a clonotype-wise summary table
#
#    Collapses the contig-wise table to one row per clonotype.
#    sample_id_columns: extra metadata columns to group by (adjust to your data).
# -----------------------------------------------------------------------------

df.TCR <- get_clonotype_table(
  tcr_df            = df.contigwise,
  sample_id_columns = c("SampleID", "Patient", "Time", "TIME_POINT_SAMPLE", "chain")
)

head(df.TCR)
