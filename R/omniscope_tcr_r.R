# =============================================================================
# LoadOmniscopeData.R
# Functions to load, QC-filter, and format Omniscope OS-T TCR data.
#
# Three granularities are produced:
#   - cell-wise      : one row per cell × chain (from input.contigs)
#   - contig-wise    : one row per unique contig_id, augmented with UMI counts (input.counts) and pairing info (input.pairs)
#   - clonotype-wise : one row per unique clontoype defined from contig-wise table
#
# =============================================================================


# -----------------------------------------------------------------------------
# Private helpers
# -----------------------------------------------------------------------------

# Load input: accepts a pre-loaded data.frame, file path(s), or a directory.
.load_csv_input <- function(x, name, verbose) {
  if (is.data.frame(x)) {
    if (verbose && nrow(x) > 0)
      message("  ", name, ": using pre-loaded data.frame (", nrow(x), " rows).")
    return(x)
  }
  if (length(x) == 1 && dir.exists(x)) {
    files <- list.files(x, pattern = "\\.(csv|CSV)$", full.names = TRUE, recursive = TRUE)
    if (length(files) == 0) { warning("No CSV files found in directory: ", x); return(data.frame()) }
    if (verbose) message("  ", name, ": loading ", length(files), " file(s) from directory.")
    return(do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE)))
  }
  if (is.character(x)) {
    missing_files <- x[!file.exists(x)]
    if (length(missing_files) > 0) warning("File(s) not found: ", paste(missing_files, collapse = ", "))
    found <- x[file.exists(x)]
    if (length(found) == 0) return(data.frame())
    if (verbose) message("  ", name, ": loading ", length(found), " file(s).")
    return(do.call(rbind, lapply(found, read.csv, stringsAsFactors = FALSE)))
  }
  stop(name, " must be a data.frame, file path(s), or directory path.")
}

# Stop with a clear message if required columns are absent.
.validate_cols <- function(df, required, name) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0)
    stop("Missing required columns in ", name, ": ", paste(missing, collapse = ", "),
         "\nAvailable: ", paste(colnames(df), collapse = ", "))
  invisible(NULL)
}

# Coerce string or logical values to proper logical type.
.to_logical <- function(x) {
  if (is.logical(x)) return(x)
  x <- trimws(toupper(as.character(x)))
  ifelse(x %in% c("TRUE", "T"), TRUE, ifelse(x %in% c("FALSE", "F"), FALSE, NA))
}

# Resolve ambiguous gene calls (e.g. "TRBV7-9*01|TRBV7-9*02" -> "TRBV7-9*01").
.resolve_gene <- function(gene_str, method = "first") {
  if (is.na(gene_str) || gene_str == "") return(NA_character_)
  if (method == "first") return(strsplit(gene_str, "\\|")[[1]][1])
  gene_str  # method = "all": keep as-is
}


# =============================================================================
# 1. Load and QC-filter contigs — cell-wise output
# =============================================================================

#' Load and QC-filter Omniscope contigs table (cell-wise)
#'
#' Reads one or more Omniscope contigs CSV files, applies QC filters, resolves
#' ambiguous gene calls, optionally removes extra chains per cell, and returns
#' a standardised cell-wise data.frame (one row per cell x chain).
#'
#' @param input.contigs data.frame, character vector of file paths, or a single
#'   directory path. Must contain: cell_id, contig_id, locus, v_call, j_call,
#'   cdr3, cdr3_aa, productive.
#' @param filter_nonproductive Logical. Remove non-productive contigs (default TRUE).
#' @param filter_stop_codon Logical. Remove stop-codon-containing contigs (default TRUE).
#' @param filter_out_of_frame Logical. Remove out-of-frame contigs (default TRUE).
#' @param require_complete_vdj Logical. Require complete VDJ annotation (default FALSE,
#'   because many valid TRB contigs lack a D-gene assignment).
#' @param remove_multi_chain Logical. Keep at most 2 chains per cell, ranked by
#'   UMI count (default TRUE). Requires umi_count column.
#' @param resolve_ambiguous_genes Character. "first" (default) keeps the first
#'   allele before a pipe; "all" keeps the raw string.
#' @param verbose Logical. Print progress and QC summary (default TRUE).
#'
#' @return data.frame with columns:
#'   barcode, contig_id, sample_id, chain, v_gene, d_gene, j_gene,
#'   cdr3_nt, cdr3_aa, junction, junction_aa, umi_count, productive,
#'   raw_v_call, raw_d_call, raw_j_call,
#'   plus any extra metadata columns present in the input (Patient, Time, etc.).

load_omniscope_contigs_cellwise <- function(
    input.contigs,
    filter_nonproductive    = TRUE,
    filter_stop_codon       = TRUE,
    filter_out_of_frame     = TRUE,
    require_complete_vdj    = FALSE,
    remove_multi_chain      = TRUE,
    resolve_ambiguous_genes = "first",
    verbose                 = TRUE
) {
  if (verbose) message("--- Loading contigs (cell-wise) ---")
  df <- .load_csv_input(input.contigs, "contigs", verbose)
  if (nrow(df) == 0) { warning("No contigs loaded."); return(data.frame()) }

  .validate_cols(df, c("cell_id", "contig_id", "locus", "v_call", "j_call",
                        "cdr3", "cdr3_aa", "productive"), "contigs")

  n_start <- nrow(df)
  if (verbose) message("  Starting: ", n_start, " contigs | ",
                       length(unique(df$cell_id)), " cells")

  # Standardise logical columns (handle TRUE/FALSE stored as strings)
  for (col in intersect(c("productive", "stop_codon", "vj_in_frame", "complete_vdj"), colnames(df)))
    df[[col]] <- .to_logical(df[[col]])

  # Replace empty strings with NA
  df[df == ""] <- NA

  # ---------------------------------------------------------------------------
  # QC filtering — each step is logged
  # ---------------------------------------------------------------------------
  filter_log <- character(0)

  .filt <- function(df, step, mask) {
    n_before <- nrow(df)
    df <- df[mask, ]
    filter_log[[step]] <<- sprintf("-%d (remaining: %d)", n_before - nrow(df), nrow(df))
    df
  }

  if (filter_nonproductive && "productive" %in% colnames(df))
    df <- .filt(df, "Remove non-productive",  !is.na(df$productive)  & df$productive)

  if (filter_stop_codon && "stop_codon" %in% colnames(df))
    df <- .filt(df, "Remove stop codons",     !is.na(df$stop_codon)  & !df$stop_codon)

  if (filter_out_of_frame && "vj_in_frame" %in% colnames(df))
    df <- .filt(df, "Remove out-of-frame",    !is.na(df$vj_in_frame) & df$vj_in_frame)

  if (require_complete_vdj && "complete_vdj" %in% colnames(df))
    df <- .filt(df, "Require complete VDJ",   !is.na(df$complete_vdj) & df$complete_vdj)

  if ("cdr3" %in% colnames(df))
    df <- .filt(df, "Remove missing CDR3",    !is.na(df$cdr3) & nchar(df$cdr3) > 0)

  if (all(c("v_call", "j_call") %in% colnames(df)))
    df <- .filt(df, "Remove missing V/J",     !is.na(df$v_call) & !is.na(df$j_call))

  if (nrow(df) == 0) { warning("No contigs remaining after QC filtering!"); return(data.frame()) }

  # ---------------------------------------------------------------------------
  # Resolve ambiguous gene calls — save raw versions first
  # ---------------------------------------------------------------------------
  for (col in intersect(c("v_call", "d_call", "j_call"), colnames(df))) {
    df[[paste0("raw_", col)]] <- df[[col]]
    df[[col]] <- sapply(df[[col]], .resolve_gene, method = resolve_ambiguous_genes,
                        USE.NAMES = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Multi-chain removal — keep top 2 chains per cell, ranked by UMI count
  # ---------------------------------------------------------------------------
  if (remove_multi_chain && "cell_id" %in% colnames(df) && "umi_count" %in% colnames(df)) {
    n_before <- nrow(df)
    df <- df %>%
      dplyr::arrange(cell_id, dplyr::desc(umi_count)) %>%
      dplyr::group_by(cell_id) %>%
      dplyr::slice_head(n = 2) %>%
      dplyr::ungroup() %>%
      as.data.frame()
    filter_log[["Keep top 2 chains per cell"]] <-
      sprintf("-%d (remaining: %d)", n_before - nrow(df), nrow(df))
  }

  # ---------------------------------------------------------------------------
  # Standardised output
  # ---------------------------------------------------------------------------
  internal_cols <- c("cell_id", "contig_id", "SampleID", "locus",
                     "v_call", "d_call", "j_call", "cdr3", "cdr3_aa",
                     "junction", "junction_aa", "umi_count", "productive",
                     "raw_v_call", "raw_d_call", "raw_j_call",
                     "stop_codon", "vj_in_frame", "complete_vdj")
  extra_cols <- setdiff(colnames(df), internal_cols)

  out <- data.frame(
    barcode     = df$cell_id,
    contig_id   = df$contig_id,
    sample_id   = if ("SampleID"    %in% colnames(df)) df$SampleID    else NA_character_,
    chain       = df$locus,
    v_gene      = if ("v_call"      %in% colnames(df)) df$v_call      else NA_character_,
    d_gene      = if ("d_call"      %in% colnames(df)) df$d_call      else NA_character_,
    j_gene      = if ("j_call"      %in% colnames(df)) df$j_call      else NA_character_,
    cdr3_nt     = if ("cdr3"        %in% colnames(df)) df$cdr3        else NA_character_,
    cdr3_aa     = if ("cdr3_aa"     %in% colnames(df)) df$cdr3_aa     else NA_character_,
    junction    = if ("junction"    %in% colnames(df)) df$junction    else NA_character_,
    junction_aa = if ("junction_aa" %in% colnames(df)) df$junction_aa else NA_character_,
    umi_count   = if ("umi_count"   %in% colnames(df)) df$umi_count   else NA_real_,
    productive  = if ("productive"  %in% colnames(df)) df$productive  else NA,
    raw_v_call  = if ("raw_v_call"  %in% colnames(df)) df$raw_v_call  else NA_character_,
    raw_d_call  = if ("raw_d_call"  %in% colnames(df)) df$raw_d_call  else NA_character_,
    raw_j_call  = if ("raw_j_call"  %in% colnames(df)) df$raw_j_call  else NA_character_,
    stringsAsFactors = FALSE
  )

  # Append extra metadata columns (Patient, Time, etc.)
  if (length(extra_cols) > 0)
    out <- cbind(out, df[, extra_cols, drop = FALSE])

  if (verbose) {
    message("\n  QC Filtering:")
    for (step in names(filter_log))
      message(sprintf("    %-35s %s", step, filter_log[[step]]))
    message(sprintf("  Final: %d contigs | %d cells | %.1f%% retained",
                    nrow(out), length(unique(out$barcode)), 100 * nrow(out) / n_start))
    chains <- table(out$chain)
    message(sprintf("  Chains: %s",
                    paste(paste0(names(chains), " (n=", chains, ")"), collapse = ", ")))
  }

  attr(out, "n_raw") <- n_start
  out
}


# =============================================================================
# 2. Collapse cell-wise contigs to a contig-wise table
# =============================================================================

#' Collapse a cell-wise contigs table to one row per unique contig_id
#'
#' Groups by contig_id and summarises:
#'   - sequence-level columns (chain, genes, CDR3, ...) are taken from the first
#'     row in each group (they should be identical across cells sharing a contig);
#'   - n_cells and cell_ids record how many and which cells carry each contig;
#'   - umi_count is averaged across cells.
#'
#' Extra metadata columns (Patient, Time, etc.) are also taken from the first row.
#'
#' @param df.cellwise Output of \code{load_omniscope_contigs_cellwise()}.
#' @param verbose Logical (default TRUE).
#'
#' @return data.frame with one row per unique contig_id. Columns:
#'   contig_id, [sequence cols], n_cells, cell_ids, umi_count (mean across cells).

contigs_to_contigwise <- function(df.cellwise, verbose = TRUE) {
  if (nrow(df.cellwise) == 0) return(data.frame())
  if (verbose) message("--- Collapsing to contig-wise table ---")

  # Part 1: column handling: 
  
  # Adding umi_count if missing
  if (!"umi_count" %in% colnames(df.cellwise)){
    df.cellwise$umi_count <- 0
  }
  
  # Columns describing the TCR sequence — expected to be constant per contig_id
  # Note: There are some small inconsistencies:
  # Most inconsistencies are due to missing data (NA values)
  # Very rare inconsistencies were found in the d_gene with some missing data
  # Some inconsistencies were found in the "junction" and  "junction_aa" columns
  # For the sake of simplicity, this is not processed here. 
  seq_cols <- intersect(
    c("v_gene", "d_gene", "j_gene",
      "cdr3_nt", "cdr3_aa", 
      "junction", "junction_aa",
      "raw_v_call", "raw_d_call", "raw_j_call"),
    colnames(df.cellwise)
  )
  
  # df.sanity.check <- df.cellwise %>% 
  #   dplyr::group_by(contig_id) %>%
  #   dplyr::summarise_at(
  #     .vars = dplyr::all_of(seq_cols), 
  #     .funs = dplyr::n_distinct, 
  #     na.rm = TRUE
  #   ) %>%
  #   as.data.frame()
  # contig_issue <- unique(df.sanity.check[apply(df.sanity.check[, seq_cols], 1, function(x) max(x, na.rm = T)) > 1, ]$contig_id)
  
  # Get groups:
  df.contigwise <- df.cellwise %>%
    dplyr::group_by(sample_id, chain, contig_id) %>%
    dplyr::summarise(
      n_cells  = dplyr::n(),
      cell_ids = paste(barcode, collapse = ","),
      umi_count = mean(umi_count, na.rm = T),
      dplyr::across(dplyr::all_of(seq_cols), dplyr::first),
      .groups  = "drop"
    ) %>%
    as.data.frame()
  
  if (verbose){
    message(sprintf("  %d unique contigs from %d cell-wise rows (avg %.1f cells/contig)",
                    nrow(df.contigwise), nrow(df.cellwise),
                    nrow(df.cellwise) / nrow(df.contigwise)))
  }
  
  df.contigwise
}


# =============================================================================
# 3. Load counts table (contig-wise)
# =============================================================================

#' Load and validate an Omniscope UMI/read counts table
#'
#' @param input.counts data.frame, file path(s), or directory of CSV files.
#'   Must contain: contig_id, counts.
#' @param verbose Logical (default TRUE).
#'
#' @return data.frame with at least contig_id and counts columns (plus any
#'   additional columns present in the file).

load_omniscope_counts <- function(input.counts, verbose = TRUE) {
  if (verbose) message("--- Loading counts table ---")
  df <- .load_csv_input(input.counts, "counts", verbose)
  if (nrow(df) == 0) return(data.frame())
  .validate_cols(df, c("contig_id", "counts"), "counts")
  if (verbose) message("  ", nrow(df), " contig count entries.")
  df
}


# =============================================================================
# 4. Load pairs table (contig-wise)
# =============================================================================

#' Load, validate, and reshape an Omniscope paired-chain table
#'
#' The raw pairs table has one row per cell, recording the alpha (A) and beta
#' (B) contig IDs. This function reshapes it so that contig_id_B becomes the
#' join key (renamed to contig_id), making it ready to left-join onto a
#' contig-wise table indexed on the beta chain.
#'
#' @param input.pairs data.frame, file path(s), or directory of CSV files.
#'   Must contain: contig_id_A, cdr3_aa_A, contig_id_B, cdr3_aa_B.
#' @param verbose Logical (default TRUE).
#'
#' @return data.frame with columns: contig_id (= former contig_id_B),
#'   contig_id_A, cdr3_aa_A, cdr3_aa_B, plus any other columns from the input.

load_omniscope_pairs <- function(input.pairs, verbose = TRUE) {
  if (verbose) message("--- Loading pairs table ---")
  df <- .load_csv_input(input.pairs, "pairs", verbose)
  if (nrow(df) == 0) return(data.frame())
  .validate_cols(df, c("contig_id_A", "cdr3_aa_A", "contig_id_B", "cdr3_aa_B"), "pairs")

  # contig_id_B is the join key for the beta-chain contig-wise table
  df <- dplyr::rename(df, contig_id = contig_id_B)

  if (verbose)
    message(sprintf("  %d pair entries (%d unique B-chain contigs).",
                    nrow(df), length(unique(df$contig_id))))
  df
}


# =============================================================================
# 5. Main orchestrator
# =============================================================================

#' Load and process Omniscope TCR data
#'
#' Orchestrates all loading steps and returns both granularities:
#'
#' \enumerate{
#'   \item Load and QC-filter contigs -> cell-wise table
#'   \item Collapse to contig-wise table
#'   \item Merge counts into contig-wise table (if provided)
#'   \item Merge pairs into contig-wise table (if provided)
#' }
#'
#' @param input.contigs data.frame, file path(s), or directory (cell-wise).
#'   Pass NULL to skip (default).
#' @param input.counts data.frame, file path(s), or directory (contig-wise).
#'   Pass NULL to skip (default).
#' @param input.pairs data.frame, file path(s), or directory (contig-wise).
#'   Pass NULL to skip (default).
#' @param filter_nonproductive,filter_stop_codon,filter_out_of_frame,
#'   require_complete_vdj,remove_multi_chain,resolve_ambiguous_genes,verbose
#'   Forwarded to \code{load_omniscope_contigs_cellwise()}.
#'
#' @return Named list:
#' \describe{
#'   \item{cellwise}{Cell-wise table (one row per cell x chain). Empty
#'     data.frame if no contigs were provided.}
#'   \item{contigwise}{Contig-wise table (one row per unique contig_id),
#'     augmented with counts and/or pairing information if provided.}
#' }
#'
#' @examples
#' res <- load_omniscope_tcr(
#'   input.contigs = "path/to/contigs.csv",
#'   input.counts  = "path/to/counts.csv",
#'   input.pairs   = "path/to/pairs.csv"
#' )
#' res$cellwise    # cell x chain table
#' res$contigwise  # clone-level table with counts and pairing

load_omniscope_tcr <- function(
    input.contigs           = NULL,
    input.counts            = NULL,
    input.pairs             = NULL,
    filter_nonproductive    = TRUE,
    filter_stop_codon       = TRUE,
    filter_out_of_frame     = TRUE,
    require_complete_vdj    = FALSE,
    remove_multi_chain      = TRUE,
    resolve_ambiguous_genes = "first",
    verbose                 = TRUE
) {
  
  # Create empty placeholders
  df.cellwise   <- data.frame()
  df.contigwise <- data.frame()

  # Initialise summary counters
  contigs_n_contigs_total  <- NA_integer_
  contigs_n_contigs_filter <- NA_integer_
  counts_n_contigs_total   <- NA_integer_
  counts_n_contigs_common  <- NA_integer_
  pairs_n_contigs_total    <- NA_integer_
  pairs_n_contigs_common   <- NA_integer_

  # 1. Load contigs -> cell-wise table
  if (!is.null(input.contigs) &&
      !(is.data.frame(input.contigs) && nrow(input.contigs) == 0)) {
    df.cellwise <- load_omniscope_contigs_cellwise(
      input.contigs,
      filter_nonproductive    = filter_nonproductive,
      filter_stop_codon       = filter_stop_codon,
      filter_out_of_frame     = filter_out_of_frame,
      require_complete_vdj    = require_complete_vdj,
      remove_multi_chain      = remove_multi_chain,
      resolve_ambiguous_genes = resolve_ambiguous_genes,
      verbose                 = verbose
    )
    contigs_n_contigs_total <- attr(df.cellwise, "n_raw")
  }

  # 2. Collapse cell-wise into contig-wise
  if (nrow(df.cellwise) > 0){
    df.contigwise <- contigs_to_contigwise(df.cellwise, verbose = verbose)
    contigs_n_contigs_filter <- nrow(df.contigwise)
  }

  # 3. Load counts and merge into contig-wise
  if (!is.null(input.counts) &&
      !(is.data.frame(input.counts) && nrow(input.counts) == 0)) {
    df.counts <- load_omniscope_counts(input.counts, verbose = verbose)
    if (nrow(df.counts) > 0) {
      counts_n_contigs_total <- nrow(df.counts)
      if (nrow(df.contigwise) > 0) {
        n_common <- length(
          intersect(df.contigwise$contig_id, df.counts$contig_id)
        )
        counts_n_contigs_common <- n_common
        if (verbose)
          message(sprintf(
            "  Merge counts -> contigwise: %d / %d contigwise match",
            n_common, nrow(df.contigwise)
          ), sprintf(" (%d total counts contigs).", counts_n_contigs_total))
        df.contigwise <- merge(df.contigwise, df.counts, by = "contig_id", all.x = TRUE)
      } else {
        # No contigs table: counts becomes the base contig-wise table
        df.contigwise <- df.counts
      }
    }
  }

  # 4. Load pairs and merge into contig-wise
  if (!is.null(input.pairs) &&
      !(is.data.frame(input.pairs) && nrow(input.pairs) == 0)) {
    df.pairs <- load_omniscope_pairs(input.pairs, verbose = verbose)
    if (nrow(df.pairs) > 0) {
      pairs_n_contigs_total <- nrow(df.pairs)
      if (nrow(df.contigwise) > 0) {
        n_common <- length(
          intersect(df.contigwise$contig_id, df.pairs$contig_id)
        )
        pairs_n_contigs_common <- n_common
        if (verbose)
          message(sprintf(
            "  Merge pairs -> contigwise: %d / %d contigwise match",
            n_common, nrow(df.contigwise)
          ), sprintf(" (%d total pairs contigs).", pairs_n_contigs_total))
        df.contigwise <- merge(df.contigwise, df.pairs, by = "contig_id", all.x = TRUE)
      } else {
        # No prior contig-wise table: pairs becomes the base
        df.contigwise <- df.pairs
      }
    }
  }

  # Build summary
  qc_summary <- list(
    contigs_n_contigs_total  = contigs_n_contigs_total,
    contigs_n_contigs_filter = contigs_n_contigs_filter,
    counts_n_contigs_total   = counts_n_contigs_total,
    counts_n_contigs_common  = counts_n_contigs_common,
    pairs_n_contigs_total    = pairs_n_contigs_total,
    pairs_n_contigs_common   = pairs_n_contigs_common
  )

  if (verbose) {
    message("\n===== Final Output =====")
    if (nrow(df.cellwise) > 0) {
      message(sprintf("  Cell-wise:   %d rows | %d unique cells",
                      nrow(df.cellwise), length(unique(df.cellwise$barcode))))
    } else {
      message("  Cell-wise:   (none)")
    }
    if (nrow(df.contigwise) > 0) {
      message(sprintf("  Contig-wise: %d unique contigs", nrow(df.contigwise)))
      if ("counts" %in% colnames(df.contigwise))
        message(sprintf("               counts:  median=%.0f, range=[%d-%d]",
                        median(df.contigwise$counts, na.rm = TRUE),
                        min(df.contigwise$counts, na.rm = TRUE),
                        max(df.contigwise$counts, na.rm = TRUE)))
      if ("contig_id_A" %in% colnames(df.contigwise))
        message(sprintf("               paired:  %d / %d contigs have a partner",
                        sum(!is.na(df.contigwise$contig_id_A)), nrow(df.contigwise)))
    } else {
      message("  Contig-wise: (none)")
    }
  }

  list(cellwise = df.cellwise, contigwise = df.contigwise,
       summary = qc_summary)
}


# =============================================================================
# Convenience: assign clonotype IDs
# =============================================================================

#' Assign clonotype IDs to a cell-wise or contig-wise TCR table
#'
#' @param tcr_df A data.frame with at least cdr3_aa (and optionally v_gene, j_gene).
#'   Works on either the cellwise or contigwise output of load_omniscope_tcr().
#' @param method "cdr3_aa" (default), "cdr3_nt", or "gene+cdr3_aa".
#'
#' @return The input data.frame with added columns: clonotype_id, clone_size.

assign_clonotypes <- function(tcr_df, method = "cdr3_aa") {
  tcr_df$clone_key <- switch(method,
    "cdr3_aa"      = tcr_df$cdr3_aa,
    "cdr3_nt"      = tcr_df$cdr3_nt,
    "gene+cdr3_aa" = paste(tcr_df$v_gene, tcr_df$j_gene, tcr_df$cdr3_aa, sep = "_"),
    stop("method must be 'cdr3_aa', 'cdr3_nt', or 'gene+cdr3_aa'")
  )

  freq   <- sort(table(tcr_df$clone_key), decreasing = TRUE)
  id_map <- setNames(seq_along(freq), names(freq))

  tcr_df$clonotype_id <- paste0("clonotype_", id_map[tcr_df$clone_key])
  tcr_df$clone_size   <- as.integer(freq[tcr_df$clone_key])
  tcr_df$clone_key    <- NULL

  return(tcr_df)
}


#' Collapse a contig-wise table to one row per clonotype
#'
#' Groups by \code{chain} and \code{clonotype_id} and summarises all
#' contig-level columns into clonotype-level aggregates. For numeric columns
#' (\code{n_cells}, \code{counts}) values are summed across contigs. For
#' character columns the number of unique non-NA values is recorded in a
#' paired \code{*_n} column and the unique values are collapsed to a
#' comma-separated string.
#'
#' Optional columns (\code{d_gene}, \code{cdr3_nt}, \code{junction},
#' \code{junction_aa}, \code{cdr3_aa_A}, \code{contig_id_A}, \code{cell_ids},
#' etc.) are gracefully handled: if absent from \code{tcr_df} a warning is
#' issued and the corresponding output columns are filled with \code{NA}.
#'
#' @param tcr_df A contig-wise data.frame, typically the \code{contigwise}
#'   element returned by \code{load_omniscope_tcr()} after
#'   \code{assign_clonotypes()} has been called. Must contain at least
#'   \code{chain} and \code{clonotype_id}.
#' @param sample_id_columns Character vector of additional column names to
#'   include in the \code{group_by} (e.g. \code{c("sample_id")} or
#'   \code{c("Patient", "Time")}). Default \code{NULL} groups only by
#'   \code{chain} and \code{clonotype_id}.
#'
#' @return A data.frame with one row per grouping combination. Columns mirror
#'   the input but with numeric columns summed and character columns collapsed;
#'   each character column gains a paired \code{<col>_n} column with the count
#'   of unique values.

get_clonotype_table <- function(tcr_df, sample_id_columns = NULL) {

  # Suppress R CMD CHECK "no visible binding" notes from data.table NSE
  v_gene <- d_gene <- j_gene <- cdr3_nt <- cdr3_aa <- cdr3_aa_A <-
    junction <- junction_aa <- contig_id <- contig_id_A <-
    cell_ids <- n_cells <- counts <- NULL

  # --- Validate required columns -------------------------------------------
  .validate_cols(tcr_df, c("chain", "clonotype_id"), "tcr_df")
  if (!is.null(sample_id_columns))
    .validate_cols(tcr_df, sample_id_columns, "tcr_df (sample_id_columns)")

  # --- Fill missing optional columns with NA, with a warning ---------------
  optional_cols <- c("n_cells", "counts",
                     "v_gene", "d_gene", "j_gene",
                     "cdr3_nt", "cdr3_aa", "cdr3_aa_A",
                     "junction", "junction_aa",
                     "contig_id", "contig_id_A", "cell_ids")
  missing_optional <- setdiff(optional_cols, colnames(tcr_df))
  if (length(missing_optional) > 0) {
    warning("get_clonotype_table: optional columns not found, ",
            "filled with NA: ",
            paste(missing_optional, collapse = ", "))
    for (col in missing_optional) tcr_df[[col]] <- NA
  }

  # --- Summarise (via data.table for speed and statistics) -----------------
  group_cols <- c("chain", "clonotype_id", sample_id_columns)
  dt <- data.table::as.data.table(tcr_df)

  df.out <- dt[, {
    .v         <- unique(v_gene[!is.na(v_gene)])
    .d         <- unique(d_gene[!is.na(d_gene)])
    .j         <- unique(j_gene[!is.na(j_gene)])
    .cdr3_nt   <- unique(cdr3_nt[!is.na(cdr3_nt)])
    .cdr3_aa   <- unique(cdr3_aa[!is.na(cdr3_aa)])
    .cdr3_aa_A <- unique(cdr3_aa_A[!is.na(cdr3_aa_A)])
    .junc      <- unique(junction[!is.na(junction)])
    .junc_aa   <- unique(junction_aa[!is.na(junction_aa)])
    .ctg       <- unique(contig_id[!is.na(contig_id)])
    .ctg_A     <- unique(contig_id_A[!is.na(contig_id_A)])
    .cells     <- unique(cell_ids[!is.na(cell_ids)])
    list(
      n_cells       = sum(n_cells, na.rm = TRUE),
      counts        = sum(counts,  na.rm = TRUE),
      v_gene_n      = length(.v),
      v_gene        = paste(.v,         collapse = ","),
      d_gene_n      = length(.d),
      d_gene        = paste(.d,         collapse = ","),
      j_gene_n      = length(.j),
      j_gene        = paste(.j,         collapse = ","),
      cdr3_nt_n     = length(.cdr3_nt),
      cdr3_nt       = paste(.cdr3_nt,   collapse = ","),
      cdr3_aa_n     = length(.cdr3_aa),
      cdr3_aa       = paste(.cdr3_aa,   collapse = ","),
      cdr3_aa_A_n   = length(.cdr3_aa_A),
      cdr3_aa_A     = paste(.cdr3_aa_A, collapse = ","),
      junction_n    = length(.junc),
      junction      = paste(.junc,      collapse = ","),
      junction_aa_n = length(.junc_aa),
      junction_aa   = paste(.junc_aa,   collapse = ","),
      contig_id_n   = length(.ctg),
      contig_id     = paste(.ctg,       collapse = ","),
      contig_id_A_n = length(.ctg_A),
      contig_id_A   = paste(.ctg_A,     collapse = ","),
      cell_id_n     = length(.cells),
      cell_id       = paste(.cells,     collapse = ",")
    )
  }, by = group_cols]

  return(as.data.frame(df.out))
}
