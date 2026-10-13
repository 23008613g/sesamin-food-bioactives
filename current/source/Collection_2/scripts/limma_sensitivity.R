#!/usr/bin/env Rscript
# Genome-wide limma sensitivity analysis; no candidate-based selection.
# Usage: Rscript limma_sensitivity.R D:/ResearchData/pfhn_genome39_2026-10-13
# All three input files must be generated under the frozen protocol first.
# This file does not transform, normalize, filter, or select input genes.

options(stringsAsFactors = FALSE, digits = 17, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply exactly one absolute stage directory.")
stage <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
out_dir <- file.path(stage, "results", "limma")
review_file <- file.path(stage, "review", "limma_execution.json")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(review_file), recursive = TRUE, showWarnings = FALSE)

# Small base-R JSON encoder: no dependency installation is needed.
json_quote <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub('"', '\\"', x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x <- gsub("\r", "\\r", x, fixed = TRUE)
  x <- gsub("\t", "\\t", x, fixed = TRUE)
  paste0('"', x, '"')
}
to_json <- function(x) {
  if (is.null(x)) return("null")
  if (is.list(x)) {
    vals <- vapply(x, to_json, character(1))
    if (!is.null(names(x))) {
      return(paste0("{", paste(paste0(json_quote(names(x)), ":", vals), collapse = ","), "}"))
    }
    return(paste0("[", paste(vals, collapse = ","), "]"))
  }
  if (length(x) != 1L) return(to_json(as.list(unname(x))))
  if (is.na(x)) return("null")
  if (is.logical(x)) return(if (x) "true" else "false")
  if (is.numeric(x)) {
    if (!is.finite(x)) return(json_quote(as.character(x)))
    return(sprintf("%.17g", x))
  }
  json_quote(as.character(x))
}
write_json <- function(x, path) writeLines(to_json(x), path, useBytes = TRUE)
utc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
execution <- list(status = "started", started_utc = utc_now(), stage = stage,
                  biological_input_read = FALSE,
                  warnings = list(),
                  analysis_role = "Genome-wide exploratory sensitivity analysis under parameters frozen before this computation; not independent biological replication.")
write_json(execution, review_file)

main <- function() {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Required existing limma installation is unavailable.")
  limma_version <- as.character(utils::packageVersion("limma"))
  if (limma_version != "3.68.4") stop("Unexpected limma version: ", limma_version)
  if (as.character(getRversion()) != "4.6.1") stop("Unexpected R version: ", getRversion())

  samples <- paste0("GSM", 1048036:1048045)
  groups <- factor(rep(c("control", "sesamin"), each = 5L), levels = c("control", "sesamin"))
  names(groups) <- samples
  design <- stats::model.matrix(~ groups)
  rownames(design) <- samples
  contrast <- matrix(c(0, 1), ncol = 1L,
                     dimnames = list(colnames(design), "sesamin_minus_control"))
  stopifnot(nrow(design) == 10L, qr(design)$rank == 2L,
            identical(rownames(design), samples), all(table(groups) == 5L))
  utils::write.csv(data.frame(sample = samples, group = as.character(groups), design,
                             check.names = FALSE), file.path(out_dir, "design.csv"), row.names = FALSE)
  utils::write.csv(data.frame(coefficient = rownames(contrast), contrast, check.names = FALSE),
                   file.path(out_dir, "contrast.csv"), row.names = FALSE)
  params <- list(model = "~ groups; control is the reference level",
                 contrast = "sesamin minus control; positive logFC means higher sesamin-group signal",
                 trend = TRUE, robust = FALSE, proportion = 0.01,
                 stdev.coef.lim = c(0.1, 4), confidence_level = 0.95,
                 adjustments = c("BH separately over every supplied gene for each variant",
                                  "BY separately over every supplied gene for each variant"),
                 transformations_in_this_script = "none",
                 filtering_in_this_script = "none",
                 interval = "Marginal 95% moderated t interval; not simultaneous or multiplicity-adjusted",
                 SE_prior_definition = "stdev.unscaled * sqrt(s2.prior); trend-dependent prior residual variance",
                 SE_moderated_definition = "stdev.unscaled * sqrt(s2.post); used for moderated t and CI",
                 inputs = "Frozen, aligned gene-level log2-signal matrices; no further log transformation")
  write_json(params, file.path(out_dir, "parameters.json"))

  assert_input <- function(dat, label) {
    expected <- c("rat_entrez", samples)
    if (!identical(names(dat), expected)) stop(label, ": unexpected columns or sample order.")
    ids <- as.character(dat$rat_entrez)
    if (nrow(dat) < 2L || anyNA(ids) || any(ids == "") || anyDuplicated(ids)) {
      stop(label, ": empty or nonunique gene identifiers.")
    }
    if (any(!grepl("^[0-9]+$", ids))) stop(label, ": rat_entrez must contain integer Entrez identifiers.")
    if (!all(vapply(dat[samples], is.numeric, logical(1)))) stop(label, ": expression columns must be numeric.")
    x <- as.matrix(dat[samples])
    storage.mode(x) <- "double"
    rownames(x) <- ids
    if (any(!is.finite(x))) stop(label, ": nonfinite values; no silent deletion is permitted.")
    stopifnot(identical(colnames(x), rownames(design)), ncol(x) == 10L)
    x
  }
  fit_matrix <- function(x) {
    base_fit <- limma::lmFit(x, design = design)
    cf <- limma::contrasts.fit(base_fit, contrasts = contrast)
    fit <- limma::eBayes(cf, trend = TRUE, robust = FALSE, proportion = 0.01,
                         stdev.coef.lim = c(0.1, 4))
    manual_diff <- rowMeans(x[, groups == "sesamin", drop = FALSE]) -
      rowMeans(x[, groups == "control", drop = FALSE])
    delta <- max(abs(as.numeric(fit$coefficients[, 1L]) - manual_diff))
    if (!is.finite(delta) || delta > 1e-10) stop("Coefficient differs from arithmetic mean contrast.")
    n <- nrow(x)
    repeat_to_n <- function(v) {
      if (!(length(v) %in% c(1L, n))) stop("Unexpected variance/df vector dimension.")
      rep_len(as.numeric(v), n)
    }
    unscaled <- as.numeric(fit$stdev.unscaled[, 1L])
    se <- unscaled * sqrt(repeat_to_n(fit$s2.post))
    df <- repeat_to_n(fit$df.total)
    q <- stats::qt(0.975, df = df)
    logfc <- as.numeric(fit$coefficients[, 1L])
    pv <- as.numeric(fit$p.value[, 1L])
    tstat <- as.numeric(fit$t[, 1L])
    if (any(!is.finite(se)) || any(se <= 0) || anyNA(pv) || any(pv < 0 | pv > 1)) {
      stop("Unexpected moderated standard errors or p values.")
    }
    t_difference <- max(abs(tstat - logfc / se))
    p_difference <- max(abs(pv - 2 * stats::pt(-abs(tstat), df = df)))
    if (t_difference > 1e-10 || p_difference > 1e-10) stop("Moderated t/p-value identity check failed.")
    result <- data.frame(
      gene_id = rownames(x),
      logFC = logfc, AveExpr = as.numeric(fit$Amean),
      mean_control = rowMeans(x[, groups == "control", drop = FALSE]),
      mean_sesamin = rowMeans(x[, groups == "sesamin", drop = FALSE]),
      sd_control = apply(x[, groups == "control", drop = FALSE], 1L, stats::sd),
      sd_sesamin = apply(x[, groups == "sesamin", drop = FALSE], 1L, stats::sd),
      n_control = 5L, n_sesamin = 5L,
      t = tstat, P.Value = pv,
      adj.P.Val = stats::p.adjust(pv, method = "BH"),
      adj.P.Val.BY = stats::p.adjust(pv, method = "BY"),
      B = as.numeric(fit$lods[, 1L]),
      stdev.unscaled = unscaled,
      sigma = repeat_to_n(fit$sigma),
      s2.prior = repeat_to_n(fit$s2.prior),
      s2.post = repeat_to_n(fit$s2.post),
      SE_unmoderated = unscaled * repeat_to_n(fit$sigma),
      SE_prior = unscaled * sqrt(repeat_to_n(fit$s2.prior)),
      SE_moderated = se,
      df.residual = repeat_to_n(fit$df.residual),
      df.prior = repeat_to_n(fit$df.prior), df.total = df,
      CI95_moderated_low = logfc - q * se,
      CI95_moderated_high = logfc + q * se,
      check.names = FALSE)
    stopifnot(nrow(result) == nrow(x), identical(result$gene_id, rownames(x)),
              all(result$CI95_moderated_low <= result$logFC),
              all(result$CI95_moderated_high >= result$logFC))
    list(fit = fit, table = result, checks = list(
      max_abs_logFC_minus_mean_difference = delta,
      max_abs_t_identity_difference = t_difference,
      max_abs_p_identity_difference = p_difference,
      all_input_genes_retained_in_original_order = TRUE))
  }

  # Explicitly synthetic deterministic software fixture: not biological evidence.
  # Varied residual patterns permit eBayes/trend estimation without relying on RNG.
  toy_n <- 40L
  toy <- t(vapply(seq_len(toy_n), function(i) {
    base <- 3 + i / 7
    effect <- rep(c(-1.2, 0, 0.8, 1.7), length.out = toy_n)[i]
    noise_control <- c(-2, -1, 0, 1, 2) * (0.04 + (i %% 9) / 70)
    noise_treat <- c(2, -2, 1, -1, 0) * (0.06 + (i %% 7) / 60)
    base + c(noise_control, effect + noise_treat)
  }, numeric(10L)))
  colnames(toy) <- samples
  rownames(toy) <- as.character(9000000L + seq_len(toy_n))
  toy_dat <- data.frame(rat_entrez = rownames(toy), toy, check.names = FALSE)
  toy_x <- assert_input(toy_dat, "synthetic fixture")
  toy_fit <- fit_matrix(toy_x)
  toy_expected <- rep(c(-1.2, 0, 0.8, 1.7), length.out = toy_n)
  stopifnot(max(abs(toy_fit$table$logFC - toy_expected)) < 1e-10,
            nrow(toy_fit$table) == toy_n,
            all(toy_fit$table$logFC[toy_expected > 0] > 0),
            all(toy_fit$table$logFC[toy_expected < 0] < 0))
  bad_order <- toy_dat[, c("rat_entrez", rev(samples))]
  bad_order_rejected <- inherits(try(assert_input(bad_order, "synthetic reversed columns"), silent = TRUE), "try-error")
  stopifnot(bad_order_rejected)
  toy_swapped <- toy_x[, c(6:10, 1:5)]
  colnames(toy_swapped) <- samples
  swapped_fit <- fit_matrix(toy_swapped)
  swap_error <- max(abs(swapped_fit$table$logFC + toy_fit$table$logFC))
  stopifnot(swap_error < 1e-10)
  toy_record <- list(status = "passed", type = "Synthetic deterministic software test only",
                     n_genes = toy_n, all_rows_retained = TRUE,
                     known_mean_differences_and_directions = TRUE,
                     reversed_column_order_rejected = bad_order_rejected,
                     group_swap_sign_reversal_max_error = swap_error,
                     numerical_identity_checks = toy_fit$checks)
  write_json(toy_record, file.path(out_dir, "synthetic_software_test.json"))

  variants <- c("deposited", "qn10", "qn20")
  input_paths <- setNames(file.path(stage, "results", paste0("rat_gene_expression_", variants, ".csv.gz")), variants)
  if (!all(file.exists(input_paths))) stop("Input missing: ", paste(input_paths[!file.exists(input_paths)], collapse = ", "))
  inputs <- lapply(variants, function(v) {
    con <- gzfile(input_paths[[v]], open = "rt")
    on.exit(close(con), add = TRUE)
    dat <- utils::read.csv(con, check.names = FALSE, colClasses = c("character", rep("numeric", 10L)))
    assert_input(dat, v)
  })
  names(inputs) <- variants
  if (!all(vapply(inputs, function(x) identical(rownames(x), rownames(inputs[[1L]])), logical(1)))) {
    stop("Frozen gene universe/order differs between input variants.")
  }
  execution$biological_input_read <<- TRUE
  input_md5 <- as.list(unname(tools::md5sum(input_paths)))
  names(input_md5) <- variants
  summary <- vector("list", length(variants))
  names(summary) <- variants
  for (v in variants) {
    x <- inputs[[v]]
    result <- fit_matrix(x)
    path <- file.path(out_dir, paste0("rat_limma_", v, ".csv.gz"))
    con <- gzfile(path, open = "wt")
    tryCatch(utils::write.csv(result$table, con, row.names = FALSE, na = "NA"),
             finally = close(con))
    saveRDS(result$fit, file.path(out_dir, paste0("rat_limma_", v, "_fit.rds")), compress = "gzip")
    diag_path <- file.path(out_dir, paste0("plotSA_", v, ".pdf"))
    grDevices::pdf(diag_path, width = 7.2, height = 5.2, useDingbats = FALSE)
    tryCatch(limma::plotSA(result$fit, main = paste("Rat liver, sesamin vs control:", v)),
             finally = grDevices::dev.off())
    summary[[v]] <- list(n_genes = nrow(x), n_control = 5L, n_sesamin = 5L,
                         n_BH_005 = sum(result$table$adj.P.Val < 0.05),
                         n_BY_005 = sum(result$table$adj.P.Val.BY < 0.05),
                         n_BH_005_positive = sum(result$table$adj.P.Val < 0.05 & result$table$logFC > 0),
                         n_BH_005_negative = sum(result$table$adj.P.Val < 0.05 & result$table$logFC < 0),
                         p_min = min(result$table$P.Value),
                         BH_min = min(result$table$adj.P.Val),
                         median_df_prior = stats::median(result$table$df.prior),
                         df_total_range = range(result$table$df.total),
                         checks = result$checks,
                         input_path = input_paths[[v]], input_md5 = input_md5[[v]],
                         output_path = path, output_md5 = unname(tools::md5sum(path)))
  }
  capture.output(utils::sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
  execution$status <<- "completed"
  execution$completed_utc <<- utc_now()
  execution$R_version <<- R.version.string
  execution$limma_version <<- limma_version
  execution$sample_order <<- samples
  execution$groups <<- as.character(groups)
  execution$parameters <<- params
  execution$synthetic_software_test <<- toy_record
  execution$variants <<- summary
  execution$limitations <<- c(
    "Empirical-Bayes results are model-based sensitivity analyses, not substitutes for exact label-permutation results.",
    "BH is performed within each preprocessing variant over the full supplied gene universe; significance does not imply robustness across variants.",
    "Individual moderated intervals are not adjusted for simultaneous inference.",
    "Animal selection and the five-versus-five design remain limitations of the deposited study.",
    "No normal human hepatocyte response or clinical effect is established by this analysis.")
  write_json(execution, review_file)
  message("Completed all three limma variants with full-gene BH/BY adjustment: ", review_file)
}

tryCatch(withCallingHandlers(main(), warning = function(w) {
  execution$warnings <<- c(execution$warnings, list(conditionMessage(w)))
}), error = function(e) {
  execution$status <- "failed"
  execution$failed_utc <- utc_now()
  execution$error <- conditionMessage(e)
  write_json(execution, review_file)
  stop(e)
})
