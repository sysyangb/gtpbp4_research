clean_positive_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[!is.finite(x) | x <= 0] <- NA_real_
  x
}

count_observed <- function(mat, samples) {
  rowSums(!is.na(mat[, samples, drop = FALSE]))
}

fit_limma_two_group <- function(log2_matrix, sample_ids, group_labels,
                                trend = TRUE, robust = TRUE,
                                lfc_threshold = 1) {
  stopifnot(all(sample_ids %in% colnames(log2_matrix)))
  y <- log2_matrix[, sample_ids, drop = FALSE]
  group <- factor(group_labels[sample_ids], levels = c("Free", "Polysome"))
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit <- limma::lmFit(y, design)
  contrast <- limma::makeContrasts(Polysome - Free, levels = design)
  fit <- limma::contrasts.fit(fit, contrast)
  ebayes_fit <- limma::eBayes(fit, trend = trend, robust = robust)
  treat_fit <- limma::treat(
    fit,
    lfc = lfc_threshold,
    trend = trend,
    robust = robust
  )
  list(ebayes = ebayes_fit, treat = treat_fit, design = design, contrast = contrast)
}

extract_limma_results <- function(fits, protein_metadata) {
  efit <- fits$ebayes
  tfit <- fits$treat
  stopifnot(nrow(protein_metadata) == nrow(efit$coefficients))
  out <- protein_metadata
  out$log2FC_Polysome_vs_Free <- as.numeric(efit$coefficients[, 1])
  out$AveExpr <- as.numeric(efit$Amean)
  out$moderated_t <- as.numeric(efit$t[, 1])
  out$P_Value <- as.numeric(efit$p.value[, 1])
  out$BH_FDR <- stats::p.adjust(out$P_Value, method = "BH")
  out$B_statistic <- if (is.null(efit$lods)) NA_real_ else as.numeric(efit$lods[, 1])
  out$TREAT_t <- as.numeric(tfit$t[, 1])
  out$TREAT_P_Value <- as.numeric(tfit$p.value[, 1])
  out$TREAT_BH_FDR <- stats::p.adjust(out$TREAT_P_Value, method = "BH")
  out
}

add_result_calls <- function(x, cfg) {
  lfc <- x$log2FC_Polysome_vs_Free
  x$TREAT_call <- dplyr::case_when(
    x$TREAT_BH_FDR < cfg$fdr_threshold & lfc > cfg$log2fc_threshold ~ "Up",
    x$TREAT_BH_FDR < cfg$fdr_threshold & lfc < -cfg$log2fc_threshold ~ "Down",
    TRUE ~ "Not_significant"
  )
  x$BH_plus_FC_call <- dplyr::case_when(
    x$BH_FDR < cfg$fdr_threshold & lfc > cfg$log2fc_threshold ~ "Up",
    x$BH_FDR < cfg$fdr_threshold & lfc < -cfg$log2fc_threshold ~ "Down",
    TRUE ~ "Not_significant"
  )
  x$legacy_rawP_plus_FC_call <- dplyr::case_when(
    x$P_Value < cfg$raw_p_threshold & lfc > cfg$log2fc_threshold ~ "Up",
    x$P_Value < cfg$raw_p_threshold & lfc < -cfg$log2fc_threshold ~ "Down",
    TRUE ~ "Not_significant"
  )
  x
}

fit_paired_two_pair_sensitivity <- function(log2_matrix, protein_metadata) {
  samples <- c("free_1", "free_2", "polysome_1", "polysome_2")
  y <- log2_matrix[, samples, drop = FALSE]
  keep <- stats::complete.cases(y)
  y <- y[keep, , drop = FALSE]
  md <- protein_metadata[keep, , drop = FALSE]
  group <- factor(c("Free", "Free", "Polysome", "Polysome"), levels = c("Free", "Polysome"))
  pair <- factor(c("P1", "P2", "P1", "P2"))
  design <- stats::model.matrix(~ pair + group)
  fit <- limma::lmFit(y, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  coef_name <- "groupPolysome"
  out <- md
  out$log2FC_Polysome_vs_Free <- as.numeric(fit$coefficients[, coef_name])
  out$moderated_t <- as.numeric(fit$t[, coef_name])
  out$P_Value <- as.numeric(fit$p.value[, coef_name])
  out$BH_FDR <- stats::p.adjust(out$P_Value, method = "BH")
  out
}
