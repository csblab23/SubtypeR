#' Boxplot of UCell signature scores with one-vs-rest significance
#'
#' Reproduces the QC boxplot from the original analysis: one box per
#' signature, ordered by median score, annotated with FDR-adjusted
#' Wilcoxon one-vs-rest p-values (each signature's scores vs. all other
#' signatures' scores pooled).
#'
#' @param results A results data.frame as returned by [assign_subtypes()] or
#'   [RunSubtypeR()]`$results`.
#' @param score_suffix Suffix identifying score columns to plot. Default
#'   `"_UCell"`.
#' @param title Plot title. Default `"UCell Signature Scores"`.
#' @param flip Flip to horizontal boxplots (useful with many signatures).
#'   Default `FALSE`.
#'
#' @return A `ggplot` object.
#' @export
plot_subtype_scores <- function(results,
                                 score_suffix = "_UCell",
                                 title = "UCell Signature Scores",
                                 flip = FALSE) {

  score_cols <- grep(paste0(score_suffix, "$"), colnames(results), value = TRUE)
  if (length(score_cols) == 0) {
    stop(sprintf("No columns matching suffix '%s' found in `results`.", score_suffix), call. = FALSE)
  }

  scores_long <- tidyr::pivot_longer(
    results[, c("Cell", score_cols)],
    cols = -("Cell"),
    names_to = "Signature",
    values_to = "Score"
  )

  scores_long <- scores_long %>%
    dplyr::group_by(.data$Signature) %>%
    dplyr::mutate(med = stats::median(.data$Score)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Signature = stats::reorder(.data$Signature, .data$med))

  one_vs_rest <- scores_long %>%
    dplyr::group_by(.data$Signature) %>%
    dplyr::group_modify(function(df, key) {
      test <- stats::wilcox.test(
        scores_long$Score ~ (scores_long$Signature == key$Signature)
      )
      data.frame(p_value = test$p.value)
    }) %>%
    dplyr::ungroup()

  one_vs_rest$p_adj <- stats::p.adjust(one_vs_rest$p_value, method = "fdr")
  one_vs_rest$label <- ifelse(
    one_vs_rest$p_adj < 1e-4,
    "p < 1e-4",
    paste0("p = ", signif(one_vs_rest$p_adj, 3))
  )

  p <- ggplot2::ggplot(scores_long, ggplot2::aes(x = .data$Signature, y = .data$Score)) +
    ggplot2::geom_boxplot(outlier.size = 0.3, color = "black", fill = "white") +
    ggplot2::geom_text(
      data = one_vs_rest,
      ggplot2::aes(x = .data$Signature, y = max(scores_long$Score) * 1.05, label = .data$label),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 14, face = "bold", color = "black"),
      axis.text.y = ggplot2::element_text(size = 14, face = "bold", color = "black"),
      axis.title = ggplot2::element_text(size = 18, face = "bold", color = "black"),
      plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5, color = "black")
    ) +
    ggplot2::labs(title = title, x = "Signature", y = "UCell Score")

  if (isTRUE(flip)) p <- p + ggplot2::coord_flip()

  p
}

#' Bar plot of subtype composition
#'
#' @param results A results data.frame containing a `Subtype` column, as
#'   returned by [assign_subtypes()] or [RunSubtypeR()]`$results`.
#' @param title Plot title. Default `"Subtype Composition"`.
#'
#' @return A `ggplot` object.
#' @export
plot_subtype_composition <- function(results, title = "Subtype Composition") {
  if (!"Subtype" %in% colnames(results)) {
    stop("`results` must contain a `Subtype` column.", call. = FALSE)
  }

  tab <- as.data.frame(table(Subtype = results$Subtype))

  ggplot2::ggplot(tab, ggplot2::aes(x = stats::reorder(.data$Subtype, -.data$Freq), y = .data$Freq)) +
    ggplot2::geom_col(fill = "grey30") +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 13, face = "bold", color = "black"),
      axis.text.y = ggplot2::element_text(size = 13, face = "bold", color = "black"),
      axis.title = ggplot2::element_text(size = 16, face = "bold", color = "black"),
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black")
    ) +
    ggplot2::labs(title = title, x = "Subtype", y = "Number of cells")
}
