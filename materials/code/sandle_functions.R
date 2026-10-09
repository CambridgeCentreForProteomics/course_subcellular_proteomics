## =============================================================================
## sandle_functions.R
##
## Streamlined helper functions for SANDLE differential localisation analysis,
## condensed from Lisa Breckels' SANDLE scripts (mainFunctions.R,
## differentialLocalisationFunctions.R, sandle_plotting.R, lisa_plotting.R).
##
## Just `source("sandle_functions.R")` - no package installation needed beyond
## MSnbase and a few CRAN packages.
##
## User-facing functions (everything else is internal):
##   split_by_replicate()        - MSnSet  -> list of per-replicate MSnSets
##   sandle_markers()            - marker-based test   -> DL probability
##   sandle_marker_free()        - marker-free test    -> p-value / FDR
##   combine_results()           - join both result tables
##   plot_dl_results()           - DL probability vs -log10(FDR)
##   plot_marker_profiles()      - fitted niche profiles (QC for marker method)
##   plot_protein_profiles()     - profiles of chosen proteins, both conditions
## =============================================================================

suppressPackageStartupMessages({
  library(MSnbase)
  library(robustbase)
  library(glasso)
  library(mvtnorm)
  library(dplyr)
  library(ggplot2)
})

## -----------------------------------------------------------------------------
## Data helpers
## -----------------------------------------------------------------------------

#' Split one MSnSet into a list with one MSnSet per replicate
split_by_replicate <- function(object, rep_col = "Replicate") {
  reps <- sort(unique(pData(object)[[rep_col]]))
  lapply(reps, function(r) object[, pData(object)[[rep_col]] == r])
}

## marker classes in a dataset (everything except "unknown")
.marker_classes <- function(object, fcol) {
  sort(setdiff(unique(as.character(fData(object)[, fcol])), "unknown"))
}

## proteins quantified in every dataset
.common_proteins <- function(objects) {
  Reduce(intersect, lapply(objects, featureNames))
}

## -----------------------------------------------------------------------------
## Method 1: marker-based SANDLE
## Fit a robust multivariate t-distribution to the markers of each niche, in
## each replicate. Every protein then gets a probability of belonging to each
## niche. Differential localisation probability = probability that a protein's
## niche in condition 1 differs from its niche in condition 2.
## -----------------------------------------------------------------------------

## robust multivariate t fit for one niche (marker profiles in rows)
.fit_robust_mvt <- function(x, lambda = 0.01) {
  mcd <- try(covMcd(x), silent = TRUE)
  if (inherits(mcd, "try-error") || anyNA(mcd$center) || anyNA(mcd$cov)) {
    loc <- colMeans(x); S <- cov(x)
  } else {
    loc <- mcd$center;  S <- mcd$cov
  }
  ## fix any zero/negative variances
  bad <- which(diag(S) <= 0)
  for (j in bad) S[j, j] <- var(x[, j])
  ## fall back to graphical-lasso covariance if S is not positive definite
  ev <- eigen(S, only.values = TRUE)$values
  if (any(is.na(ev)) || any(ev < 0)) {
    S <- glasso(cov(scale(x, scale = FALSE)), rho = lambda,
                penalize.diagonal = FALSE)$w
  }
  S <- S + 1e-16 * diag(ncol(x))
  ## estimate degrees of freedom by maximum likelihood
  nll <- function(df) -sum(dmvt(x, delta = loc, sigma = S, df = df, log = TRUE))
  df  <- suppressWarnings(optim(par = ncol(x) + 2, fn = nll, method = "CG")$par)
  list(location = loc, covariance = S, df = df)
}

## fit all niches in one replicate and return protein x niche probabilities
.sandle_one <- function(object, fcol) {
  X       <- exprs(object)
  classes <- .marker_classes(object, fcol)
  mk      <- as.character(fData(object)[, fcol])
  fits    <- lapply(classes, function(k) .fit_robust_mvt(X[mk == k, , drop = FALSE]))
  names(fits) <- classes
  weights <- as.numeric(table(mk[mk != "unknown"])[classes])
  weights <- weights / sum(weights)

  ## posterior niche probabilities (log-sum-exp for numerical stability)
  loglik <- sapply(seq_along(fits), function(k)
    dmvt(X, delta = fits[[k]]$location, sigma = fits[[k]]$covariance,
         df = fits[[k]]$df, log = TRUE) + log(weights[k]))
  loglik <- loglik - apply(loglik, 1, max)
  probs  <- exp(loglik) / rowSums(exp(loglik))
  dimnames(probs) <- list(rownames(X), classes)
  list(probs = probs, fits = fits)
}

#' Marker-based SANDLE differential localisation
#'
#' @param control,treated lists of MSnSets (one per replicate)
#' @param fcol            feature-data column holding the markers
#' @return data.frame: protein, dl_prob, is_marker, niche_control, niche_treated
#'         (plus the per-replicate fits as attribute "fits", for QC plots)
sandle_markers <- function(control, treated, fcol = "markers") {
  objects <- c(control, treated)
  cond    <- rep(1:2, c(length(control), length(treated)))

  res  <- lapply(objects, .sandle_one, fcol = fcol)
  prot <- .common_proteins(objects)

  ## average niche probabilities over the replicates of each condition
  pbar <- lapply(1:2, function(c)
    Reduce(`+`, lapply(res[cond == c], function(r) r$probs[prot, ])) / sum(cond == c))

  ## P(niche differs) = 1 - sum_k P(niche k | ctrl) * P(niche k | trt)
  ## (exact value of the Monte-Carlo estimate in the original SANDLE code)
  dl <- 1 - rowSums(pbar[[1]] * pbar[[2]])

  markers <- unique(unlist(lapply(objects, function(o)
    featureNames(o)[fData(o)[, fcol] != "unknown"])))

  out <- data.frame(protein       = prot,
                    dl_prob       = unname(dl),
                    is_marker     = prot %in% markers,
                    niche_control = colnames(pbar[[1]])[max.col(pbar[[1]])],
                    niche_treated = colnames(pbar[[2]])[max.col(pbar[[2]])])
  out$dl_prob[out$is_marker] <- 0   # markers are assumed not to move
  attr(out, "fits") <- lapply(res, `[[`, "fits")
  out
}

## -----------------------------------------------------------------------------
## Method 2: marker-free test
## For each matched pair of replicates (control rep i, treated rep i) take the
## difference of the two profiles. If the protein does not move, the
## differences should be zero in every fraction. Test this per protein with an
## F-test (lm(diff ~ fraction) vs lm(diff ~ 0)) and correct with BH FDR.
## No markers and no model fitting required.
## -----------------------------------------------------------------------------

#' Marker-free (paired) differential localisation test
#'
#' @param control,treated lists of MSnSets; control[[i]] is paired with treated[[i]]
#' @return data.frame: protein, p_value, fdr, neg_log10_fdr, F_stat, RSS_diff
sandle_marker_free <- function(control, treated) {
  stopifnot(length(control) == length(treated))
  prot <- .common_proteins(c(control, treated))

  ## array of differences: proteins x fractions x pairs
  D <- simplify2array(lapply(seq_along(control), function(i)
    exprs(control[[i]])[prot, ] - exprs(treated[[i]])[prot, ]))
  n <- dim(D)[3]          # number of pairs
  d <- dim(D)[2]          # number of fractions

  ## vectorised version of anova(lm(diff ~ fraction), lm(diff ~ 0)) per protein
  RSS0  <- apply(D^2, 1, sum)                             # diff ~ 0
  M     <- apply(D, c(1, 2), mean)                        # fraction means
  RSS1  <- RSS0 - n * rowSums(M^2)                        # diff ~ fraction
  df1   <- d; df2 <- n * d - d
  Fstat <- ((RSS0 - RSS1) / df1) / (RSS1 / df2)
  p     <- pf(Fstat, df1, df2, lower.tail = FALSE)

  data.frame(protein       = prot,
             p_value       = p,
             fdr           = p.adjust(p, method = "BH"),
             neg_log10_fdr = -log10(p.adjust(p, method = "BH")),
             F_stat        = Fstat,
             RSS_diff      = RSS0 - RSS1,
             row.names     = NULL)
}

## -----------------------------------------------------------------------------
## Combining and plotting
## -----------------------------------------------------------------------------

#' Join marker-based and marker-free results by protein
combine_results <- function(markers_res, marker_free_res) {
  inner_join(marker_free_res, markers_res, by = "protein")
}

#' Scatter plot of DL probability (marker method) vs -log10 FDR (marker-free)
plot_dl_results <- function(res, prob_cut = 0.9, fdr_cut = 0.1, highlight = NULL) {
  res$hit <- res$dl_prob > prob_cut & res$fdr < fdr_cut
  g <- ggplot(res, aes(dl_prob, neg_log10_fdr)) +
    geom_point(aes(colour = hit, size = RSS_diff), alpha = 0.6) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
    geom_vline(xintercept = prob_cut, linetype = "dashed") +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#2171B5"),
                        labels = c("no", "yes"), name = "Both criteria") +
    scale_size_continuous(range = c(0.3, 4), guide = "none") +
    xlim(0, 1) + theme_bw(base_size = 12) +
    labs(x = "Differential localisation probability (marker method)",
         y = expression(-log[10]~FDR~"(marker-free method)"))
  if (!is.null(highlight)) {
    g <- g + geom_point(data = filter(res, protein %in% highlight),
                        shape = 21, size = 3, colour = "red", stroke = 1)
  }
  g
}

#' Fitted niche profiles (mean +/- 1 sd) for one replicate - QC of marker fits
plot_marker_profiles <- function(markers_res, replicate = 1) {
  fits <- attr(markers_res, "fits")[[replicate]]
  df <- do.call(rbind, lapply(names(fits), function(k) {
    f <- fits[[k]]; s <- sqrt(diag(f$covariance))
    data.frame(niche = k, fraction = seq_along(f$location),
               mean = f$location, lo = f$location - s, hi = f$location + s)
  }))
  ggplot(df, aes(fraction, mean, group = niche, colour = niche)) +
    geom_line() + geom_point() + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
    facet_wrap(~ niche, scales = "free_y") + theme_classic() +
    scale_x_continuous(breaks = seq_len(max(df$fraction))) +
    theme(legend.position = "none") +
    labs(x = "Fraction", y = "Intensity")
}

#' Profiles of selected proteins in both conditions (thin = replicates,
#' thick = mean, ribbon = +/- 1 sd)
plot_protein_profiles <- function(control, treated, proteins,
                                  labels = c("control", "treated")) {
  long <- function(objs, lab) do.call(rbind, lapply(seq_along(objs), function(i) {
    x <- exprs(objs[[i]])[intersect(proteins, featureNames(objs[[i]])), , drop = FALSE]
    data.frame(protein   = rep(rownames(x), ncol(x)),
               fraction  = rep(seq_len(ncol(x)), each = nrow(x)),
               value     = as.vector(x),
               replicate = i, condition = lab)
  }))
  df <- rbind(long(control, labels[1]), long(treated, labels[2]))
  df$condition <- factor(df$condition, labels)
  sm <- df %>% group_by(protein, condition, fraction) %>%
    summarise(m = mean(value), s = sd(value), .groups = "drop")
  ggplot() +
    geom_line(data = df, aes(fraction, value, colour = condition,
                             group = interaction(condition, replicate)), alpha = 0.25) +
    geom_ribbon(data = sm, aes(fraction, ymin = m - s, ymax = m + s, fill = condition),
                alpha = 0.2) +
    geom_line(data = sm, aes(fraction, m, colour = condition), linewidth = 1.1) +
    facet_wrap(~ protein, scales = "free_y") +
    scale_x_continuous(breaks = seq_len(max(df$fraction))) +
    scale_colour_manual(values = c("#636363", "#E6550D")) +
    scale_fill_manual(values = c("#636363", "#E6550D")) +
    theme_minimal() + labs(x = "Fraction", y = "Abundance", colour = NULL, fill = NULL)
}
