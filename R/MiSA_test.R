
### function to filter/truncate d ###

#' Prepare d-curves for hypothesis testing
#'
#' Prepares \eqn{d}-curves for downstream hypothesis testing by filtering
#' out images with limited spatial range and truncating the remaining
#' curves to a common maximum radius.
#'
#'
#' @param est_d A list of estimated \eqn{d}-curves and group labels for each
#'   taxon and image, as returned by \code{\link{estimate_d}}.
#' @param r_thresh (optional) Numeric value giving the radius threshold used to filter
#'   out images with limited spatial range. Images whose \eqn{d}-curve does
#'   not have valid values up to \code{r_thresh} are excluded from further
#'   analysis. Defaults to \code{NULL}, in which case images whose
#'   \eqn{d}-curve does not reach 0.5 times the maximum radius considered while
#'   estimating \eqn{d}-curve. Users may set \code{r_thresh} based on the
#'   scientific question of interest; for example, \code{r_thresh = 500}
#'   will exclude \eqn{d}-curves that do not have valid estimates up to
#'   \eqn{r = 500}.
#'
#' @return A list of filtered and truncated \eqn{d}-curves, ready for input
#'   into downstream hypothesis testing functions.
#'
#' @examples
#' \dontrun{
#' data_pre <- MiSA_data(est_d)
#' data_pre <- MiSA_data(est_d, r_thresh = 500)
#' }
#'
#' @export
MiSA_data <- function(est_d, r_thresh=NULL){

  MiSA_data_list <- lapply(est_d, trim_dmat_eachTaxon, r_thresh=r_thresh)
  return(MiSA_data_list)
}

### drop the taxon that not valid for test ###

#' Check taxa for validity prior to hypothesis testing
#'
#' Checks each taxon against a set of minimum requirements
#' for hypothesis testing, and drops any taxon that does not meet these
#' requirements.
#'
#' A taxon is dropped if any of the following conditions are not met:
#' \enumerate{
#'   \item At least two groups remain after filtering.
#'   \item The total sample size is greater than 2.
#'   \item The corresponding \code{d} has more than one valid estimate.
#' }
#'
#'
#' @param data_list A named list, with one element per taxon, as returned by
#'   \code{\link{MiSA_data}}.
#'
#' @return A named list with the same structure as \code{data}, containing
#'   only the taxa that passed all validity checks.
#'
#'@examples
#' \dontrun{
#' data_checked <- MiSA_check(data_list = data_pre)
#' }
#'
#' @export
MiSA_check <- function(data_list){

  for (taxon in names(data_list)){
    each_data <- data_list[[taxon]]

    group.filtered <- each_data$group
    d <- each_data$d

    cond1 <- length(unique(group.filtered)) == 2
    cond2 <- length(group.filtered) >= 3
    cond3 <- ncol(d) > 1

    passes_check <- cond1 && cond2 && cond3

    if (!passes_check) {
      reasons <- c(
        if (!cond1) "only one group remains after filtering, not valid for hypothesis testing",
        if (!cond2) "total sample size is not greater than 2, not sufficient for testing",
        if (!cond3) "d has one or fewer valid estimate"
      )
      message(sprintf("Dropped '%s': %s", taxon, paste(reasons, collapse = "; ")))

      data_list[[taxon]] <- NULL
    }
  }

  return (data_list)

}

### perform test using MiSA ###


#' Test for differential spatial patterns using MiSA
#'
#' Performs hypothesis testing to detect differential spatial patterns
#' between two groups, for each taxon in \code{MiSA_data}, using one of
#' four available MiSA tests.
#'
#' Four tests are available via the \code{method} argument:
#' \describe{
#'   \item{\code{"orig"}}{Evaluates the raw statistic \eqn{d(r)} (capturing
#'     both the magnitude and direction of deviation from CSR) at a set of
#'     pre-selected radii, and combines the resulting p-values across radii
#'     using the harmonic mean p-value (HMP) method.}
#'   \item{\code{"abs"}}{Similar to \code{"orig"}, but compares the absolute
#'     value of \eqn{d(r)} (capturing only the magnitude of deviation) at the
#'     selected radii instead.}
#'   \item{\code{"sup"}}{Evaluates the maximum absolute deviation,
#'     \eqn{\max_r |d(r)|}, capturing signals that may arise at any spatial
#'     scale.}
#'   \item{\code{"area"}}{Evaluates the integrated magnitude of deviation of
#'     \eqn{d(r)} across the full range of radii, capturing signals that may
#'     arise at any spatial scale.}
#' }
#'
#' For \code{"orig"} and \code{"abs"}, the radii used can be specified via
#' \code{omnibus_r} or \code{omnibus_num}. If \code{omnibus_r} is specified,
#' the function always uses \code{omnibus_r}, regardless of the input for
#' \code{omnibus_num}. If \code{omnibus_r} is left at its default
#' (\code{NULL}) and \code{omnibus_num} is specified, the function uses
#' \code{omnibus_num} evenly spaced radii instead. If both \code{omnibus_r}
#' and \code{omnibus_num} are left at their defaults (\code{NULL}), the
#' function sets \code{omnibus_num = min(5, c - 1)}, where \code{c} is the
#' number of radii available for that taxon.
#'
#' This function uses whichever group is entered as \code{groupA} (during
#' data preparation) as the reference group. A negative effect size
#' indicates that the summary statistic is smaller in the reference group
#' than in the comparison group, indicating that the comparison group is
#' more spatially structured (e.g., exhibiting stronger clustering or
#' regularity) relative to the reference group.
#'
#' @param misa_data_list A named list, with one element per taxon, as
#'   returned by \code{\link{MiSA_data}} and \code{\link{MiSA_check}}. Each
#'   element must contain \code{group}, a vector of group labels, and
#'   \code{d}, a matrix of estimated \eqn{d}-curves.
#' @param method Character string specifying which test to perform. One of
#'   \code{"orig"}, \code{"abs"}, \code{"sup"}, or \code{"area"}. Defaults to
#'   \code{"area"}.
#' @param omnibus_r Optional numeric vector specifying the specific radii at
#'   which to evaluate the test statistic, used only when
#'   \code{method = "orig"} or \code{method = "abs"}. Defaults to
#'   \code{NULL}.
#' @param omnibus_num Optional integer specifying the number of radii,
#'   evenly spaced across the available range, at which to evaluate the
#'   test statistic, used only when \code{method = "orig"} or
#'   \code{method = "abs"} and \code{omnibus_r} is not specified. Defaults
#'   to \code{NULL}.
#'
#' @return A list containing taxon-specific p-values and a community-level
#'   p-value assessing overall spatial dysbiosis. For \code{method = "area"}
#'   or \code{method = "sup"}, taxon-specific effect sizes are additionally
#'   returned. For \code{method = "orig"} or \code{method = "abs"}, a
#'   \code{per_r_pval} element is additionally returned, describing the
#'   testing results at each of the specified radii.
#'
#' @examples
#' \dontrun{
#' result <- MiSA(misa_data_list = data_checked, method = "area")
#' result <- MiSA(misa_data_list = data_checked, method = "orig", omnibus_r = c(300, 600, 900))
#' result <- MiSA(misa_data_list = data_checked, method = "abs", omnibus_num = 6)
#' }
#'
#' @export
MiSA <- function(misa_data_list, method="area", omnibus_r=NULL, omnibus_num = NULL){

  if (length(method) > 1) {
    stop("Please select only one method from: 'area', 'sup', 'orig', 'abs'.")
  }

  if (! method %in% c("area", "sup","orig","abs")){
    stop(sprintf("Invalid method '%s'. Please select one of: 'area', 'sup', 'orig', 'abs'.", method))
  }

  if (method == "area"){

    res_area <- lapply(misa_data_list, function(data) {
      result <- MiSA_area_eachTaxon(data$d, data$group)
      as.data.frame(result)
    })

    results_area.df <- do.call(rbind, res_area)

    p.v <- results_area.df$p.value
    p.comb <- harmonicmeanp::p.hmp(p.v, L=length(p.v))

    results <- list(taxon_results = results_area.df, community_pval = p.comb)
  }

  if (method=="sup"){

    res_sup <- lapply(misa_data_list, function(data) {
      result <- MiSA_sup_eachTaxon(data$d, data$group)
      as.data.frame(result)
    })

    results_sup.df <- do.call(rbind, res_sup)

    p.v <- results_sup.df$p.value
    p.comb <- harmonicmeanp::p.hmp(p.v, L=length(p.v))

    results <- list(taxon_results = results_sup.df, community_pval = p.comb)
  }

  if (method == "orig"){

    p_orig <- lapply(seq_along(misa_data_list), function(i) {
      data <- misa_data_list[[i]]
      res <- MiSA_orig_eachTaxon(data$d, data$group, omnibus_r = omnibus_r, omnibus_num = omnibus_num)
      return(res)
    })

    names(p_orig) <- names(misa_data_list)

    p.v <- sapply(p_orig, function(x) x[[1]])   # combined p for each taxon
    p.comb <- harmonicmeanp::p.hmp(p.v, L = length(p.v)) # community level p

    per_r_list <- lapply(p_orig, function(x) x[-1])
    names(per_r_list) <- names(p_orig)

    results <- list(
      taxon_results = data.frame(p.value = p.v),
      community_pval = p.comb,
      per_r_pval = per_r_list
    )
  }

  if (method == "abs"){

    p_abs <- lapply(seq_along(misa_data_list), function(i) {
      data <- misa_data_list[[i]]
      res <- MiSA_abs_eachTaxon(data$d, data$group, omnibus_r=omnibus_r, omnibus_num = omnibus_num)
      return(res)
    })

    names(p_abs) <- names(misa_data_list)

    p.v <- sapply(p_abs, function(x) x[[1]])   # combined p for each taxon
    p.comb <- harmonicmeanp::p.hmp(p.v, L = length(p.v)) # community level p

    per_r_list <- lapply(p_abs, function(x) x[-1])
    names(per_r_list) <- names(p_abs)

    results <- list(
      taxon_results = data.frame(p.value = p.v),
      community_pval = p.comb,
      per_r_pval = per_r_list
    )

  }

  return(results)
}



