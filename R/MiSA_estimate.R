#####  Main Functions for MiSA  ######

#' Plot cell centroids for a single image, colored by taxon
#'
#' Creates a scatter plot of microbial cell centroids for a single image,
#' with points colored by taxon identity.
#'
#' @param image.df A data frame representing a single image, containing the
#'   columns \code{x}, \code{y}, and \code{taxon_name}, corresponding to the
#'   x and y spatial coordinates and the taxon identity of each cell,
#'   respectively.
#' @param image.title (optional) Character string giving the plot title. Defaults to
#'   \code{NULL}, in which case the title \code{"Scatter Plot by Taxon Name"}
#'   is used.
#'
#' @return A \code{ggplot} object showing the spatial distribution of cell
#'   centroids, colored by taxon.
#'
#' @examples
#' \dontrun{
#' plot_centroid(healthImg.list[['image1']])
#' plot_centroid(diseaseImg.list[['image1']], image.title = "Disease Image 1")
#' }
#'
#' @export
plot_centroid <- function(image.df, image.title=NULL){

  # Validate data frame columns
  required_cols <- c("x", "y", "taxon_name")
  if (!all(required_cols %in% colnames(image.df))) {
    stop("Each image data frame must contain the columns: x, y, taxon_name")
  }

  if (is.null(image.title)){
    image.title <- "Scatter Plot by Taxon Name"
  }
  p1 <- ggplot2::ggplot(image.df, ggplot2::aes(x = x, y = y, color = taxon_name)) +
    ggplot2::geom_point(size = 0.05) +  # Adjust size of points
    ggplot2::labs(title = image.title, x = NULL, y = NULL, color = "Taxon") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8),
      legend.key.size = grid::unit(10, "pt")
    )

  return(p1)
}


#' Estimate L curves for two groups of images
#'
#' Estimates Ripley's \eqn{L}-function for each taxon of interest, across all
#' images in two groups (e.g., health vs. disease), over a pre-specified
#' sequence of radii.
#'
#' @param groupA_image_list A named list of data frames representing the
#'   images in group A. Each data frame must contain the columns \code{x},
#'   \code{y}, and \code{taxon_name}, representing the x and y spatial
#'   coordinates and the corresponding taxon identity, respectively.
#' @param groupB_image_list A named list of data frames representing the
#'   images in group B, with the same structure as \code{groupA_image_list}.
#' @param taxa A character vector giving the taxa of interest for which the
#'   \eqn{L}-function will be estimated.
#' @param r_vec A numeric vector giving the sequence of radii at which the
#'   \eqn{L}-function will be evaluated. This set of radii is used
#'   consistently across all images.
#' @param groupA_window_list (optional) A named list specifying pre-defined
#'   study regions (observation windows) for the images in
#'   \code{groupA_image_list}. Each element must be an object of class
#'   \code{owin}, as required by the \code{window} argument of
#'   \code{spatstat.geom::ppp()}. Names must match the corresponding image
#'   names in \code{groupA_image_list}. A window may be specified for some
#'   images and left as \code{NULL} for others. If \code{NULL} is provided
#'   (either for a single image or the entire list, the default), the
#'   function derives a convex hull from the data as the default window.
#' @param groupB_window_list (optional) A named list specifying pre-defined
#'   study regions for the images in \code{groupB_image_list}, with the same
#'   structure and default behavior as \code{groupA_window_list}.
#' @param border_correct (optional) Character string specifying the edge correction
#'   method used when estimating the \eqn{L}-function. One of \code{"border"}
#'   (the default), \code{"iso"}, or \code{"trans"}. For spatial microbiome
#'   data, \code{"border"} is recommended, as the other two options
#'   substantially increase computation time.
#'
#' @return A list containing the estimated \eqn{L} curves for each
#'   image and taxon.
#'
#' @examples
#' \dontrun{
#' r_all <- seq(0, 1500, by = 6)
#' taxa_all <- unique(healthImg.list[[1]]$taxon_name)
#'
#' L_result <- estimate_L(groupA_image_list=healthImg.list, groupB_image_list=diseaseImg.list,
#' taxa=taxa_all, r_vec=r_all,groupA_window_list = NULL, groupB_window_list = NULL,
#' border_correct = "border")
#' }
#'
#' @export
estimate_L <- function(groupA_image_list, groupB_image_list, taxa, r_vec,
                       groupA_window_list = NULL, groupB_window_list = NULL, border_correct= "border") {

  if (! border_correct %in% c("iso","trans","border")){
    warning(sprintf("Invalid border_correct value ('%s'); using default border correction instead.", border_correct))
    border_correct <- "border"
  }

  if (border_correct %in% c("iso","trans")){
    warning("'border_correct' is set to a non-default value, which may increase computation time. ",
            "It is recommended to use the default: border_correct = 'border'.")
  }

  if (!is.list(groupA_image_list) || is.null(names(groupA_image_list))) {
    stop("'groupA_image_list' must be a named list of data frames.")
  }

  if (!is.list(groupB_image_list) || is.null(names(groupB_image_list))) {
    stop("'groupB_image_list' must be a named list of data frames.")
  }

  # Process each image
  L.A <- lapply(names(groupA_image_list), function(img_name) {

    img <- groupA_image_list[[img_name]]

    # Get window: user-defined, or derive from data if not provided
    win <- groupA_window_list[[img_name]]

    if (is.null(win)) {
      cat("Estimating L for Group A,", img_name, "- using default window (convex hull).\n")
    } else {
      cat("Estimating L for Group A,", img_name, "- using user-provided window.\n")
    }

    # Process single image
    estL_single_image(img=img, win=win, taxa=taxa, r_vec=r_vec, border_correct=border_correct)
  })

  names(L.A) <- names(groupA_image_list)

  L.B <- lapply(names(groupB_image_list), function(img_name) {

    img <- groupB_image_list[[img_name]]

    # Get window: user-defined, or derive from data if not provided
    win <- groupB_window_list[[img_name]]

    if (is.null(win)) {
      cat("Estimating L for Group B,", img_name, "- using default window (convex hull).\n")
    } else {
      cat("Estimating L for Group B,", img_name, "- using user-provided window.\n")
    }

    # Process single image
    estL_single_image(img=img, win=win, taxa=taxa, r_vec=r_vec, border_correct=border_correct)
  })

  names(L.B) <- names(groupB_image_list)

  return(list(groupA_L = L.A, groupB_L = L.B))
}


#' Estimate d-curves from L-function estimates for two groups
#'
#' Computes the spatial summary function \eqn{d(r)} for each taxon and image,
#' based on previously estimated \eqn{L}-function curves for two groups
#' (e.g., health vs. disease). The function \eqn{d(r)} quantifies the
#' deviation of each point pattern from complete spatial randomness (CSR).
#'
#' @param groupA_L A list of estimated \eqn{L}-function curves for group A,
#'   as returned by \code{\link{estimate_L}}.
#' @param groupB_L A list of estimated \eqn{L}-function curves for group B,
#'   as returned by \code{\link{estimate_L}}.
#' @param taxa A character vector giving the taxa of interest for which
#'   \eqn{d(r)} will be computed.
#' @param r_vec A numeric vector giving the sequence of radii at which
#'   \eqn{d(r)} will be evaluated. This should match the sequence of radii
#'   used when estimating \code{groupA_L} and \code{groupB_L}.
#'
#' @return A list containing the estimated \eqn{d}-curve and the
#'   corresponding group label (\code{groupA} or \code{groupB}) for each
#'   image, for each taxon.
#'
#' @examples
#' \dontrun{
#' r_all <- seq(0, 1500, by = 6)
#' taxa_all <- unique(healthImg.list[[1]]$taxon_name)
#' L_result <- estimate_L(healthImg.list, diseaseImg.list,
#'                         taxa = taxa_all, r_vec = r_all)
#' d_result <- estimate_d(groupA_L = est_L$groupA_L, groupB_L = est_L$groupB_L,
#'                         taxa = taxa_all, r_vec = r_all)
#' }
#'
#' @export
estimate_d <- function(groupA_L, groupB_L, taxa, r_vec){

  L1.all <- groupA_L
  L2.all <- groupB_L

  ## calculate d curves
  result.list <- list()
  for(taxon in taxa){
    cat("calculating d for", taxon, "\n")
    taxon_border.df <- data.frame()
    group <- c()

    for (i in 1:length(L1.all)){
      L.res <- L1.all[[i]]
      taxon.res <- L.res[[taxon]]
      if (is.null(taxon.res)) next
      taxon_border.df <- rbind(taxon_border.df, taxon.res)
      group <- append(group,"groupA")
    }

    for (i in 1:length(L2.all)){
      L.res <- L2.all[[i]]
      taxon.res <- L.res[[taxon]]
      if (is.null(taxon.res)) next
      taxon_border.df <- rbind(taxon_border.df, taxon.res)
      group <- append(group,"groupB")
    }

    diff_L <- apply(taxon_border.df, 1, function(a){a-r_vec})
    d <- t(diff_L)
    colnames(d) <- r_vec

    d.mat <- data.frame(group=group)
    d.mat$d <- d

    result.list[[taxon]] <- d.mat
  }

  return(result.list)
}

