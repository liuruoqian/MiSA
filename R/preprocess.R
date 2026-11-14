
#' Kernel density estimate of cell locations
#'
#' Estimates the spatial density of segmented cells in an image, used to
#' identify the microbial colonization region before building the observation window.
#'
#' @param image A data frame of segmented cells with columns \code{x}, \code{y} and \code{taxon_name}.
#'
#' @return An object with class \code{"density"} whose underlying structure is a list.
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_image)
#' img_kde <- kde_est(example_image)}
kde_est <- function(image){
  max_coord <- max(image$x, image$y)
  img_pp <- spatstat.geom::ppp(x=image$x, y=image$y, window=spatstat.geom::square(max_coord), marks = image$sciname)

  # bandwidth
  bw2 <- spatstat.explore::bw.diggle(img_pp, warn=FALSE)

  # kernel density estimation
  image_kde <- stats::density(img_pp, sigma = bw2, edge = TRUE)

  return (image_kde)
}


#' Build observation windows from alpha shapes of high-density regions
#'
#' Builds closed observation windows from alpha shapes of high-density cell
#' regions. Separated regions are kept apart, and each is returned as its own
#' \code{owin} object, sorted by area from largest to smallest.
#'
#' @param image_kde Output of \code{\link{kde_est}}.
#' @param alpha Positive number controlling how tightly the window follows
#'   the shape of the high-density cell region, in the same units as the coordinates.
#'   Smaller values give a tighter, more detailed boundary and may split the
#'   region into several pieces; larger values give a smoother boundary that
#'   bridges small gaps. As a rough guide, empty gaps wider than about
#'   \code{2 * alpha} separate regions or form holes. Default is 200.
#' @param cutoff Density threshold used to remove low-density areas. Default is 0.002.
#' @param fill_holes Logical. If \code{TRUE} (default), holes within each region are
#'   filled, so each window is bounded by its outer boundary only.
#'   If \code{FALSE}, holes are kept and excluded from the window.
#'
#' @return A named list of \code{owin} objects (\code{region1}, \code{region2}, ...).
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_image)
#' wins <- alpha_windows(kde_est(example_image))
#' wins <- alpha_windows(kde_est(example_image), alpha = 100, cutoff = 2e-3, fill_holes = FALSE)
#' plot(wins$region1)}
alpha_windows <- function(image_kde, alpha = 200, cutoff = 2e-3, fill_holes = TRUE) {

  image_kde.filtered <- image_kde$v
  image_kde.filtered[which(image_kde$v < cutoff)] <- 0

  idx <- which(image_kde.filtered != 0, arr.ind = TRUE)

  # idx is a matrix with row and col indices of nonzero entries
  coords <- data.frame(
    x.orig = image_kde$xcol[idx[, "col"]],
    y.orig = image_kde$yrow[idx[, "row"]]
  )

  x <- coords$x.orig
  y <- coords$y.orig

  tris <- deldir::triang.list(deldir::deldir(x, y, suppressMsge = TRUE))

  keep <- vapply(tris, function(t) {
    a <- sqrt((t$x[1] - t$x[2])^2 + (t$y[1] - t$y[2])^2)
    b <- sqrt((t$x[2] - t$x[3])^2 + (t$y[2] - t$y[3])^2)
    c <- sqrt((t$x[3] - t$x[1])^2 + (t$y[3] - t$y[1])^2)
    tri_area <- abs((t$x[2] - t$x[1]) * (t$y[3] - t$y[1]) -
                  (t$x[3] - t$x[1]) * (t$y[2] - t$y[1])) / 2
    tri_area > 0 && (a * b * c) / (4 * tri_area) <= alpha     # circumradius <= alpha
  }, logical(1))
  if (!any(keep)) stop("No triangles kept: alpha is too small.")

  polys <- lapply(tris[keep], function(t)
    sf::st_polygon(list(cbind(c(t$x, t$x[1]), c(t$y, t$y[1])))))

  # this step we calculate eps, how much the window enlarged in every direction
  # we set it no larger than 2e-5 * the minimum dist between two grid
  s <- image_kde$xstep
  eps <- s * 2e-5

  shape  <- sf::st_buffer(sf::st_make_valid(sf::st_union(sf::st_sfc(polys))), dist = eps)
  pieces <- sf::st_cast(shape, "POLYGON")                 # one polygon per region

  if (fill_holes) {
    pieces <- sf::st_sfc(lapply(pieces, function(p){sf::st_polygon(p[1])}),
    crs = sf::st_crs(pieces))
  }

  pieces <- pieces[order(-as.numeric(sf::st_area(pieces)))]

  wins <- lapply(seq_along(pieces), function(i) spatstat.geom::as.owin(pieces[i]))
  names(wins) <- paste0("region", seq_along(wins))

  print(data.frame(region = names(wins),
                   region_area   = round(vapply(wins, spatstat.geom::area, numeric(1)))))
  return (wins)
}


#' Plot Observation Window Regions
#'
#' Draws the borders of all regions returned by \code{\link{alpha_windows}}
#' on a single plot, with each region labeled by its name.
#'
#' @param wins A named list of \code{owin} objects, as returned by
#'   \code{\link{alpha_windows}}.
#' @param main Plot title. Default is \code{"Alpha-shape regions"}.
#'
#' @return Invisibly returns the union of all regions as a single \code{owin}.
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_image)
#' wins <- alpha_windows(kde_est(example_image))
#' plot_regions(wins)}
plot_regions <- function(wins, main = "Alpha-shape regions") {
  W_all <- do.call(spatstat.geom::union.owin, unname(wins))
  plot(W_all, main = main)
  for (i in seq_along(wins)) {
    p <- sf::st_coordinates(sf::st_point_on_surface(sf::st_as_sfc(wins[[i]])))
    graphics::text(p[1, "X"], p[1, "Y"], labels = names(wins)[i])
  }
  invisible(W_all)
}


