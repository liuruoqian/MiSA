#' Sample image-based spatial microbiome data for the health group
#'
#' A named list of data frames containing example image-based spatial microbiome
#' data for the health group, used to demonstrate the MiSA
#' workflow.
#'
#' @format A named list, with one element per image. Each element is a
#'   data frame with columns:
#' \describe{
#'   \item{x}{x spatial coordinate of the microbial cell}
#'   \item{y}{y spatial coordinate of the microbial cell}
#'   \item{taxon_name}{taxon identity of the microbial cell}
#' }
#' @source A simulated dataset generated using oral biofilm data as a
#'   template.
"healthImg.list"

#' Sample image-based spatial microbiome data for disease group
#'
#' A named list of data frames containing example image-based spatial microbiome
#' data for the disease (e.g., mucositis) group, used to
#' demonstrate the MiSA workflow.
#'
#' @format A named list, with one element per image. Each element is a
#'   data frame with columns:
#' \describe{
#'   \item{x}{x spatial coordinate of the microbial cell}
#'   \item{y}{y spatial coordinate of the microbial cell}
#'   \item{taxon_name}{taxon identity of the microbial cell}
#' }
#' @source A simulated dataset generated using oral biofilm data as a
#'   template.
"diseaseImg.list"


#' Example image for preprocessing
#'
#' An example image used to demonstrate preprocessing step for
#' identifying the region of study in each image.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{x}{x spatial coordinate of the microbial cell}
#'   \item{y}{y spatial coordinate of the microbial cell}
#'   \item{taxon_name}{taxon identity of the microbial cell}
#' }
#' @source Derived from data described in Grodner et al. (2024), \url{https://doi.org/10.1101/2024.07.24.604873}.
#' @examples
#' head(example_image)
"example_image"

