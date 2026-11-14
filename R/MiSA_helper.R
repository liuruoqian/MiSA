
###############  MiSA: helper functions  ###############

### estimate L curves for one single image (all taxa) ###

estL_single_image <- function(img, win, taxa, r_vec, border_correct) {

  # Validate img data columns
  required_cols <- c("x", "y", "taxon_name")
  if (!all(required_cols %in% colnames(img))) {
    stop("Each image data frame must contain the columns: x, y, taxon_name")
  }

  img <- img[,c("x","y","taxon_name")]

  # window is derived from the input data (convex hull of the points)
  if (is.null(win)) {
    win <- spatstat.geom::convexhull.xy(x = img$x, y = img$y)
  }

  L.res <- list()
  for (taxon in taxa){
    pp.df <- subset(img, taxon_name == taxon)

    # filter out low abundant taxa
    if (nrow(pp.df) < 10){
      cat(taxon, "has < 10 cells; excluded from analysis.\n")
      next
    }

    pp <- spatstat.geom::ppp(x=pp.df$x, y=pp.df$y, window=win)

    L <- spatstat.explore::Lest(pp, r=r_vec, correction = border_correct)
    L.res[[taxon]] <- L[[border_correct]]
  }

  return (L.res)
}

###  Prepare Data for Hypothesis Testing  ###

trim_dmat_eachTaxon <- function(d.mat, r_thresh){

  group <- d.mat$group
  d <- d.mat$d
  r <- as.numeric(colnames(d))

  if (is.null(r_thresh)){
    r_thresh <- 0.5 * max(r)
  }

  # exclude curves that don't reach r_thresh
  d.filtered <- d[rowSums(!is.na(d))>=max(which(r <= r_thresh)),]
  group.filtered <- group[rowSums(!is.na(d))>=max(which(r <= r_thresh))]

  # Truncation:
  # For the retained curves, truncate each to the maximum r shared across all of them.
  min_l <- min(rowSums(!is.na(d.filtered)))

  d.truncated <- d.filtered[,1:min_l]
  theo_truncated <- r[1:min_l]
  r_truncated <- r[1:min_l]

  # groupA is used as the reference level; this affects the sign of the effect size
  group.filtered <- factor(group.filtered, levels = c("groupA", "groupB"))
  res <- list(d=d.truncated, group=group.filtered)

  return (res)
}


###  MiSA testing helper functions  ###

test_pure_dist <- function(r, diff_mat, group){
  d <- diff_mat[r, ] # select one row in diff_L mat
  d.new <- d[!is.na(d)] # exclude NAs

  group.new <- group[!is.na(d)] # exclude group values for which d is NA
  wil <- suppressWarnings(stats::wilcox.test(d.new ~ group.new)) # use updated d and group
  return(wil$p.value)
}

test_abs_dist <- function(r, diff_mat, group){
  d <- diff_mat[r, ] # select one row in diff_L mat
  d.new <- d[!is.na(d)] # exclude NAs

  group.new <- group[!is.na(d)] # exclude group values for which d is NA
  wil <- suppressWarnings(stats::wilcox.test(abs(d.new) ~ group.new)) # use updated d and group
  return(wil$p.value)
}

MiSA_area_eachTaxon <- function(d.truncated, group.filtered){

  ## ------------- compute stat for effect size ---------

  n1 <- table(group.filtered)[1]
  n2 <- table(group.filtered)[2]
  mu_W <- n1 * n2 / 2
  sigma_W <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)

  ## ------------ absolute area ------------------

  abs_area <- apply(d.truncated, 1, function(a){sum(abs(a))})
  wilcox_area <- suppressWarnings(stats::wilcox.test(abs_area ~ group.filtered))
  p_area <- wilcox_area$p.value

  W_area <- wilcox_area$statistic
  Z_area <- (W_area - mu_W) / sigma_W
  r_area <- Z_area / sqrt(n1 + n2)

  area_res <- list(p.value = p_area, effect.size = r_area)

  return (area_res)
}

MiSA_sup_eachTaxon <- function(d.truncated, group.filtered){

  ## ------------- compute stat for effect size ---------

  n1 <- table(group.filtered)[1]
  n2 <- table(group.filtered)[2]
  mu_W <- n1 * n2 / 2
  sigma_W <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)

  sup_dist <- apply(d.truncated, 1, function(a){max(abs(a))})
  wilcox_sup <- suppressWarnings(stats::wilcox.test(sup_dist ~ group.filtered))
  p_sup <- wilcox_sup$p.value

  W_sup <- wilcox_sup$statistic
  Z_sup <- (W_sup - mu_W) / sigma_W
  r_sup <- Z_sup / sqrt(n1 + n2)

  sup_res <- list(p.value = p_sup, effect.size = r_sup)

  return (sup_res)
}

MiSA_orig_eachTaxon <- function(d, group.filtered, omnibus_r, omnibus_num){

  # Check whether d is a named matrix
  if (is.null(colnames(d))) {
    stop("'d' matrix in MiSA_data must have column names.")
  }

  # Ensure the first column of d is all zeros (representing r = 0);
  # if not, prepend a zero column named "0".
  if (!all(d[, 1] == 0)) {
    d <- cbind("0" = 0, d)
  }

  col_radii <- as.numeric(colnames(d))

  if (!is.null(omnibus_r)){
    # check whether the omnibus_r input is valid
    if (all(omnibus_r %in% col_radii)){
      # find the column indices in d
      idx <- match(omnibus_r, col_radii)
    } else{
      invalid_r <- omnibus_r[!omnibus_r %in% col_radii]
      cat("The following omnibus_r values were not found in the available radii:",
          paste(invalid_r, collapse = ", "), "- using default instead.\n")

      omnibus_num <- min(5, ncol(d)-1)
      idx <- round(seq(1, ncol(d), length.out= (omnibus_num+1))[-1])
    }
  } else {
    if (is.null(omnibus_num)){
      omnibus_num <- min(5, ncol(d)-1)
    } else{
      if (omnibus_num >= ncol(d)){
        cat("omnibus_num exceeds the number of available radii; using all available radii instead.\n")
        omnibus_num <- ncol(d)-1 # exclude col of r=0
      }
    }

    idx <- round(seq(1, ncol(d), length.out= (omnibus_num+1))[-1])
  }

  p.v <- sapply(idx, test_pure_dist, diff_mat=t(d), group=group.filtered)
  names(p.v) <- paste("r",colnames(d)[idx])

  p.pure_dist <- harmonicmeanp::p.hmp(p.v, L=length(p.v))

  # orig_res <- data.frame("hmp combined p"=p.pure_dist, t(p.v))
  orig_res <- c(p.pure_dist, p.v)

  return (orig_res)
}

MiSA_abs_eachTaxon <- function(d, group.filtered, omnibus_r, omnibus_num){
  # Check whether d is a named matrix
  if (is.null(colnames(d))) {
    stop("'d' matrix in MiSA_data must have column names.")
  }

  # Ensure the first column of d is all zeros (representing r = 0);
  # if not, prepend a zero column named "0".
  if (!all(d[, 1] == 0)) {
    d <- cbind("0" = 0, d)
  }

  col_radii <- as.numeric(colnames(d))

  if (!is.null(omnibus_r)){
    # check whether the omnibus_r input is valid
    if (all(omnibus_r %in% col_radii)){
      # find the column indices in d
      idx <- match(omnibus_r, col_radii)
    } else{
      invalid_r <- omnibus_r[!omnibus_r %in% col_radii]
      cat("The following omnibus_r values were not found in the available radii:",
          paste(invalid_r, collapse = ", "), "- using default instead.\n")

      omnibus_num <- min(5, ncol(d)-1)
      idx <- round(seq(1, ncol(d), length.out= (omnibus_num+1))[-1])
    }
  } else {
    if (is.null(omnibus_num)){
      omnibus_num <- min(5, ncol(d)-1)
    } else{
      if (omnibus_num >= ncol(d)){
        cat("omnibus_num exceeds the number of available radii; using all available radii instead.\n")
        omnibus_num <- ncol(d)-1 # exclude col of r=0
      }
    }

    idx <- round(seq(1, ncol(d), length.out= (omnibus_num+1))[-1])
  }

  p.v2 <- sapply(idx, test_abs_dist, diff_mat=t(d), group=group.filtered)
  names(p.v2) <- paste("r",colnames(d)[idx])

  p.abs_dist <- harmonicmeanp::p.hmp(p.v2, L=length(p.v2))

  #abs_res <- data.frame("hmp combined p"=p.abs_dist, t(p.v2))
  abs_res <- c(p.abs_dist, p.v2)

  return (abs_res)
}



