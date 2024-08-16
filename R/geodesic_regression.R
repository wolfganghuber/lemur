
#######################
# Grassmann Manifold  #
#######################


#' Solve d(P, exp_p(V * x))^2 for V
#'
#' @returns A three-dimensional array with the coefficients `V`.
#'
#' @keywords internal
grassmann_geodesic_regression <- function(coordsystems, design, base_point, weights = 1, tangent_regression = FALSE){
  # Validate input
  n_obs <- nrow(design)
  n_coef <- ncol(design)
  n_amb <- nrow(base_point)
  n_emb <- ncol(base_point)

  coordsystems <- if(is.list(coordsystems)){
    coordsystems
  }else if(is.array(coordsystems)){
    stopifnot(length(dim(coordsystems)) == 3)
    destack_slice(coordsystems)
  }else{
    stop("Cannot handle coordsystems of type: ", toString(class(coordsystems), width = 100))
  }
  stopifnot(length(coordsystems) == n_obs)
  stopifnot(all(vapply(coordsystems, \(emb) nrow(emb) == n_amb && ncol(emb) == n_emb, FUN.VALUE = logical(1L))))
  # stopifnot(all(vapply(coordsystems, \(emb) is_grassmann_element(emb), FUN.VALUE = logical(1L))))
  weights <- rep_len(weights, n_obs)



  # Initialize with tangent regression (if possible)
  tangent_vecs <- lapply(coordsystems, \(emb) as.vector(grassmann_log(base_point, emb)))
  merged_vecs <- stack_cols(tangent_vecs)
  tangent_fit <- if(nrow(merged_vecs) == 0){
    matrix(nrow = 0, ncol = ncol(design))
  }else{
    t(lm.wfit(design, t(merged_vecs), w = weights)$coefficients)
  }
  coef <- stack_slice(lapply(seq_len(ncol(tangent_fit)), \(idx) matrix(tangent_fit[,idx], nrow = n_amb, ncol = n_emb)))
  dimnames(coef) <- list(NULL, NULL, colnames(tangent_fit))


  if(tangent_regression){
    coef
  }else{
    # warning("Refine regression using Riemannian optimization. (Not yet implemented)")
    coef
  }
}

#' Solve ||Y - exp_p(V * x) Y ||^2_2 for V
#'
#' @returns A three-dimensional array with the coefficients `V`.
#'
#' @keywords internal
grassmann_lm <- function(data, design, base_point, tangent_regression = FALSE){
  nas <- apply(data, 2, anyNA) | apply(design, 1, anyNA)
  data <- data[,!nas,drop=FALSE]
  design <- design[!nas,,drop=FALSE]

  n_obs <- nrow(design)
  n_coef <- ncol(design)
  n_amb <- nrow(base_point)
  n_emb <- ncol(base_point)

  # Initialize with tangent regression
  mm_groups <- get_groups(design)
  groups <- unique(mm_groups)
  reduced_design <- mply_dbl(groups, \(gr) design[which(mm_groups == gr)[1],], ncol = ncol(design))
  if(any(table(mm_groups) < n_emb)){
    problematic_mat <- cbind(n_occurrences = c(table(mm_groups)), reduced_design)
    stop("Too few datapoints in some design matrix group.\n\n", glmGamPoi:::format_matrix(problematic_mat),
         "\nEach row must occurr at least n_embedding=", n_emb, " times.\n")
  }
  group_planes <- lapply(groups, \(gr) pca(data[,mm_groups == gr,drop=FALSE], n = n_emb, center = FALSE)$coordsystem)
  group_sizes <- vapply(groups, \(gr) sum(mm_groups == gr), FUN.VALUE = 0L)
  coef <- grassmann_geodesic_regression(group_planes, design = reduced_design, base_point = base_point, weights = group_sizes, tangent_regression = TRUE)
  if(tangent_regression){
    coef
  }else{
    # warning("Refine regression using Riemannian optimization. (Not yet implemented)")
    coef
  }
}


get_groups <- function (design_matrix) {
  vctrs::vec_group_id(unclass(design_matrix))
}




