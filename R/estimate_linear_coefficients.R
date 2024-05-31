

estimate_linear_coefficient <- function(Y, design_matrix, method = c("linear", "mean", "cluster_median", "zero")){
  method <- match.arg(method)

  if(method == "linear"){
    linear_fit <- lm.fit(design_matrix, t(Y))
    t(linear_fit$coefficients)
  }else if(method == "mean"){
    # Check for intercept column / columns
    ones <- rep(1, nrow(design_matrix))
    intercept_fit <- lm.fit(design_matrix, ones)
    if(! sum(intercept_fit$residuals^2) < 1e-12){
      stop("The design matrix does not have an intercept. Cannot apply a single mean offset. Please change",
           "'linear_coefficient_estimator' to 'linear', 'cluster_median', or 'zero'.")
    }
    means <- MatrixGenerics::rowMeans2(Y)
    matrix(means, ncol = 1) %*% matrix(intercept_fit$coefficients, nrow = 1)
  }else if(method == "zero"){
    matrix(0, nrow = nrow(Y), ncol = ncol(design_matrix))
  }else if(method == "cluster_median"){
    min_cluster_membership <- 0.01
    pca <- pca(Y, n = 20)
    harm_obj <- harmony_init(pca$embedding, design_matrix, nclust = 30, verbose = FALSE)
    harm_obj <- harmony_max_div_clustering(harm_obj)
    Yt <- as.matrix(t(Y))
    coef <- do.call(cbind, lapply(seq_len(nrow(harm_obj$R)), \(cl){
      threshold <- min(min_cluster_membership, max(harm_obj$R) * 0.5)
      sel <- harm_obj$R[cl, ] > threshold
      tryCatch({
        fit <- lm.wfit(design_matrix[sel,,drop=FALSE], y = Yt[sel,,drop=FALSE], w = harm_obj$R[cl, sel])
        as.numeric(t(fit$coefficients))
      }, error = function(e){
        rep(NA_real_, nrow(Y) * ncol(design_matrix))
      })
    }))
    wmed <- matrixStats::rowWeightedMedians(coef, w = rowSums(harm_obj$R), na.rm = TRUE, interpolate = FALSE)
    matrix(wmed, nrow = nrow(Y), ncol = ncol(design_matrix))
  }
}


