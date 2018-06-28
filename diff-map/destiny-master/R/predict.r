#' @include diffusionmap.r
NULL

#' Predict new data points using an existing DiffusionMap. The resulting matrix can be used in \link[=plot.DiffusionMap]{the plot method for the DiffusionMap}
#' 
#' @param dm        A \code{\link{DiffusionMap}} object.
#' @param new_data  New data points to project into the diffusion map.
#'                  Can be a \link[base]{matrix}, \link[base]{data.frame},
#'                  \link[Biobase:class.ExpressionSet]{ExpressionSet},
#'                  or \link[SingleCellExperiment]{SingleCellExperiment}.
#' @param verbose   Show progress messages?
#' @param ...       Passed to \code{\link[proxy:dist]{proxy::dist}(new_data, data, dm@distance, ...)}.
#' 
#' @return A \eqn{nrow(new\_data) \times ncol(eigenvectors(dif))} matrix of projected diffusion components for the new data.
#' 
#' @examples
#' data(guo)
#' g1 <- guo[, guo$num_cells != 32L]
#' g2 <- guo[, guo$num_cells == 32L]
#' dm <- DiffusionMap(g1)
#' dc2 <- dm_predict(dm, g2)
#' plot(dm, new_dcs = dc2)
#' 
#' @importFrom methods is
#' @importFrom Matrix Diagonal colSums rowSums
#' @importFrom proxy dist
#' @export
dm_predict <- function(dm, new_data, ..., verbose = FALSE) {
	if (!is(dm, 'DiffusionMap')) stop('dm needs to be a DiffusionMap')
	
	data <- dataset_extract_doublematrix(dataset(dm), dm@vars)
	new_data <- dataset_extract_doublematrix(new_data, dm@vars)
	if (!dataset_n_features(data) == dataset_n_features(new_data))
		stop('new data needs to have the same features as the one used to create the diffusion map')
	
	censor <- !is.null(dm@censor_val) || any(is.na(data)) || any(is.na(new_data))
	#censor imples euclidean distance
	
	#TODO: local sigma
	sigma <- optimal_sigma(dm)
	if (!censor) {
		if (verbose) cat('Creating distance matrix without censoring\n')
		measure <- switch(dm@distance, euclidean = 'Euclidean', cosine = 'correlation', rankcor = 'rankcor', stop('Unknown distance measure'))
		d <-   # no auto-simil2dist to not use 1-abs(simil)
			if (measure == 'Euclidean') proxy::dist(new_data, data, measure, ...)
			else 1 - proxy::simil(new_data, data, measure, ...)
		d2 <- unclass(d ^ 2) # matrix (dense)
		stopifnot(nrow(d2) == nrow(new_data))
		stopifnot(ncol(d2) == nrow(data))
		
		#TODO: zeros not on diag
		trans_p <- exp(-d2 / (2 * sigma ^ 2))
		trans_p[d2 == 0] <- 0
	} else {
		if (verbose) cat('Creating distance matrix with censoring\n')
		trans_p <- predict_censoring(data, new_data, dm@censor_val, dm@censor_range, dm@missing_range, sigma)
	}
	
	#trans_p:  columns: old data, rows: new data
	d_new <- rowSums(trans_p, na.rm = TRUE)
	norm_p <- get_norm_p(trans_p, dm@d, d_new, dm@density_norm)
	rm(trans_p)  # free memory
	
	d_norm_new <- rowSums(norm_p)
	
	d_rot     <- Diagonal(x = dm@d_norm  ^ -.5)
	d_rot_new <- Diagonal(x = d_norm_new ^ -.5)
	M_new <- d_rot_new %*% norm_p %*% d_rot
	
	if (verbose) cat('Transforming data\n')
	t(t(M_new %*% eigenvectors(dm)) / eigenvalues(dm))
}
