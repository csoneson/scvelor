#' Wrapper for scVelo
#'
#' This function provides a wrapper for running the scVelo pipeline
#' to estimate RNA velocities based on spliced and unspliced counts.
#' The following sequence of functions will be run, after the creation
#' of an AnnData object:
#' * pp.filter_and_normalize
#' * pp.moments
#' * tl.recover_dynamics (only for the dynamical model)
#' * tl.velocity
#' * tl.velocity_graph
#' * tl.velocity_pseudotime
#' * tl.latent_time (only for the dynamical model)
#' * tl.velocity_confidence
#'
#' @param spliced,unspliced Count matrices with spliced and unspliced
#'   counts. These matrices must have the same number of rows (representing
#'   genes) and columns (representing cells).
#' @param output_anndata Logical, whether to return the AnnData object
#'   with all estimated values.
#' @param params List of lists, allowing the user to provide arguments to
#'   the functions in scVelo. The elements of \code{params} are named
#'   after the respective scVelo function (e.g., \code{filter_and_normalize},
#'   see full list above). Each element is a list, defining the arguments
#'   to provide to the function. The only argument that must be defined
#'   is \code{params$velocity$mode}, which has to be one of 'steady_state',
#'   'stochastic' or 'dynamical'.
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @references
#' Bergen et al (2019)
#'
#' @importFrom basilisk basiliskRun basiliskStart basiliskStop
#' @importFrom reticulate import
#'
scvelor <- function(spliced, unspliced,
                    output_anndata=FALSE,
                    params=list(velocity=list(mode="stochastic"))) {
  if (is.null(params$velocity$mode)) {
    stop("The velocity mode must be specified (in params$velocity$mode)")
  }

  proc <- basilisk::basiliskStart(scvelo)
  on.exit(basilisk::basiliskStop(proc))

  scvres <- basilisk::basiliskRun(env=scvelo, fun=function(spliced, unspliced, params) {
    ## Import python modules
    and <- reticulate::import("anndata")
    scv <- reticulate::import("scvelo")

    ## Create AnnData object (note that matrices should have cells as rows)
    adata <- and$AnnData(t(spliced),
                         layers=list(spliced=t(spliced), unspliced=t(unspliced)))
    adata$obs_names <- colnames(spliced)
    adata$var_names <- rownames(spliced)

    ## Filter and normalize
    do.call(scv$pp$filter_and_normalize, c(list(data=adata),
                                           params$filter_and_normalize))

    ## Estimate moments
    do.call(scv$pp$moments, c(list(data=adata),
                              params$moments))

    ## Estimate velocity
    if (params$velocity$mode=="dynamical") {
      do.call(scv$tl$recover_dynamics, c(list(data=adata),
                                         params$recover_dynamics))
    }
    do.call(scv$tl$velocity, c(list(data=adata),
                               params$velocity))

    ## Build velocity graph
    do.call(scv$tl$velocity_graph, c(list(data=adata),
                                     params$velocity_graph))

    ## Estimate pseudotime
    do.call(scv$tl$velocity_pseudotime, c(list(adata=adata),
                                          params$velocity_pseudotime))

    ## Estimate latent time
    if (params$velocity$mode=="dynamical") {
      do.call(scv$tl$latent_time, c(list(data=adata),
                                    params$latent_time))
    }

    ## Estimate velocity confidence
    do.call(scv$tl$velocity_confidence, c(list(data=adata),
                                          params$velocity_confidence))

    ## Return values
    output <- list(obs=adata$obs,
                   var=adata$var)
    if (output_anndata) {
      output$adata=adata
    }
    output
  }, spliced=spliced, unspliced=unspliced, params=params)

  scvres
}
