## ==========================================================================
## transformList allow us to programmatically access transformations.
## They normaly contain items of class transformMap.
## ==========================================================================






## ==========================================================================
## colnames method: This gives us the parameter names we want to transform
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("colnames",
          signature=signature(x="transformList"),
          definition=function(x, do.NULL=TRUE, prefix="col")
      {
          unique(sapply(x@transforms, slot, "input"))
      })



## ==========================================================================
## Concatenate two transformLists (or one list and a transformList)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("c",
          signature=signature(x="transformList"),
          definition=function(x, ..., recursive=FALSE)
      {
          ## Try to coerce everyone to a transformList first
          all.t <- lapply(list(...), as, "transformList")
          params <- c(sapply(x@transforms, slot, "output"),
                      unlist(lapply(all.t, function(x)
                                    sapply(x@transforms, slot, "output"))))
          if(length(params) > length(unique(params)))
              stop("All output parameters must be unique when combining ",
                   "transforms.", call.=FALSE)
          new("transformList", transforms=c(x@transforms,
                               unlist(sapply(all.t, slot, "transforms"))))
      })


## ==========================================================================
## Get/Set of functions of a transformList or a transformMap
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Get the parameters of the function called by the transformation
#'
#' @param trans a transformation list (such as the value returned by
#'   transformList or estimateLogicle functions) or a transform map (such as the
#'   value returned by arcsinhTransform or logicleTransform functions).
#' @returns a list of parameters  of the function called by each  transformation
#'   of trans.
#' @export
transform_getParams <- function(trans)
{
  getParams <- function(transMap) {
    func <- transMap@f
    env_func <- environment(func)
    vars <- ls(env_func)
    vals <- sapply(vars, function(v) get(v, envir = env_func), simplify = FALSE)
    return(vals)
  }
  if (inherits(trans, "transformMap"))
    return(getParams(trans))
  if (inherits(trans, "transformList"))
    return(lapply(trans@transforms, getParams))
  stop("Not implemented for an object of class ", class(trans))
}

#' Set the parameters of the function called by the transformation
#'
#' @param trans a transformation list (such as the value returned by
#'   transformList or estimateLogicle functions) or a transform map (such as the
#'   value returned by arcsinhTransform or logicleTransform functions).
#' @param params a list of parameters as returned by transfom_getParams
#'   function.
#' @returns a transformation list or a single transformation map according to
#'   trans.
#' @export
transform_setParams <- function(trans, params)
{
  setParams <- function(transMap, params_map) {
    func <- transMap@f
    env_func <- environment(func)
    vars <- ls(env_func)
    varp <- names(params_map)
    comm <- intersect(vars, varp)
    if (length(comm)) {
      # make a copy of the environment to avoid overwriting
      env_copy <- as.environment(as.list(env_func, all.names=TRUE))
      # update the copy
      for (v in comm) {
        assign(v, unname(params_map[v]), envir = env_copy)
      }
      # assign the copy to the function
      environment(transMap@f) <- env_copy
    }
    return(transMap)
  }
  if (inherits(trans, "transformMap")) {
    # TODO: check params match trans
    return(setParams(trans, params))
  }
  if (inherits(trans, "transformList")) {
    # TODO: check params match trans
    trans_maps <- trans@transforms
    trans_mrks <- names(trans_maps)
    params_mrks <- names(params)
    comm_mrks <- intersect(trans_mrks, params_mrks)
    if (length(comm_mrks)) {
      for (mrk in comm_mrks) {
        trans_maps[[mrk]] <- setParams(trans_maps[[mrk]], params[[mrk]])
      }
    }
    trans@transforms <- trans_maps
    return(trans)
  }
  stop("Not implemented for an object of class ", class(trans))
}
