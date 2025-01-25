#' Computes a transform using the 'logicle_transform' function but keeping zero
#' at zero
#'
#' Using the logicle transformation, the zero intensity is transformed into w
#' which is usually not zero. The logicleZ removes this offset after the
#' standard logicle transformation in order to keep zero at zero. See Logicle
#' transformation \code{\link{logicleTransform}} for all parameters. By default,
#' the offset (z_post) is set the value of w.
#'
#' The logicleZ allows shifting input intensity before transformation and
#' scaling transformed intensity. The pre-shifting coefficient (z_pre) aims to
#' bring the center of the negative peak at zero before any computation. The
#' user can estimate it as the raw intensity of the center of the negative peak.
#' z_pre is zero by default (no shifting). The post-multiplicative coefficient
#' (m_post) is one by default (no post-scaling).
#'
#' The estimateLogicle function estimates the parameter w (width of the linear
#' zone) from the negative (left) part of a marker's negative peak. It
#' calculates the 5th percentile of intensities below the zero threshold. The
#' estimateLogicleZ function allows to specify the center of the negative peak
#' and subtract this value from all input intensities before performing the
#' classic estimate. The center of the negative peak can be introduced in the
#' calculations either as z_est or as z_pre (previously defined). When
#' specifying z_est, z_est is subtracted from raw intensities only for
#' estimating w. When specifying z_pre, z_pre is subtracted from raw intensities
#' for estimating w and during the transformation. NB: z_pre and z_est should
#' not be both different from zero.
#'
#' @aliases logicleZTransform estimateLogicleZ
#' @usage logicleZTransform(transformationId="defaultLogicleZTransform", w =
#'   0.5, t = 262144, m = 4.5, a = 0)
#' @usage estimateLogicle(x, channels, ...)
#' @param transformationId A name to assign to the transformation. Used by the
#'   transform/filter routines.
#' @param w w of logicleTransform.
#' @param t t of logicleTransform.
#' @param m m of logicleTransform.
#' @param a a of logicleTransform.
#' @param z_pre intensity to remove before any computation to shift the peak at
#'   zero. This is 0 by default.
#' @param z_post value to remove after the transformation to bring zero back at
#'   zero. This is w by default.
#' @param m_post multiplicative coefficient to apply on the transformed
#'   intensities. This is 1 by default.
#' @param x Input flow frame for which the logicleZ transformations are to be
#'   estimated.
#' @param channels channels or markers for which the logicleZ transformation is
#'   to be estimated.
#' @param ... other arguments:
#'
#'   z_est intensity of the zero for estimating the negative part of the
#'   negative peak. This is 0 by default. NB: z_pre and z_est should not be both
#'   different from zero.
#'
#'   q: a numeric type specifying quantile value, default is 0.05.
#' @seealso \code{\link[flowCore]{inverseLogicleZTransform}},
#'   \code{\link[flowCore]{estimateLogicle} }
#' @references Parks D.R., Roederer M., Moore W.A.(2006) A new "logicle" display
#'   method avoids deceptive effects of logarithmic scaling for low signals and
#'   compensated data. CytometryA, 96(6):541-51.
#' @keywords methods
#' @examples
#'
#' # GvHD is not a recent dataset to show the interest of this transformation
#'
#' @export
logicleZTransform <- function(
    transformationId="defaultLogicleZTransform",
    w = 0.5, t = 262144, m = 4.5, a = 0, z_pre = 0, z_post = w, m_post = 1)
{
  new("transform", .Data = function(x) {
    if (z_pre) x <- x - z_pre
    x <- logicle_transform(as.double(x), as.double(t), as.double(w), as.double(m), as.double(a), FALSE)
    if (z_post) x <- x - z_post
    if (m_post) x <- x * m_post
    return(x)
  },
  transformationId = transformationId)
}
#' @export
estimateLogicleZ <- function(x, channels, ...)
  UseMethod("estimateLogicleZ")
#' @export
estimateLogicleZ.flowFrame <- function(x, channels, ...)
{
  trans <- .estimateLogicleZ(x, channels, ...)
  channels <- names(trans)
  transformList(channels, trans)
}
.estimateLogicleZ <- function(x, channels, ...)
{
  if(!is(x,"flowFrame")&&!is(x,"cytoframe"))
    stop("x has to be an object of class \"flowFrame\"")
  if(missing(channels))
    stop("Please specify the channels to be transformed")
  channels <- sapply(channels, function(channel)
    getChannelMarker(x, channel)[["name"]], USE.NAMES = FALSE)
  sapply(channels, function(p) {
    .lgclZTrans(x, p, ...)
  })
}

#' It is mainly trying to estimate w (linearization width in asymptotic decades)
#' value based on given m and data range
#' @param obj flowFrame
#' @param p channel name
#' @param m full length of transformed display in decodes
#' @param t top of the scale of data value
#' @param a additional negative range to be included in display in decades
#' @param q quantile of negative data value (used to adjust w calculation)
#' @param type character either "instrument" or "data". The data range.
#' @param z_pre intensity to remove before any computation to shift the peak at
#'   zero; this changes the data and its range. NB: z_pre and z_est should not
#'   be both different from zero.
#' @param z_est intensity of the zero for estimating the negative part of the
#'   negative peak; this solely changes the data for the estimation of w. NB:
#'   z_pre and z_est should not be both different from zero.
#' @noRd
.lgclZTrans  <- function(
    obj, p, t , m, a = 0, q = 0.05, type = "instrument",
    z_pre = 0, z_est = 0)
{
  type <- match.arg(type, c("instrument", "data"))
  transId <- paste(p, "logicleZTransform", sep = "_")

  rng <- range(obj)
  dat <- exprs(obj)[,p]

  # shift the intensity to bring the negative peak at zero; this changes the
  # data range
  if (z_pre) dat <- dat - z_pre

  # standard estimation
  if(missing(t)){
    if(type == "instrument")
      t <- rng[,p][2]
    else
      t <- max(dat)
  }

  if(missing(m)){
    if(type == "instrument")
      m <- 4.5
    else
      m <- log10(t) + 1
  }

  # shift the intensity to bring the negative peak at zero; this solely changes
  # the data for the estimation of w
  if (z_est) dat <- dat - z_est

  # standard estimation
  dat <- dat[dat<0]
  w <- 0
  if(length(dat)) {
    r <- .Machine$double.eps + unname(quantile(dat, q))
    w <- (m-log10(t/abs(r))) / 2
    if(w < 0)
      stop("w is negative! Try to increase 'm'")
  }
  # the logicleZ sets z_post = w by default
  logicleZTransform(transformationId = transId, w = w, t = t, m = m, a = a, z_pre = z_pre)
}


#' Computes the inverse of the transform defined by the 'logicleZTransform'
#' function or the transformList generated by 'estimateLogicleZ' function
#'
#' inverseLogicleZTransform can be use to compute the inverse of the LogicleZ
#' transformation. The parameters w, t, m, a for calculating the inverse are
#' obtained from the 'trans' input passed to the 'inverseLogicleZTransform'
#' function.
#'
#' @usage inverseLogicleZTransform(trans, transformationId, ...)
#' @param trans An object of class 'transform' created using the
#' 'logicleZTransform' function or class 'transformList' created by
#' 'estimateLogicleZ'.
#' @param transformationId A name to assigned to the inverse transformation.
#' Used by the transform routines.
#' @param ...  not used.
#' @seealso \code{\link[flowCore]{logicleZTransform}}
#' @keywords methods
#' @examples
#'
#' # GvHD is not a recent dataset to show the interest of this transformation
#'
#' @export
inverseLogicleZTransform <- function(trans, transformationId, ...)
  UseMethod("inverseLogicleZTransform")
#' @export
inverseLogicleZTransform.default <- function(trans, transformationId, ...)
{
  stop("trans has to be an object of class \"transform\" created using the ",
       "\"logicleZTransform\" function\n or a 'transformList' created by ",
       "'estimateLogicleZ'\n")
}
#' @export
inverseLogicleZTransform.transform <- function(trans, transformationId, ...)
{
  k <- .inverseLogicleZTransform(trans@.Data)
  if(missing(transformationId))
    k@transformationId <- paste( "inverse", trans@transformationId, sep ="_")
  k
}
.inverseLogicleZTransform <- function(func)
{
  # check needed parameters
  pars <- c("w", "t", "m", "a", "z_pre", "z_post", "m_post")
  if(!all(pars %in% ls(environment(func))))
    stop("\"trans\" is not a valid object produced using the \"logicleZ\" function")
  # copy parameters to the current environment
  for (v in pars)
    assign(v, environment(func)[[v]])
  # clean up current environment as it is inherited
  rm(pars, func)
  # return the object without creating a variable
  new("transform", .Data = function(x) {
    if (m_post) x <- x / m_post
    if (z_post) x <- x + z_post
    x <- logicle_transform(as.double(x), as.double(t), as.double(w), as.double(m), as.double(a), TRUE)
    if (z_pre) x <- x + z_pre
    return(x)
  })
}
#' @export
inverseLogicleZTransform.transformList <- function(trans, transformationId, ...)
{
  invs <- sapply(trans@transforms, function(obj){
    .inverseLogicleZTransform(obj@f)
  })
  channels <- names(invs)
  if(missing(transformationId))
    transformationId <- paste( "inverse", trans@transformationId, sep ="_")
  transformList(channels, invs, transformationId = transformationId)
}
