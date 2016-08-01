#' Log-catches-at-age of Gulf menhaden in the U.S. Gulf of Mexico, 1964-2004.
#'
#' A data set containing Y = log(X + 0.01), where X is the estimated annual 
#' catch-at-age, in million fish, of Brevoortia patronus, from the reduction
#' fishery. Derived from Table 2 of Vaughan et al., 2007, Fish. Res. 83, 263-275.
#'
#' @format A data frame with 41 rows and 6 columns
#' \describe{
#'   \item{age0}{annual estimated log-catch at age zero}
#'   \item{age1}{annual estimated log-catch at age one}
#'   ...
#' }
#' @source \url{http://www.sciencedirect.com/science/article/pii/S0165783606003638}
"Y"

#' Nominal fishing effort of the gulf menhaden reduction fishery in the U.S.
#' Gulf of Mexico, 1964-2004. 
#'
#' A data set containing nominal effort, in 1000 vessel-tonne-weeks, associated
#' with the US menhaden reduction fishery in the Gulf of Mexico.
#' Extracted from Table 2 of Vaughan et al., 2007, Fish. Res. 83, 263-275.
#'
#' @format A data frame with 41 rows and 1 columns
#' \describe{
#'   \item{effort}{annual estimated nominal effort}
#' }
#' @source \url{http://www.sciencedirect.com/science/article/pii/S0165783606003638}
"E"