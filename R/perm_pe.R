#'Exact permutation p-values
#'
#'Calculates exact p-values for permutation tests when permutations are randomly
#'drawn with replacement. This implementation is based on
#'(slightly adapted) the implementation by Belinda Phipson and Gordon Smyth
#'from the R package \code{statmod}
#'
#'
#'@author Belinda Phipson and Gordon Smyth (adapted by Boris Hejblum)
#'
#'@param nPermSupObs number of permutations that yielded test statistics at
#'least as extreme as the observed data. Can be a vector or an array of values.
#'@param nPermEff number of permutations effectively computed.
#'@param totalPossibleNPerm total number of permutations possible.
#'
#'@references Phipson B, and Smyth GK (2010). Permutation p-values should never
#'be zero: calculating exact p-values when permutations are randomly drawn.
#'\emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9,
#'Issue 1, Article 39.
#'\url{http://www.statsci.org/smyth/pubs/PermPValuesPreprint.pdf}
#'
#'@importFrom statmod gauss.quad.prob
#'@importFrom stats pbinom
#'@seealso statmod::permp
#'
#'@return a vector (or an array, similar to \code{nperm_supobs})
#'of exact p-values
#'
#'@examples
#'permPvals(10, 100, 1000)
#'
#'@export

permPvals <- function(nPermSupObs, nPermEff, totalPossibleNPerm) {

    ## (slightly) adapted from the implementation by Belinda Phipson and
    ## Gordon Smyth from the R package statmod

    nbnodes_quad <- 128
    z <- statmod::gauss.quad.prob(n = nbnodes_quad, l = 0,
                                  u = 0.5/totalPossibleNPerm)
    npvals <- length(nPermSupObs)
    prob <- rep(z$nodes, npvals)
    x2 <- rep(nPermSupObs, each = nbnodes_quad)
    Y <- matrix(data = stats::pbinom(x2, prob = prob, size = nPermEff),
                nrow = nbnodes_quad,
                ncol = npvals)
    int <- 0.5/totalPossibleNPerm * colSums(z$weights * Y)

    pe <- (nPermSupObs + 1)/(nPermEff + 1) - int

    return(pe)
}
