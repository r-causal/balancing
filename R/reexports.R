# Re-export the broom-ecosystem generics so that tidy() and glance() are
# available with `library(balancing)` alone. balancing registers methods on both
# generics for its result class. The propensity accessors are re-exported so
# that inspecting a bw weight vector's estimand and causal-weight status needs no
# second attachment.

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom propensity is_causal_wt
#' @export
propensity::is_causal_wt

#' @importFrom propensity estimand
#' @export
propensity::estimand
