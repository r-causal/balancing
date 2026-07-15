# The propensity accessors are re-exported so that inspecting a bw weight
# vector's estimand and causal-weight status needs no second attachment. The
# shared ipw() generic is re-exported so that a fitted balancing object drives
# the effect-estimation workflow with an unqualified call after
# `library(balancing)`.

#' @importFrom causalgenerics ipw
#' @export
causalgenerics::ipw

#' @importFrom propensity is_causal_wt
#' @export
propensity::is_causal_wt

#' @importFrom propensity estimand
#' @export
propensity::estimand
