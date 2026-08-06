# The causalgenerics accessors are re-exported so that inspecting a bw weight
# vector's estimand and causal-weight status needs no second attachment. The
# shared ipw() generic is re-exported so that a fitted balancing object drives
# the effect-estimation workflow with an unqualified call after
# `library(balancing)`. The two generics that move an ipw result between its
# marginal and conditional readings join it for the same reason: the result
# comes from an unqualified call, so it is flipped by one too, and so is the
# entrypoint that pools a set of them across multiply imputed datasets.

#' @importFrom causalgenerics ipw
#' @export
causalgenerics::ipw

#' @importFrom causalgenerics as_marginal
#' @export
causalgenerics::as_marginal

#' @importFrom causalgenerics as_conditional
#' @export
causalgenerics::as_conditional

#' @importFrom causalgenerics pool_ipw
#' @export
causalgenerics::pool_ipw

#' @importFrom causalgenerics is_causal_wt
#' @export
causalgenerics::is_causal_wt

#' @importFrom causalgenerics estimand
#' @export
causalgenerics::estimand
