
// clang-format sorts includes unless SortIncludes: Never. However, the ordering
// does matter here. So, we need to disable clang-format for safety.

// clang-format off
#include <stdint.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>
// clang-format on

#include "rust/api.h"

static uintptr_t TAGGED_POINTER_MASK = (uintptr_t)1;

SEXP handle_result(SEXP res_) {
    uintptr_t res = (uintptr_t)res_;

    // An error is indicated by tag.
    if ((res & TAGGED_POINTER_MASK) == 1) {
        // Remove tag
        SEXP res_aligned = (SEXP)(res & ~TAGGED_POINTER_MASK);

        // Currently, there are two types of error cases:
        //
        //   1. Error from Rust code
        //   2. Error from R's C API, which is caught by R_UnwindProtect()
        //
        if (TYPEOF(res_aligned) == CHARSXP) {
            // In case 1, the result is an error message that can be passed to
            // Rf_errorcall() directly.
            Rf_errorcall(R_NilValue, "%s", CHAR(res_aligned));
        } else {
            // In case 2, the result is the token to restart the
            // cleanup process on R's side.
            R_ContinueUnwind(res_aligned);
        }
    }

    return (SEXP)res;
}

SEXP savvy_eval_psi_cbps__impl(SEXP c_arg__coefs, SEXP c_arg__covs, SEXP c_arg__treat, SEXP c_arg__s_weights, SEXP c_arg__estimand, SEXP c_arg__link) {
    SEXP res = savvy_eval_psi_cbps__ffi(c_arg__coefs, c_arg__covs, c_arg__treat, c_arg__s_weights, c_arg__estimand, c_arg__link);
    return handle_result(res);
}

SEXP savvy_eval_psi_entropy__impl(SEXP c_arg__coefs, SEXP c_arg__covs, SEXP c_arg__group_idx, SEXP c_arg__targets, SEXP c_arg__base_weights, SEXP c_arg__s_weights, SEXP c_arg__n_eff, SEXP c_arg__esteq_scale) {
    SEXP res = savvy_eval_psi_entropy__ffi(c_arg__coefs, c_arg__covs, c_arg__group_idx, c_arg__targets, c_arg__base_weights, c_arg__s_weights, c_arg__n_eff, c_arg__esteq_scale);
    return handle_result(res);
}

SEXP savvy_eval_psi_ipt__impl(SEXP c_arg__coefs, SEXP c_arg__covs, SEXP c_arg__treat_idx, SEXP c_arg__focal, SEXP c_arg__s_weights, SEXP c_arg__estimand, SEXP c_arg__link) {
    SEXP res = savvy_eval_psi_ipt__ffi(c_arg__coefs, c_arg__covs, c_arg__treat_idx, c_arg__focal, c_arg__s_weights, c_arg__estimand, c_arg__link);
    return handle_result(res);
}

SEXP savvy_solve_cbps__impl(SEXP c_arg__covs_mod, SEXP c_arg__covs_bal, SEXP c_arg__treat, SEXP c_arg__s_weights, SEXP c_arg__estimand, SEXP c_arg__link, SEXP c_arg__over, SEXP c_arg__twostep, SEXP c_arg__options) {
    SEXP res = savvy_solve_cbps__ffi(c_arg__covs_mod, c_arg__covs_bal, c_arg__treat, c_arg__s_weights, c_arg__estimand, c_arg__link, c_arg__over, c_arg__twostep, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_cbps_cont__impl(SEXP c_arg__covs, SEXP c_arg__expo, SEXP c_arg__s_weights, SEXP c_arg__options) {
    SEXP res = savvy_solve_cbps_cont__ffi(c_arg__covs, c_arg__expo, c_arg__s_weights, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_cbps_multi__impl(SEXP c_arg__covs, SEXP c_arg__treat_idx, SEXP c_arg__focal, SEXP c_arg__s_weights, SEXP c_arg__estimand, SEXP c_arg__link, SEXP c_arg__options) {
    SEXP res = savvy_solve_cbps_multi__ffi(c_arg__covs, c_arg__treat_idx, c_arg__focal, c_arg__s_weights, c_arg__estimand, c_arg__link, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_energy__impl(SEXP c_arg__covs, SEXP c_arg__treat, SEXP c_arg__s_weights, SEXP c_arg__distance, SEXP c_arg__estimand, SEXP c_arg__improved, SEXP c_arg__moment_covs, SEXP c_arg__targets, SEXP c_arg__tols, SEXP c_arg__min_weight, SEXP c_arg__weight_penalty, SEXP c_arg__options) {
    SEXP res = savvy_solve_energy__ffi(c_arg__covs, c_arg__treat, c_arg__s_weights, c_arg__distance, c_arg__estimand, c_arg__improved, c_arg__moment_covs, c_arg__targets, c_arg__tols, c_arg__min_weight, c_arg__weight_penalty, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_energy_cont__impl(SEXP c_arg__covs, SEXP c_arg__treat, SEXP c_arg__s_weights, SEXP c_arg__distance, SEXP c_arg__dimension_adj, SEXP c_arg__min_weight, SEXP c_arg__weight_penalty, SEXP c_arg__d_covs, SEXP c_arg__d_treat, SEXP c_arg__bal_covs, SEXP c_arg__bal_tols, SEXP c_arg__options) {
    SEXP res = savvy_solve_energy_cont__ffi(c_arg__covs, c_arg__treat, c_arg__s_weights, c_arg__distance, c_arg__dimension_adj, c_arg__min_weight, c_arg__weight_penalty, c_arg__d_covs, c_arg__d_treat, c_arg__bal_covs, c_arg__bal_tols, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_energy_multi__impl(SEXP c_arg__covs, SEXP c_arg__treat_idx, SEXP c_arg__focal, SEXP c_arg__s_weights, SEXP c_arg__distance, SEXP c_arg__estimand, SEXP c_arg__improved, SEXP c_arg__moment_covs, SEXP c_arg__targets, SEXP c_arg__tols, SEXP c_arg__min_weight, SEXP c_arg__weight_penalty, SEXP c_arg__options) {
    SEXP res = savvy_solve_energy_multi__ffi(c_arg__covs, c_arg__treat_idx, c_arg__focal, c_arg__s_weights, c_arg__distance, c_arg__estimand, c_arg__improved, c_arg__moment_covs, c_arg__targets, c_arg__tols, c_arg__min_weight, c_arg__weight_penalty, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_entropy__impl(SEXP c_arg__covs, SEXP c_arg__group_idx, SEXP c_arg__targets, SEXP c_arg__base_weights, SEXP c_arg__s_weights, SEXP c_arg__tols, SEXP c_arg__n_eff, SEXP c_arg__options) {
    SEXP res = savvy_solve_entropy__ffi(c_arg__covs, c_arg__group_idx, c_arg__targets, c_arg__base_weights, c_arg__s_weights, c_arg__tols, c_arg__n_eff, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_entropy_cont__impl(SEXP c_arg__covs, SEXP c_arg__targets, SEXP c_arg__tols, SEXP c_arg__dist_ind, SEXP c_arg__base_weights, SEXP c_arg__s_weights, SEXP c_arg__n_eff, SEXP c_arg__options) {
    SEXP res = savvy_solve_entropy_cont__ffi(c_arg__covs, c_arg__targets, c_arg__tols, c_arg__dist_ind, c_arg__base_weights, c_arg__s_weights, c_arg__n_eff, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_ipt__impl(SEXP c_arg__covs, SEXP c_arg__treat, SEXP c_arg__s_weights, SEXP c_arg__estimand, SEXP c_arg__link, SEXP c_arg__options) {
    SEXP res = savvy_solve_ipt__ffi(c_arg__covs, c_arg__treat, c_arg__s_weights, c_arg__estimand, c_arg__link, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_ipt_multi__impl(SEXP c_arg__covs, SEXP c_arg__treat_idx, SEXP c_arg__focal, SEXP c_arg__s_weights, SEXP c_arg__estimand, SEXP c_arg__link, SEXP c_arg__options) {
    SEXP res = savvy_solve_ipt_multi__ffi(c_arg__covs, c_arg__treat_idx, c_arg__focal, c_arg__s_weights, c_arg__estimand, c_arg__link, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_sbw__impl(SEXP c_arg__treat, SEXP c_arg__s_weights, SEXP c_arg__estimand, SEXP c_arg__norm, SEXP c_arg__moment_covs, SEXP c_arg__targets, SEXP c_arg__tols, SEXP c_arg__min_weight, SEXP c_arg__options) {
    SEXP res = savvy_solve_sbw__ffi(c_arg__treat, c_arg__s_weights, c_arg__estimand, c_arg__norm, c_arg__moment_covs, c_arg__targets, c_arg__tols, c_arg__min_weight, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_sbw_cont__impl(SEXP c_arg__treat, SEXP c_arg__covs, SEXP c_arg__s_weights, SEXP c_arg__norm, SEXP c_arg__tols, SEXP c_arg__min_weight, SEXP c_arg__options) {
    SEXP res = savvy_solve_sbw_cont__ffi(c_arg__treat, c_arg__covs, c_arg__s_weights, c_arg__norm, c_arg__tols, c_arg__min_weight, c_arg__options);
    return handle_result(res);
}

SEXP savvy_solve_sbw_multi__impl(SEXP c_arg__treat_idx, SEXP c_arg__focal, SEXP c_arg__s_weights, SEXP c_arg__estimand, SEXP c_arg__norm, SEXP c_arg__moment_covs, SEXP c_arg__targets, SEXP c_arg__tols, SEXP c_arg__min_weight, SEXP c_arg__options) {
    SEXP res = savvy_solve_sbw_multi__ffi(c_arg__treat_idx, c_arg__focal, c_arg__s_weights, c_arg__estimand, c_arg__norm, c_arg__moment_covs, c_arg__targets, c_arg__tols, c_arg__min_weight, c_arg__options);
    return handle_result(res);
}

SEXP savvy_thread_info__impl(void) {
    SEXP res = savvy_thread_info__ffi();
    return handle_result(res);
}


static const R_CallMethodDef CallEntries[] = {
    {"savvy_eval_psi_cbps__impl", (DL_FUNC) &savvy_eval_psi_cbps__impl, 6},
    {"savvy_eval_psi_entropy__impl", (DL_FUNC) &savvy_eval_psi_entropy__impl, 8},
    {"savvy_eval_psi_ipt__impl", (DL_FUNC) &savvy_eval_psi_ipt__impl, 7},
    {"savvy_solve_cbps__impl", (DL_FUNC) &savvy_solve_cbps__impl, 9},
    {"savvy_solve_cbps_cont__impl", (DL_FUNC) &savvy_solve_cbps_cont__impl, 4},
    {"savvy_solve_cbps_multi__impl", (DL_FUNC) &savvy_solve_cbps_multi__impl, 7},
    {"savvy_solve_energy__impl", (DL_FUNC) &savvy_solve_energy__impl, 12},
    {"savvy_solve_energy_cont__impl", (DL_FUNC) &savvy_solve_energy_cont__impl, 12},
    {"savvy_solve_energy_multi__impl", (DL_FUNC) &savvy_solve_energy_multi__impl, 13},
    {"savvy_solve_entropy__impl", (DL_FUNC) &savvy_solve_entropy__impl, 8},
    {"savvy_solve_entropy_cont__impl", (DL_FUNC) &savvy_solve_entropy_cont__impl, 8},
    {"savvy_solve_ipt__impl", (DL_FUNC) &savvy_solve_ipt__impl, 6},
    {"savvy_solve_ipt_multi__impl", (DL_FUNC) &savvy_solve_ipt_multi__impl, 7},
    {"savvy_solve_sbw__impl", (DL_FUNC) &savvy_solve_sbw__impl, 9},
    {"savvy_solve_sbw_cont__impl", (DL_FUNC) &savvy_solve_sbw_cont__impl, 7},
    {"savvy_solve_sbw_multi__impl", (DL_FUNC) &savvy_solve_sbw_multi__impl, 10},
    {"savvy_thread_info__impl", (DL_FUNC) &savvy_thread_info__impl, 0},
    {NULL, NULL, 0}
};

void R_init_balancing(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

    // Functions for initialization, if any.

}
