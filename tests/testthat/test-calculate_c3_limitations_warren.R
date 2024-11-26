# Get test curves to use
source('one_curve_c3_aci.R')

test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- fit_c3_aci(
        one_curve_bad,
        Ca_atmospheric = 420,
        optim_fun = optimizer_nmkb(1e-7)
    )

    limit_res_bad <- expect_no_error(
        calculate_c3_limitations_warren(fit_res_bad$fits)
    )

    expect_true(all(is.na(limit_res_bad[, 'lm_warren'])))
    expect_true(all(is.na(limit_res_bad[, 'An_inf_gmc'])))
})

test_that('fit results have not changed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        optim_fun = optimizer_nmkb(1e-7),
        fit_options = list(alpha_old = 0, gmc = 1.0),
        calculate_confidence_intervals = FALSE,
        remove_unreliable_param = 0
    )

    limit_res <- expect_silent(
        calculate_c3_limitations_warren(fit_res$fits)
    )

    expect_equal(
        as.numeric(limit_res[1, c('Cc_inf_gmc', 'Cc_inf_gsc', 'An_inf_gmc', 'An_inf_gsc', 'lm_warren', 'ls_warren')]),
        c(38.20251, 29.80533, -5.91612, -8.09788, 0.20223, 0.41717),
        tolerance = 1e-5
    )
})
