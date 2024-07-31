# Get test curves to use
source('one_curve_c4_aci.R')

test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- expect_silent(
        fit_c4_aci_hyperbola(
            one_curve_bad,
            OPTIM_FUN = optimizer_nmkb(1e-7),
            hard_constraints = 2,
            calculate_confidence_intervals = TRUE
        )
    )

    expect_equal(unique(fit_res_bad$fits[, 'c4_assimilation_hyperbola_msg']), 'Ci must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c4_assimilation_hyperbola_msg'], 'Ci must be >= 0')
    expect_true(all(is.na(fit_res_bad$fits[, c('A_fit', 'Ainitial', 'Amax')])))
    expect_true(all(is.na(fit_res_bad$fits_interpolated[, c('An', 'Ainitial', 'Amax')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('c4_curvature', 'c4_slope', 'rL', 'Vmax', 'AIC')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('c4_curvature_upper', 'c4_slope_upper', 'rL_upper', 'Vmax_upper')])))
})

test_that('Ci limits can be bypassed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c4_aci_hyperbola(
            one_curve_bad,
            OPTIM_FUN = optimizer_nmkb(1e-7),
            hard_constraints = 0,
            calculate_confidence_intervals = TRUE
        )
    )

    expect_equal(unique(fit_res$fits[, 'c4_assimilation_hyperbola_msg']), '')
    expect_equal(fit_res$parameters[, 'c4_assimilation_hyperbola_msg'], '')
    expect_true(all(!is.na(fit_res$fits[, c('A_fit')])))
})

test_that('fit results have not changed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci_hyperbola(
        one_curve,
        OPTIM_FUN = optimizer_nmkb(1e-7),
        calculate_confidence_intervals = TRUE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('c4_curvature', 'c4_slope', 'rL', 'Vmax', 'AIC')]),
        c(0.697636, 1.010108, 1.322194, 65.126933, 72.003838),
        tolerance = 1e-5
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('c4_curvature_upper', 'c4_slope_upper', 'rL_upper', 'Vmax_upper')]),
        c(0.8047086, 1.1376861, 2.8630977, 67.7229577),
        tolerance = 1e-5
    )
})
