# Get test curves to use
source('one_curve_c4_aci.R')

# Choose test tolerance
TOLERANCE <- 1e-4

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
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('c4_curvature', 'c4_slope', 'rL', 'Vmax', 'AIC')]),
        c(0.697636, 1.010108, 1.322194, 65.126933, 72.003838),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('c4_curvature_upper', 'c4_slope_upper', 'rL_upper', 'Vmax_upper')]),
        c(0.8047086, 1.1376861, 2.8630977, 67.7229577),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 4, 9)
    )
})

test_that('removing and excluding points produce the same fit results', {
    pts_to_remove <- c(3, 5, 13)

    one_curve_remove <- remove_points(
        one_curve,
        list(seq_num = pts_to_remove),
        method = 'remove'
    )

    one_curve_exclude <- remove_points(
        one_curve,
        list(seq_num = pts_to_remove),
        method = 'exclude'
    )

    expect_equal(nrow(one_curve_remove), 10)
    expect_equal(nrow(one_curve_exclude), 13)

    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_remove <- fit_c4_aci_hyperbola(
        one_curve_remove,
        OPTIM_FUN = optimizer_nmkb(1e-7)
    )

    set.seed(1234)

    fit_res_exclude <- fit_c4_aci_hyperbola(
        one_curve_exclude,
        OPTIM_FUN = optimizer_nmkb(1e-7)
    )

    # Check that results haven't changed
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('c4_curvature', 'c4_slope', 'rL', 'Vmax', 'AIC')]),
        c(0.657, 1.205, 3.004, 67.773, 51.612),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('npts', 'nparam', 'dof')]),
        c(10, 4, 6)
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('RSS', 'RMSE')]),
        c(45.874, 2.142),
        tolerance = TOLERANCE
    )

    # Check that remove/exclude results are the same
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('c4_curvature', 'c4_slope', 'rL', 'Vmax', 'AIC')]),
        as.numeric(fit_res_exclude$parameters[1, c('c4_curvature', 'c4_slope', 'rL', 'Vmax', 'AIC')]),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('npts', 'nparam', 'dof')]),
        as.numeric(fit_res_exclude$parameters[1, c('npts', 'nparam', 'dof')])
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('RSS', 'RMSE')]),
        as.numeric(fit_res_exclude$parameters[1, c('RSS', 'RMSE')]),
        tolerance = TOLERANCE
    )
})

