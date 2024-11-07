# Get test curves to use
source('one_curve_c4_aci.R')

# Calculate PCm
one_curve <- apply_gm(one_curve, 'C4')
one_curve_bad <- apply_gm(one_curve_bad, 'C4')

# Choose test tolerance
TOLERANCE <- 1e-4

test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- expect_silent(
        fit_c4_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            OPTIM_FUN = optimizer_nmkb(1e-7),
            hard_constraints = 2,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 2
        )
    )

    expect_equal(unique(fit_res_bad$fits[, 'c4_assimilation_msg']), 'PCm must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c4_assimilation_msg'], 'PCm must be >= 0')
    expect_true(all(is.na(fit_res_bad$fits[, c('A_fit', 'Apr', 'Apc', 'Ar')])))
    expect_true(all(is.na(fit_res_bad$fits_interpolated[, c('An', 'Apr', 'Apc', 'Ar')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25_upper', 'Vpmax_at_25_upper', 'RL_at_25_upper')])))
})

test_that('PCm limits can be bypassed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c4_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            OPTIM_FUN = optimizer_nmkb(1e-7),
            hard_constraints = 0,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 2
        )
    )

    expect_equal(unique(fit_res$fits[, 'c4_assimilation_msg']), '')
    expect_equal(fit_res$parameters[, 'c4_assimilation_msg'], '')
    expect_true(all(!is.na(fit_res$fits[, c('A_fit')])))
})

test_that('fit results have not changed (Vcmax)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Vcmax_at_25 = 'fit', Vpr = 1000, Jmax_at_25 = 1000),
        OPTIM_FUN = optimizer_nmkb(1e-7),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC')]),
        c(3.630116e+01, 1.804791e+02, 1.069116e-08, 8.026640e+01),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'Vpmax_at_25_upper', 'RL_at_25_upper')]),
        c(38.434695, 214.046523, 1.568026),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 3, 10)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpmax_trust', 'Vcmax_trust', 'Vpr_trust', 'Jmax_trust')]),
        c(2, 2, 0, 0)
    )
})

test_that('fit results have not changed (Vpr)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Vcmax_at_25 = 1000, Vpr = 'fit', Jmax_at_25 = 1000),
        OPTIM_FUN = optimizer_nmkb(1e-7),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpr', 'Vpmax_at_25', 'RL_at_25', 'AIC')]),
        c(58.1571, 133.8038, 0.0000, 86.3427),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpr_upper', 'Vpmax_at_25_upper', 'RL_at_25_upper')]),
        c(62.43, 156.94, 2.76),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 3, 10)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpmax_trust', 'Vcmax_trust', 'Vpr_trust', 'Jmax_trust')]),
        c(2, 1, 2, 0)
    )
})

test_that('fit results have not changed (Jmax)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Vcmax_at_25 = 1000, Vpr = 1000, Jmax_at_25 = 'fit'),
        OPTIM_FUN = optimizer_nmkb(1e-7),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Jmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC')]),
        c(5.215746e+02, 1.338467e+02, 1.475187e-08, 8.675720e+01),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Jmax_at_25_upper', 'Vpmax_at_25_upper', 'RL_at_25_upper')]),
        c(573.15632, 157.30750, 2.24672),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 3, 10)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpmax_trust', 'Vcmax_trust', 'Vpr_trust', 'Jmax_trust')]),
        c(2, 1, 0, 2)
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

    fit_res_remove <- fit_c4_aci(
        one_curve_remove,
        Ca_atmospheric = 420,
        OPTIM_FUN = optimizer_nmkb(1e-7)
    )

    set.seed(1234)

    fit_res_exclude <- fit_c4_aci(
        one_curve_exclude,
        Ca_atmospheric = 420,
        OPTIM_FUN = optimizer_nmkb(1e-7)
    )

    # Check that results haven't changed
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC')]),
        c(37.01, 211.08, 0.70, 57.67),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('npts', 'nparam', 'dof')]),
        c(10, 3, 7)
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('RSS', 'RMSE')]),
        c(102.677, 3.204),
        tolerance = TOLERANCE
    )

    # Check that remove/exclude results are the same
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC')]),
        as.numeric(fit_res_exclude$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC')]),
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

