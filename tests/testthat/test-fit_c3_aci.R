# Get test curves to use
source('one_curve_c3_aci.R')

# Specify an infinite mesophyll conductance (so `Cc` = `Ci`)
one_curve <- set_variable(
  one_curve,
  'gmc', 'mol m^(-2) s^(-1) bar^(-1)', value = Inf
)

one_curve_bad <- set_variable(
  one_curve_bad,
  'gmc', 'mol m^(-2) s^(-1) bar^(-1)', value = Inf
)

# Calculate Cc
one_curve <- apply_gm(one_curve)
one_curve_bad <- apply_gm(one_curve_bad)

# Choose test tolerance
TOLERANCE <- 1e-4

test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- expect_silent(
        fit_c3_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            OPTIM_FUN = optimizer_nmkb(1e-7),
            hard_constraints = 2,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 2
        )
    )

    expect_equal(unique(fit_res_bad$fits[, 'c3_assimilation_msg']), 'Cc must be >= 0. RL_at_25 must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c3_assimilation_msg'], 'Cc must be >= 0. RL_at_25 must be >= 0')
    expect_true(all(is.na(fit_res_bad$fits[, c('A_fit', 'Ac', 'Aj', 'Ap')])))
    expect_true(all(is.na(fit_res_bad$fits_interpolated[, c('An', 'Ac', 'Aj', 'Ap')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp', 'AIC')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_upper')])))
})

test_that('Parameter reliability settings are checked', {
    expect_error(
        fit_c3_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            OPTIM_FUN = optimizer_nmkb(1e-7),
            hard_constraints = 2,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 10
        ),
        '`remove_unreliable_param` must be 0, 1, or 2'
    )
})

test_that('Cc limits can be bypassed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c3_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            OPTIM_FUN = optimizer_nmkb(1e-7),
            hard_constraints = 0,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 2
        )
    )

    expect_equal(unique(fit_res$fits[, 'c3_assimilation_msg']), '')
    expect_equal(fit_res$parameters[, 'c3_assimilation_msg'], '')
    expect_true(all(!is.na(fit_res$fits[, c('A_fit')])))
})

test_that('fit results have not changed (no alpha)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 0, alpha_g = 0, alpha_s = 0),
        OPTIM_FUN = optimizer_nmkb(1e-7),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp', 'AIC')]),
        c(145.3336224, 232.8361365, 0.3557059, NA, 59.1303101),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_upper')]),
        c(152.831071, 238.947894, 1.034651, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 0)
    )
})

test_that('fit results have not changed (alpha_old)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 'fit', alpha_g = 0, alpha_s = 0),
        OPTIM_FUN = optimizer_deoptim(100),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 0
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp', 'AIC')]),
        c(145.33294, 232.83111, 0.35509, 38.1683, 61.13031),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_upper')]),
        c(152.8274, 238.9449, 1.0343, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 0)
    )
})

test_that('fit results have not changed (alpha_g and alpha_s)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 0, alpha_g = 'fit', alpha_s = 'fit'),
        OPTIM_FUN = optimizer_deoptim(100),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 1
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp', 'AIC')]),
        c(160.572, 254.807, 1.084, 20.0224, 56.184),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_upper')]),
        c(166.1669, 262.4657, 1.6696, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 1)
    )
})
