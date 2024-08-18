# Get test curves to use
source('one_curve_c3_aci.R')

# Choose test tolerance
TOLERANCE <- 1e-4

test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- expect_silent(
        fit_c3_variable_j(
            one_curve_bad,
            Ca_atmospheric = 420,
            OPTIM_FUN = optimizer_deoptim(200),
            hard_constraints = 2,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 2
        )
    )

    expect_equal(unique(fit_res_bad$fits[, 'c3_assimilation_msg']), 'alpha_old must be >= 0 and <= 1')
    expect_equal(unique(fit_res_bad$fits[, 'c3_variable_j_msg']), 'Ci must be >= 0. tau must be >= 0 and <= 1')
    expect_equal(fit_res_bad$parameters[, 'c3_assimilation_msg'], 'alpha_old must be >= 0 and <= 1')
    expect_equal(fit_res_bad$parameters[, 'c3_variable_j_msg'], 'Ci must be >= 0. tau must be >= 0 and <= 1')
    expect_true(all(is.na(fit_res_bad$fits[, c('A_fit', 'Ac', 'Aj', 'Ap', 'gmc', 'Cc')])))
    expect_true(all(is.na(fit_res_bad$fits_interpolated[, c('An', 'Ac', 'Aj', 'Ap', 'gmc', 'Cc')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp', 'tau', 'AIC')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp', 'tau_upper')])))
})

test_that('Ci and Cc limits can be bypassed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c3_variable_j(
            one_curve_bad,
            Ca_atmospheric = 420,
            OPTIM_FUN = optimizer_deoptim(200),
            hard_constraints = 0,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 2
        )
    )

    expect_equal(unique(fit_res$fits[, 'c3_assimilation_msg']), '')
    expect_equal(unique(fit_res$fits[, 'c3_variable_j_msg']), '')
    expect_equal(fit_res$parameters[, 'c3_assimilation_msg'], '')
    expect_equal(fit_res$parameters[, 'c3_variable_j_msg'], '')
    expect_true(all(!is.na(fit_res$fits[, c('A_fit', 'gmc', 'Cc')])))
})

test_that('fit results have not changed (no alpha)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_variable_j(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 0, alpha_g = 0, alpha_s = 0),
        OPTIM_FUN = optimizer_deoptim(200),
        require_positive_gmc = 'all',
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'tau', 'Tp', 'AIC')]),
        c(240.718, 254.101, 1.885, 0.405, NA, 38.416),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'tau_upper', 'Tp_upper')]),
        c(247.455, 256.611, 1.892, 0.409, Inf),
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

    fit_res <- fit_c3_variable_j(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 'fit', alpha_g = 0, alpha_s = 0),
        OPTIM_FUN = optimizer_deoptim(200),
        require_positive_gmc = 'all',
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'tau', 'Tp', 'AIC')]),
        c(243.821, 256.166, 1.901, 0.409, NA, 40.429),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'tau_upper', 'Tp_upper')]),
        c(250.388, 258.740, 1.912, 0.412, Inf),
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

    fit_res <- fit_c3_variable_j(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 0, alpha_g = 'fit', alpha_s = 'fit'),
        OPTIM_FUN = optimizer_deoptim(200),
        require_positive_gmc = 'all',
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'tau', 'Tp', 'AIC')]),
        c(223.833, 264.193, 1.798, 0.422, NA, 43.052),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'tau_upper', 'Tp_upper')]),
        c(230.9429, 267.1610, 2.4567, 0.4251, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 0)
    )
})
