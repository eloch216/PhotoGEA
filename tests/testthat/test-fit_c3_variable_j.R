# Get test curves to use
source('one_curve_c3_aci.R')

# Load helping function
source('get_duplicated_colnames.R')

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
            optim_fun = optimizer_deoptim(200),
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
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'tau', 'AIC')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_at_25_upper', 'tau_upper')])))
})

test_that('Ci and Cc limits can be bypassed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c3_variable_j(
            one_curve_bad,
            Ca_atmospheric = 420,
            optim_fun = optimizer_deoptim(200),
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

test_that('Gamma_star can be passed via fit_options', {
    one_curve_no_gstar <- one_curve
    one_curve_no_gstar[, 'Gamma_star'] <- NULL

    expect_silent(
        fit_c3_variable_j(
            one_curve_no_gstar,
            fit_options = list(Gamma_star = 38.6),
            optim_fun = optimizer_deoptim(200),
            calculate_confidence_intervals = FALSE
        )
    )
})

test_that('fit results have not changed (no alpha)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_variable_j(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 0, alpha_g = 0, alpha_s = 0),
        optim_fun = optimizer_deoptim(200),
        require_positive_gmc = 'all',
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2,
        check_j = FALSE
    )

    fit_res$parameters <- calculate_temperature_response(
        fit_res$parameters,
        jmax_temperature_param_bernacchi,
        'TleafCnd_avg'
    )

    fit_res$parameters <- calculate_jmax(
        fit_res$parameters,
        0.6895,
        0.97875
    )

    expect_equal(
        get_duplicated_colnames(fit_res$fits),
        character(0)
    )

    expect_equal(
        get_duplicated_colnames(fit_res$parameters),
        character(0)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'tau', 'Tp_at_25', 'AIC', 'TleafCnd_avg', 'Jmax_at_25')]),
        c(240.718, 254.101, 1.885, 0.405, NA, 40.416, 30.1448308, 255.3210905),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'tau_upper', 'Tp_at_25_upper')]),
        c(247.455, 256.611, 1.892, 0.409, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 5, 8)
    )

    expect_equal(
        as.character(fit_res$fits[, 'limiting_process']),
        c('Ac', 'Ac', 'Ac', 'Ac', 'Ac', 'Ac', 'Ac', 'Ac', 'Aj', 'Aj', 'Aj', 'Aj', 'Aj')
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(8, 5, 0))

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
        optim_fun = optimizer_deoptim(200),
        require_positive_gmc = 'all',
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2,
        check_j = FALSE
    )

    expect_equal(
        get_duplicated_colnames(fit_res$fits),
        character(0)
    )

    expect_equal(
        get_duplicated_colnames(fit_res$parameters),
        character(0)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'tau', 'Tp_at_25', 'AIC')]),
        c(243.821, 256.166, 1.901, 0.409, NA, 42.429),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'tau_upper', 'Tp_at_25_upper')]),
        c(250.388, 258.740, 1.912, 0.412, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 6, 7)
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(8, 5, 0))

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
        optim_fun = optimizer_deoptim(200),
        require_positive_gmc = 'all',
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2,
        check_j = FALSE
    )

    expect_equal(
        get_duplicated_colnames(fit_res$fits),
        character(0)
    )

    expect_equal(
        get_duplicated_colnames(fit_res$parameters),
        character(0)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'tau', 'Tp_at_25', 'AIC')]),
        c(223.833, 264.193, 1.798, 0.422, NA, 45.052),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'tau_upper', 'Tp_at_25_upper')]),
        c(230.9429, 267.1610, 2.4567, 0.4251, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 7, 6)
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(8, 5, 0))

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 0)
    )
})

test_that('fit results have not changed (pseudo-FvCB)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_variable_j(
        one_curve,
        Ca_atmospheric = 420,
        optim_fun = optimizer_deoptim(200),
        use_min_A = TRUE,
        check_j = FALSE
    )

    expect_equal(
        get_duplicated_colnames(fit_res$fits),
        character(0)
    )

    expect_equal(
        get_duplicated_colnames(fit_res$parameters),
        character(0)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'tau', 'Tp_at_25', 'AIC')]),
        c(319.868, 313.808, 2.441, 0.500, NA, 49.966),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 6, 7)
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(7, 6, 0))

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'tau_upper', 'Tp_at_25_upper')]),
        c(331.800, 317.787, 2.445, 0.505, Inf),
        tolerance = TOLERANCE
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

    fit_res_remove <- fit_c3_variable_j(
        one_curve_remove,
        Ca_atmospheric = 420,
        optim_fun = optimizer_deoptim(200)
    )

    set.seed(1234)

    fit_res_exclude <- fit_c3_variable_j(
        one_curve_exclude,
        Ca_atmospheric = 420,
        optim_fun = optimizer_deoptim(200)
    )

    # Check that results haven't changed
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'tau', 'Tp_at_25', 'AIC')]),
        c(268.81, 277.21, 1.97, 0.44, NA, 41.90),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('npts', 'nparam', 'dof')]),
        c(10, 6, 4)
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('RSS', 'RMSE')]),
        c(9.534, 0.976),
        tolerance = TOLERANCE
    )

    # Check that remove/exclude results are the same
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'tau', 'Tp_at_25', 'AIC')]),
        as.numeric(fit_res_exclude$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'tau', 'Tp_at_25', 'AIC')]),
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

