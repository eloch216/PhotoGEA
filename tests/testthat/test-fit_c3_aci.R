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
        fit_c3_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            optim_fun = optimizer_nmkb(1e-7),
            hard_constraints = 2,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 2
        )
    )

    expect_equal(unique(fit_res_bad$fits[, 'c3_assimilation_msg']), 'Cc must be >= 0. RL_at_25 must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c3_assimilation_msg'], 'Cc must be >= 0. RL_at_25 must be >= 0')
    expect_true(all(is.na(fit_res_bad$fits[, c('A_fit', 'Ac', 'Aj', 'Ap')])))
    expect_true(all(is.na(fit_res_bad$fits_interpolated[, c('An', 'Ac', 'Aj', 'Ap')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'AIC')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_at_25_upper')])))
})

test_that('Parameter reliability settings are checked', {
    expect_error(
        fit_c3_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            optim_fun = optimizer_nmkb(1e-7),
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
            optim_fun = optimizer_nmkb(1e-7),
            hard_constraints = 0,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 2
        )
    )

    expect_equal(unique(fit_res$fits[, 'c3_assimilation_msg']), '')
    expect_equal(fit_res$parameters[, 'c3_assimilation_msg'], '')
    expect_true(all(!is.na(fit_res$fits[, c('A_fit')])))
})

test_that('Gamma_star can be passed as a column', {
    one_curve_with_gstar <- set_variable(
        one_curve,
        'Gamma_star_at_25',
        'micromol mol^(-1)',
        value = 38.6
    )

    expect_silent(
        fit_c3_aci(
            one_curve_with_gstar,
            fit_options = list(Gamma_star_at_25 = 'column'),
            optim_fun = optimizer_nmkb(1e-7),
            calculate_confidence_intervals = FALSE
        )
    )
})

test_that('Bad optional arguments are caught', {
    expect_error(
        fit_c3_aci(
            one_curve,
            bad_arg_1 = TRUE,
            bad_arg_2 = 45
        ),
        'The following optional arguments are not supported: bad_arg_1, bad_arg_2'
    )
})

test_that('fit results have not changed (no alpha)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 0, alpha_g = 0, alpha_s = 0, alpha_t = 0),
        optim_fun = optimizer_nmkb(1e-7),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
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
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'AIC', 'TleafCnd_avg', 'Jmax_at_25')]),
        c(145.3336224, 232.8361365, 0.3557059, NA, 61.1303101, 30.1448308, 233.8417400),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_at_25_upper', 'Jmax_at_25_upper')]),
        c(152.831071, 238.947894, 1.034651, Inf, 240.012420),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('operating_Ci', 'operating_Cc', 'operating_An', 'operating_An_model')]),
        c(294.70316, 294.70316, 37.51608, 37.85419),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 4, 9)
    )

    expect_equal(
        as.character(fit_res$fits[, 'limiting_process']),
        c('Ac', 'Ac', 'Ac', 'Ac', 'Ac', 'Ac', 'Ac', 'Ac', 'Ac', 'Aj', 'Aj', 'Aj', 'Aj')
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(9, 4, 0))

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 0)
    )

    expect_equal(
        fit_res$parameters[1, 'c3_optional_arguments'],
        ''
    )
})

test_that('fit results have not changed (alpha_old)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 'fit', alpha_g = 0, alpha_s = 0, alpha_t = 0),
        optim_fun = optimizer_deoptim(100),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 0
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
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'AIC')]),
        c(145.33294, 232.83111, 0.35509, 38.1683, 63.13031),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_at_25_upper')]),
        c(152.8274, 238.9449, 1.0343, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 5, 8)
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(9, 4, 0))

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 0)
    )

    expect_equal(
        fit_res$parameters[1, 'c3_optional_arguments'],
        ''
    )
})

test_that('fit results have not changed (alpha_g and alpha_s)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 0, alpha_g = 'fit', alpha_s = 'fit', alpha_t = 0),
        optim_fun = optimizer_deoptim(100),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 1
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
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'AIC')]),
        c(160.572, 254.807, 1.084, 20.0224, 58.184),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_at_25_upper')]),
        c(166.1669, 262.4657, 1.6696, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 6, 7)
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(8, 4, 1))

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 1)
    )

    expect_equal(
        fit_res$parameters[1, 'c3_optional_arguments'],
        ''
    )
})

test_that('fit results have not changed (alpha_g, alpha_s, and alpha_t)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 0, alpha_g = 'fit', alpha_s = 'fit', alpha_t = 'fit'),
        optim_fun = optimizer_deoptim(100),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 1
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
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'AIC')]),
        c(159.271704, 253.817406, 1.559827, NA, 61.614932),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_at_25_upper')]),
        c(167.240166, 261.007015, 1.996194, Inf),
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

    expect_equal(
        fit_res$parameters[1, 'c3_optional_arguments'],
        ''
    )
})

test_that('fit results have not changed (gmc with temperature dependence)', {
    # Redo the temperature calculations
    one_curve_t <- calculate_temperature_response(one_curve, c3_temperature_param_sharkey)

    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve_t,
        Ca_atmospheric = 420,
        fit_options = list(gmc_at_25 = 'fit', Gamma_star_at_25 = 36.94438, Kc_at_25 = 269.3391, Ko_at_25 = 163.7146),
        optim_fun = optimizer_deoptim(100),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 0
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
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'alpha_old', 'gmc_at_25', 'AIC')]),
        c(134.329626, 236.395612, 2.227642, 45.224713, 0.292880, 9.720798, 67.454335),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_at_25_upper', 'alpha_old_upper', 'gmc_at_25_upper')]),
        c(141.4003596, 242.9683848, 2.9652952, Inf, 0.9999473, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_tl_avg', 'J_tl_avg', 'RL_tl_avg', 'Tp_tl_avg', 'gmc_tl_avg')]),
        c(210.40741, 319.61101, 3.06379, 58.43600, 13.58767),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_tl_avg_lower', 'J_tl_avg_lower', 'RL_tl_avg_lower', 'Tp_tl_avg_lower', 'gmc_tl_avg_lower')]),
        c(199.507529, 310.836478, 2.024810, 21.909000, 2.266986),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 6, 7)
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(9, 4, 0))

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 0)
    )

    expect_equal(
        fit_res$parameters[1, 'c3_optional_arguments'],
        ''
    )
})

test_that('fit results have not changed (Kc)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Kc_at_25 = 'fit'),
        optim_fun = optimizer_deoptim(100)
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
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'Kc_at_25', 'AIC')]),
        c(68.0228361, 239.3879473, 0.6131452, NA, 133.3231031, 53.7896820),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_at_25_upper', 'Kc_at_25_upper')]),
        c(69.468809, 245.162135, 1.130932, Inf, 139.323207),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 6, 7)
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(11, 2, 0))

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 0)
    )

    expect_equal(
        fit_res$parameters[1, 'c3_optional_arguments'],
        ''
    )
})

test_that('fit results have not changed (Ko)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Ko_at_25 = 'fit'),
        optim_fun = optimizer_deoptim(100)
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
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'Ko_at_25', 'AIC')]),
        c(116.0374754, 233.7308314, 0.5750712, NA, 970.4221562, 62.6153505),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_at_25_upper', 'Ko_at_25_upper')]),
        c(121.296141, 239.470000, 1.167091, Inf, 1712.602637),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 6, 7)
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(9, 4, 0))

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_trust', 'J_trust', 'Tp_trust')]),
        c(2, 2, 0)
    )

    expect_equal(
        fit_res$parameters[1, 'c3_optional_arguments'],
        ''
    )
})


test_that('fit results have not changed (pseudo-FvCB)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(
        one_curve,
        Ca_atmospheric = 420,
        optim_fun = optimizer_nmkb(1e-7),
        use_min_A = TRUE
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
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'AIC')]),
        c(138.37, 226.84, -0.75, NA, 68.62),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 5, 8)
    )

    lim_info <-
        as.numeric(fit_res$parameters[1, c('n_Ac_limiting', 'n_Aj_limiting', 'n_Ap_limiting')])

    expect_equal(sum(lim_info), nrow(one_curve))

    expect_equal(lim_info, c(7, 6, 0))

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'RL_at_25_upper', 'Tp_at_25_upper')]),
        c(147.65, 235.43, 0.09, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        fit_res$parameters[1, 'c3_optional_arguments'],
        'use_min_A = TRUE'
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

    fit_res_remove <- fit_c3_aci(
        one_curve_remove,
        Ca_atmospheric = 420,
        optim_fun = optimizer_nmkb(1e-7)
    )

    set.seed(1234)

    fit_res_exclude <- fit_c3_aci(
        one_curve_exclude,
        Ca_atmospheric = 420,
        optim_fun = optimizer_nmkb(1e-7)
    )

    # Check that results haven't changed
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'AIC')]),
        c(145.1034379, 234.0332835, 0.4648526, NA, 52.7739136),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('npts', 'nparam', 'dof')]),
        c(10, 5, 5)
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('RSS', 'RMSE')]),
        c(34.539355, 1.858477),
        tolerance = TOLERANCE
    )

    # Check that remove/exclude results are the same
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'AIC')]),
        as.numeric(fit_res_exclude$parameters[1, c('Vcmax_at_25', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'AIC')]),
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
