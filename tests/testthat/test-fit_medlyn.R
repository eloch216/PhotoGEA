# Get test curves to use
source('one_curve_ball_berry.R')

# Choose test tolerance
TOLERANCE <- 1e-4

test_that('fit results have not changed', {
    fit_res <- fit_medlyn(one_curve)

    expect_equal(
        as.numeric(fit_res$parameters[1, c('medlyn_g0', 'medlyn_g0_err', 'medlyn_g1', 'medlyn_g1_err')]),
        c(0.1035, 0.0338, 0.9325, 0.4318),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('RSS', 'RMSE')]),
        c(0.006109658, 0.029543329),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(7, 2, 5)
    )
})

test_that('removing and excluding points produce the same fit results', {
    pts_to_remove <- c(3, 5)

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

    expect_equal(nrow(one_curve_remove), 5)
    expect_equal(nrow(one_curve_exclude), 7)

    fit_res_remove <- fit_medlyn(one_curve_remove)
    fit_res_exclude <- fit_medlyn(one_curve_exclude)

    # Check that results haven't changed
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('medlyn_g0', 'medlyn_g0_err', 'medlyn_g1', 'medlyn_g1_err')]),
        c(0.1066, 0.0444, 0.9538, 0.5670),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('npts', 'nparam', 'dof')]),
        c(5, 2, 3)
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('RSS', 'RMSE')]),
        c(0.005727538, 0.033845350),
        tolerance = TOLERANCE
    )

    # Check that remove/exclude results are the same
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('medlyn_g0', 'medlyn_g0_err', 'medlyn_g1', 'medlyn_g1_err')]),
        as.numeric(fit_res_exclude$parameters[1, c('medlyn_g0', 'medlyn_g0_err', 'medlyn_g1', 'medlyn_g1_err')]),
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
