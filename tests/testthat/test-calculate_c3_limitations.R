example_exdf <- exdf(
    data.frame(
        Cc             = seq(10, 400, length.out = 3),
        Gamma_star     = 40,
        gmc            = 0.5,
        gsc            = 0.8,
        Kc             = 270,
        Ko             = 160,
        total_pressure = 1,
        Vcmax_tl       = 150,
        J_tl           = 200
    ),
    units = data.frame(
        Cc             = 'micromol mol^(-1)',
        Gamma_star     = 'micromol mol^(-1)',
        gmc            = 'mol m^(-2) s^(-1) bar^(-1)',
        gsc            = 'mol m^(-2) s^(-1)',
        Kc             = 'micromol mol^(-1)',
        Ko             = 'mmol mol^(-1)',
        total_pressure = 'bar',
        Vcmax_tl       = 'micromol m^(-2) s^(-1)',
        J_tl           = 'micromol m^(-2) s^(-1)'
    )
)

test_that('j limitations are not calculated by default', {
    limit_res <- calculate_c3_limitations(example_exdf)
    expect_true(is.exdf(limit_res))
    expect_true('dAdC_rubisco' %in% colnames(limit_res))
    expect_true('ls_rubisco' %in% colnames(limit_res))
    expect_true('lm_rubisco' %in% colnames(limit_res))
    expect_true('lb_rubisco' %in% colnames(limit_res))
    expect_false('dAdC_j' %in% colnames(limit_res))
    expect_false('ls_j' %in% colnames(limit_res))
    expect_false('lm_j' %in% colnames(limit_res))
    expect_false('lb_j' %in% colnames(limit_res))
})

test_that('limitations add to 1', {
    limit_res <- calculate_c3_limitations(example_exdf, j_column_name = 'J_tl')

    limit_res[, 'l_rubisco'] <- limit_res[, 'ls_rubisco'] + limit_res[, 'lm_rubisco'] + limit_res[, 'lb_rubisco']
    limit_res[, 'l_j'] <- limit_res[, 'ls_j'] + limit_res[, 'lm_j'] + limit_res[, 'lb_j']

    expect_equal(limit_res[, 'l_rubisco'], c(1, 1, 1))
    expect_equal(limit_res[, 'l_j'], c(1, 1, 1))
})
