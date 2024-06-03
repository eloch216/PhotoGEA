test_curve <- exdf(
    data.frame(
        A = c(0, 1),
        Ca = c(0, 10),
        Cc = c(0, 1e2),
        Ci = c(0, 1e3),
        PCm = c(0, 1e4)
    ),
    units = data.frame(
        A = 'micromol m^(-2) s^(-1)',
        Ca = 'micromol mol^(-1)',
        Cc = 'micromol mol^(-1)',
        Ci = 'micromol mol^(-1)',
        PCm = 'microbar',
        stringsAsFactors = FALSE
    )
)

test_that('default c3 operating point estimation works', {
    c3_res <- estimate_operating_point(test_curve, 4.2)
    expect_true(is.exdf(c3_res))
    expect_equal(c3_res[, 'operating_An'], 0.42)
    expect_equal(c3_res[, 'operating_Ci'], 420)
    expect_equal(c3_res[, 'operating_Cc'], 42)
    expect_false('operating_PCm' %in% colnames(c3_res))
})

test_that('default c4 operating point estimation works', {
    c4_res <- estimate_operating_point(test_curve, 4.2, type = 'c4')
    expect_true(is.exdf(c4_res))
    expect_equal(c4_res[, 'operating_An'], 0.42)
    expect_equal(c4_res[, 'operating_Ci'], 420)
    expect_equal(c4_res[, 'operating_PCm'], 4200)
    expect_false('operating_Cc' %in% colnames(c4_res))
})

test_that('list c3 operating point estimation works', {
    c3_res <- estimate_operating_point(test_curve, 4.2, return_list = TRUE)
    expect_true(is.list(c3_res))
    expect_true(is.exdf(c3_res$operating_exdf))
    expect_equal(c3_res$operating_An, 0.42)
    expect_equal(c3_res$operating_Ci, 420)
    expect_equal(c3_res$operating_Cc, 42)
    expect_false('operating_PCm' %in% names(c3_res))
})

test_that('list c4 operating point estimation works', {
    c4_res <- estimate_operating_point(test_curve, 4.2, type = 'c4', return_list = TRUE)
    expect_true(is.list(c4_res))
    expect_true(is.exdf(c4_res$operating_exdf))
    expect_equal(c4_res$operating_An, 0.42)
    expect_equal(c4_res$operating_Ci, 420)
    expect_equal(c4_res$operating_PCm, 4200)
    expect_false('operating_Cc' %in% names(c4_res))
})

test_that('extrapolation is not allowed', {
    expect_warning(
        bad_res <- estimate_operating_point(test_curve, 100)
    )
    expect_true(is.na(bad_res[, 'operating_Ci']))
    expect_true(is.na(bad_res[, 'operating_An']))
})
