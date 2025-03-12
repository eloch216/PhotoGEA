test_that('C3 limiting processes are identified', {
    example_curve <- exdf(
      data.frame(
        A_fit = c(1.0, 2.0, 3.0, 4.0, 4.0),
        Ac    = c(1.0, 2.0, 5.0, 8.0, 9.0),
        Aj    = c(2.0, 2.5, 3.0, 4.0, 8.0),
        Ap    = c(NA,  NA,  4.0, 4.0, 4.0)
      ),
      units = data.frame(
        A_fit = 'micromol m^(-2) s^(-1)',
        Ac    = 'micromol m^(-2) s^(-1)',
        Aj    = 'micromol m^(-2) s^(-1)',
        Ap    = 'micromol m^(-2) s^(-1)',
        stringsAsFactors = FALSE
      )
    )

    example_curve <- expect_silent(
        identify_c3_limiting_processes(example_curve)
    )

    expect_equal(
        as.logical(example_curve[, 'Ac_limiting']),
        c(TRUE, TRUE, FALSE, FALSE, FALSE)
    )

    expect_equal(
        as.logical(example_curve[, 'Aj_limiting']),
        c(FALSE, FALSE, TRUE, TRUE, FALSE)
    )

    expect_equal(
        as.logical(example_curve[, 'Ap_limiting']),
        c(FALSE, FALSE, FALSE, TRUE, TRUE)
    )

    expect_equal(
        as.character(example_curve[, 'limiting_process']),
        c('Ac', 'Ac', 'Aj', 'co-limited (Aj and Ap)', 'Ap')
    )
})
