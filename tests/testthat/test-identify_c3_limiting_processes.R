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

test_that('Assimilation rates and carboxylation rates are consistent', {
    inputs <- exdf(data.frame(
        Cc = c(1, 60),
        oxygen = 21,
        Tleaf = 30,
        total_pressure = 1
    ))

    inputs <- document_variables(
        inputs,
        c('', 'Cc',             'micromol mol^(-1)'),
        c('', 'oxygen',         'percent'),
        c('', 'Tleaf',          'degrees C'),
        c('', 'total_pressure', 'bar')
    )

    inputs <- calculate_temperature_response(inputs, c3_temperature_param_sharkey, 'Tleaf')

    fvcb_res <- identify_c3_limiting_processes(
        calculate_c3_assimilation(inputs, 0, 0, 0, 0, '', 150, '', '', 1, 12, 120),
        a_column_name = 'An'
    )

    min_a_res <- identify_c3_limiting_processes(
        calculate_c3_assimilation(inputs, 0, 0, 0, 0, '', 150, '', '', 1, 12, 120, use_min_A = TRUE),
        a_column_name = 'An'
    )

    # In the FvCB result, the first point is Rubisco-limited, and so is the
    # second point
    expect_equal(fvcb_res[1, 'limiting_process'], 'Ac')
    expect_equal(fvcb_res[1, 'An'], fvcb_res[1, 'Ac'])
    expect_equal(fvcb_res[1, 'Vc'], fvcb_res[1, 'Wc'])
    expect_true(fvcb_res[1, 'Ac'] != fvcb_res[1, 'Aj'])
    expect_true(fvcb_res[1, 'Wc'] != fvcb_res[1, 'Wj'])

    expect_equal(fvcb_res[2, 'limiting_process'], 'Ac')
    expect_equal(fvcb_res[2, 'An'], fvcb_res[2, 'Ac'])
    expect_equal(fvcb_res[2, 'Vc'], fvcb_res[2, 'Wc'])
    expect_true(fvcb_res[2, 'Ac'] != fvcb_res[2, 'Aj'])
    expect_true(fvcb_res[2, 'Wc'] != fvcb_res[2, 'Wj'])

    # In the min-A result, the first point is RuBP-regeneration-limited, and the
    # second point is Rubisco-limited
    expect_equal(min_a_res[1, 'limiting_process'], 'Aj')
    expect_equal(min_a_res[1, 'An'], min_a_res[1, 'Aj'])
    expect_equal(min_a_res[1, 'Vc'], min_a_res[1, 'Wj'])
    expect_true(min_a_res[1, 'Ac'] != min_a_res[1, 'Aj'])
    expect_true(min_a_res[1, 'Wc'] != min_a_res[1, 'Wj'])

    expect_equal(min_a_res[2, 'limiting_process'], 'Ac')
    expect_equal(min_a_res[2, 'An'], min_a_res[2, 'Ac'])
    expect_equal(min_a_res[2, 'Vc'], min_a_res[2, 'Wc'])
    expect_true(min_a_res[2, 'Ac'] != min_a_res[2, 'Aj'])
    expect_true(min_a_res[2, 'Wc'] != min_a_res[2, 'Wj'])
})
