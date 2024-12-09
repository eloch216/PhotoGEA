test_that('missing data frame columns are detected', {
    expect_error(
        check_required_variables(
            data.frame(A = seq(1, 3)),
            'B'
        ),
        'The following required columns are not present: B'
    )
})

test_that('present data frame columns are detected', {
    expect_silent(
        check_required_variables(
            data.frame(A = seq(1, 3)),
            'A'
        )
    )
})

test_that('NA data frame columns are detected', {
    expect_error(
        check_required_variables(
            data.frame(A = c(NA, NA, NA)),
            'A'
        ),
        'The following required columns are all NA: A'
    )

    expect_silent(
        check_required_variables(
            data.frame(A = c(NA, NA, NA)),
            'A',
            check_NA = FALSE
        )
    )
})

test_that('missing exdf columns are detected', {
    expect_error(
        check_required_variables(
            exdf(data.frame(A = seq(1, 3), B = seq(2, 4))),
            list('C' = NA)
        ),
        'The following required columns are not present: C'
    )
})

test_that('missing exdf units are detected', {
    expect_error(
        check_required_variables(
            exdf(data.frame(A = seq(1, 3), B = seq(2, 4))),
            list('B' = 'B units')
        ),
        'The `B` column must have units of `B units`, but its units are `NA`'
    )
})

test_that('incorrect exdf units are detected', {
    expect_error(
        check_required_variables(
            exdf(data.frame(A = seq(1, 3), B = seq(2, 4)), units = data.frame(A = 'B units', B = 'B units', stringsAsFactors = FALSE)),
            list('A' = 'A units')
        ),
        'The `A` column must have units of `A units`, but its units are `B units`'
    )
})

test_that('present exdf columns and units are detected', {
    expect_silent(
        check_required_variables(
            exdf(data.frame(A = seq(1, 3), B = seq(2, 4)), units = data.frame(A = 'A units', B = 'B units', stringsAsFactors = FALSE)),
            list('A' = 'A units', B = NA)
        )
    )
})

test_that('NA exdf columns are detected', {
    expect_error(
        check_required_variables(
            exdf(data.frame(A = c(NA, NA, NA), B = c(NA, NA, NA)), units = data.frame(A = 'A units', B = 'B units', stringsAsFactors = FALSE)),
            list('A' = 'A units', B = NA)
        ),
        'The following required columns are all NA: A, B'
    )

    expect_silent(
        check_required_variables(
            exdf(data.frame(A = c(NA, NA, NA), B = c(NA, NA, NA)), units = data.frame(A = 'A units', B = 'B units', stringsAsFactors = FALSE)),
            list('A' = 'A units', B = NA),
            check_NA = FALSE
        )
    )
})
