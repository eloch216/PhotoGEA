test_that('missing data frame columns are detected', {
    expect_error(
        check_required_variables(
            data.frame(A = seq(1, 3)),
            'B'
        ),
        'The following columns are undefined: B'
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

test_that('missing exdf columns are detected', {
    expect_error(
        check_required_variables(
            exdf(data.frame(A = seq(1, 3), B = seq(2, 4))),
            list('C' = NA)
        ),
        'The following columns are undefined: C'
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
            exdf(data.frame(A = seq(1, 3), B = seq(2, 4)), units = data.frame(A = 'B units', B = 'B units')),
            list('A' = 'A units')
        ),
        'The `A` column must have units of `A units`, but its units are `B units`'
    )
})

test_that('present exdf columns and units are detected', {
    expect_silent(
        check_required_variables(
            exdf(data.frame(A = seq(1, 3), B = seq(2, 4)), units = data.frame(A = 'A units', B = 'B units')),
            list('A' = 'A units', B = NA)
        )
    )
})
