# Define a simple exdf object
test_exdf <- exdf(
    data.frame(species = c('maize', 'maize', 'soybean'), gsw = c(0.2, 0.3, 0.4)),
    units = data.frame(species = '', gsw = 'mol m^(-2) s^(-1)', stringsAsFactors = FALSE)
)

test_that('new columns can be added with set_variable', {
    # Set mesophyll conductance to infinity
    gmc_v1 <- set_variable(
        test_exdf,
        'gmc',
        'mol m^(-2) s^(-1) bar^(-1)',
        value = Inf
    )

    expect_equal(
        gmc_v1[, 'gmc'],
        c(Inf, Inf, Inf)
    )

    expect_equal(gmc_v1$units$gmc, 'mol m^(-2) s^(-1) bar^(-1)')

    expect_equal(gmc_v1$categories$gmc, 'NA')

    # Set mesophyll conductance based on species
    gmc_v2 <- set_variable(
        test_exdf,
        'gmc',
        'mol m^(-2) s^(-1) bar^(-1)',
        'testing',
        id_column = 'species',
        value_table = list(maize = 0.9, soybean = 1.1)
    )

    expect_equal(
        gmc_v2[, 'gmc'],
        c(0.9, 0.9, 1.1)
    )

    expect_equal(gmc_v2$units$gmc, 'mol m^(-2) s^(-1) bar^(-1)')

    expect_equal(gmc_v2$categories$gmc, 'testing')
})

test_that('common mistakes are caught', {
    expect_error(
        set_variable(
            test_exdf,
            'gmc',
            'mol m^(-2) s^(-1) bar^(-1)',
            id_column = 'species',
            value_table = Inf
        ),
        '`value_table` must be a list'
    )

    expect_error(
        set_variable(
            test_exdf,
            'gmc',
            'mol m^(-2) s^(-1) bar^(-1)',
            id_column = 'species',
            value_table = list(0.9, soybean = 1.1)
        ),
        'All elements of `value_table` must have names'
    )

    expect_error(
        set_variable(
            test_exdf,
            'gmc',
            'mol m^(-2) s^(-1) bar^(-1)',
            id_column = 'species',
            value_table = list(0.9, 1.1)
        ),
        'All elements of `value_table` must have names'
    )

    expect_error(
        set_variable(
            test_exdf,
            'gmc',
            'mol m^(-2) s^(-1) bar^(-1)',
            value_table = list(maize = 0.9, soybean = 1.1)
        ),
        'An `id_column` must be supplied when using a `value_table`'
    )

})
