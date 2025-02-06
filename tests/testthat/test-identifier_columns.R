test_that('single identifiers can be found', {
    # Create a simple exdf object
    simple_exdf <- exdf(
        data.frame(A = c(3, 2, 7, 9), species = c('a', 'a', 'a', 'a')),
        data.frame(A = 'm', species = 'none'),
        data.frame(A = 'Cat1', species = 'meas')
    )

    # Find its (exdf) identifier columns and make sure they are correct
    exdf_id_col <- expect_silent(identifier_columns(simple_exdf))

    expect_equal(
        exdf_id_col,
        exdf(
            data.frame(species = 'a'),
            data.frame(species = 'none'),
            data.frame(species = 'meas')
        )
    )

    # Find its (data frame) identifier columns and make sure they are correct
    df_id_col <- expect_silent(identifier_columns(simple_exdf$main_data))

    expect_equal(
        df_id_col,
        data.frame(species = 'a')
    )
})
