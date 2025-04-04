test_that('Unknown quantity names are detected', {
    expect_error(
        PhotoGEA:::unit_dictionary('bad_quantity'),
        'Units were requested for a quantity named `bad_quantity`, but it is not included in the unit dictionary'
    )
})

test_that('Units are returned properly', {
    expect_equal(
        PhotoGEA:::unit_dictionary('Ci'),
        'micromol mol^(-1)'
    )
})

test_that('Normalized units are formatted correctly', {
    expect_equal(
        PhotoGEA:::normalized_units('test_name'),
        'normalized to test_name at 25 degrees C'
    )
})
